#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pathlib
import re
import urllib.request
from typing import Any

from workflow.core.runtime import (
    artifacts_exist,
    build_docker_run_command,
    complete_artifact,
    ensure_dir,
    ensure_docker_ready,
    hash_strings,
    load_json,
    parse_memory_to_mb,
    plan_worker_pool,
    remove_if_exists,
    run_command,
    run_tasks_in_parallel,
    stage_payload,
    write_json,
)


STAGE_VERSION = 6
UNION_SITE_VERSION = 2
GATK_IMAGE = "broadinstitute/gatk:4.6.1.0"
PLINK2_IMAGE = "ghcr.io/pgscatalog/plink2:2.00a5.10"
OVERLAP_RE = re.compile(r"Score (\S+) fails minimum matching threshold \(([\d.]+)% variants match\)")


class IndexedFasta:
    def __init__(self, fasta_path: pathlib.Path, fai_path: pathlib.Path) -> None:
        self.handle = open(fasta_path, "rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with open(fai_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                contig, length, offset, line_bases, line_width = raw.rstrip("\n").split("\t")[:5]
                self.index[contig] = (int(length), int(offset), int(line_bases), int(line_width))

    def get_base(self, contig: str, pos: int) -> str:
        length, offset, line_bases, line_width = self.index[contig]
        if pos < 1 or pos > length:
            raise ValueError(f"Position out of range for {contig}:{pos}")
        zero = pos - 1
        byte_offset = offset + (zero // line_bases) * line_width + (zero % line_bases)
        self.handle.seek(byte_offset)
        return self.handle.read(1).decode("ascii").upper()


def parse_overlap_diagnostics(log_text: str) -> list[dict[str, str]]:
    matches: list[dict[str, str]] = []
    for line in log_text.splitlines():
        found = OVERLAP_RE.search(line)
        if found:
            matches.append({"score": found.group(1), "match_percent": found.group(2)})
    seen = set()
    deduped = []
    for item in matches:
        key = (item["score"], item["match_percent"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(item)
    return deduped


def safe_sampleset_name(sample_id: str) -> str:
    out = re.sub(r"[^A-Za-z0-9-]+", "-", sample_id).strip("-")
    out = re.sub(r"-{2,}", "-", out)
    if not out or out.lower() == "reference":
        return "sample"
    return out


def plink_sex_code(sample_sex: str | None) -> str | None:
    if not sample_sex:
        return None
    normalized = str(sample_sex).strip().lower()
    mapping = {
        "male": "M",
        "female": "F",
    }
    if normalized not in mapping:
        raise ValueError("analysis.prs.sample_sex must be 'male' or 'female'.")
    return mapping[normalized]


def normalize_contig(raw: str, prefix: str) -> str:
    raw = raw.strip()
    if prefix == "chr" and not raw.startswith("chr"):
        return f"chr{raw}"
    if not prefix and raw.startswith("chr"):
        return raw[3:]
    return raw


def parse_base_tokens(raw: str | None) -> set[str]:
    if not raw:
        return set()
    tokens = re.split(r"[\s,;/|]+", str(raw).upper().strip())
    return {token for token in tokens if token in {"A", "C", "G", "T"}}


def resolve_score_row_alleles(row: dict[str, str], ref: str) -> tuple[dict[str, Any] | None, str | None]:
    effect_alleles = parse_base_tokens(row.get("effect_allele"))
    if len(effect_alleles) != 1:
        return None, "unsupported_effect_allele"

    explicit_other = parse_base_tokens(row.get("other_allele"))
    if explicit_other:
        source = "explicit_other_allele"
        other_alleles = explicit_other
    else:
        inferred_other = parse_base_tokens(row.get("hm_inferOtherAllele"))
        if inferred_other:
            source = "hm_inferOtherAllele"
            other_alleles = inferred_other
        else:
            return None, "missing_other_allele"

    effect = next(iter(effect_alleles))
    all_alleles = set(other_alleles)
    all_alleles.add(effect)
    if ref not in all_alleles:
        return None, "reference_not_in_score_alleles"

    alts = all_alleles - {ref}
    if not alts:
        return None, "no_alt_after_reference_match"

    return {
        "source": source,
        "effect_allele": effect,
        "other_alleles": sorted(other_alleles),
        "alts": sorted(alts),
        "multiallelic": len(alts) > 1,
    }, None


def fetch_scorefile_url(pgs_id: str) -> str:
    with urllib.request.urlopen(f"https://www.pgscatalog.org/rest/score/{pgs_id}", timeout=60) as response:
        payload = json.load(response)
    return payload["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]


def build_union_sites(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    reference_fasta: pathlib.Path,
    reference_fai: pathlib.Path,
) -> dict[str, Any]:
    pgs_ids = context["manifest"]["analysis"]["prs"]["pgs_ids"]
    prefix = context["contig_prefix"]
    score_hash = hash_strings([f"union-v{UNION_SITE_VERSION}", *pgs_ids])[:16]
    targets_root = ensure_dir(stage_root / score_hash / "targets")
    union_manifest_path = targets_root / "union_manifest.json"
    reuse_key = {
        "union_site_version": UNION_SITE_VERSION,
        "pgs_ids": pgs_ids,
        "reference_bundle_id": context["reference_bundle"].get("bundle_id"),
        "reference_fasta": str(reference_fasta),
        "contig_prefix": prefix,
    }
    if union_manifest_path.exists():
        payload = load_json(union_manifest_path)
        artifacts = payload.get("artifacts") or {}
        if payload.get("reuse_key") == reuse_key and artifacts_exist([artifacts.get("targets_vcf"), artifacts.get("targets_vcf_index")]):
            return payload

    score_cache_root = ensure_dir(pathlib.Path(context["runtime"]["tool_cache_root"]) / "pgs_scorefiles" / "GRCh38")
    fasta = IndexedFasta(reference_fasta, reference_fai)
    contig_order = {contig: idx for idx, contig in enumerate(fasta.index.keys())}
    site_map: dict[tuple[str, int], dict[str, object]] = {}
    per_score_counts: dict[str, int] = {}
    accepted_rows_by_source = {
        "explicit_other_allele": 0,
        "hm_inferOtherAllele": 0,
    }
    multiallelic_rows = 0
    skipped = 0
    skipped_reasons: dict[str, int] = {}

    for pgs_id in pgs_ids:
        url = fetch_scorefile_url(pgs_id)
        score_path = score_cache_root / pathlib.Path(url).name
        if not score_path.exists():
            urllib.request.urlretrieve(url, score_path)

        per_score_counts[pgs_id] = 0
        with gzip.open(score_path, "rt", encoding="utf-8", errors="replace") as handle:
            header = None
            for raw in handle:
                if raw.startswith("#"):
                    continue
                header = raw.rstrip("\n").split("\t")
                break
            if not header:
                raise RuntimeError(f"Could not read header from {score_path}")

            reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
            for row in reader:
                chrom_raw = row.get("hm_chr") or row.get("chr_name")
                pos_raw = row.get("hm_pos") or row.get("chr_position")
                if not chrom_raw or not pos_raw:
                    skipped += 1
                    skipped_reasons["missing_position"] = skipped_reasons.get("missing_position", 0) + 1
                    continue

                chrom = normalize_contig(chrom_raw, prefix)
                pos = int(pos_raw)
                try:
                    ref = fasta.get_base(chrom, pos)
                except Exception:
                    skipped += 1
                    skipped_reasons["reference_lookup_failed"] = skipped_reasons.get("reference_lookup_failed", 0) + 1
                    continue

                resolved, reason = resolve_score_row_alleles(row, ref)
                if resolved is None:
                    skipped += 1
                    skipped_reasons[reason] = skipped_reasons.get(reason, 0) + 1
                    continue

                key = (chrom, pos)
                entry = site_map.setdefault(key, {"ref": ref, "alts": set(), "pgs_ids": set()})
                entry["alts"].update(resolved["alts"])
                entry["pgs_ids"].add(pgs_id)
                per_score_counts[pgs_id] += 1
                accepted_rows_by_source[resolved["source"]] = accepted_rows_by_source.get(resolved["source"], 0) + 1
                if resolved["multiallelic"]:
                    multiallelic_rows += 1

    plain_vcf = targets_root / "pgs_union_sites.vcf"
    with open(plain_vcf, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, meta in fasta.index.items():
            out.write(f"##contig=<ID={contig},length={meta[0]}>\n")
        out.write("##INFO=<ID=PGS_IDS,Number=.,Type=String,Description=\"PGS IDs contributing to this union site\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos in sorted(site_map.keys(), key=lambda item: (contig_order[item[0]], item[1])):
            entry = site_map[(chrom, pos)]
            alt = ",".join(sorted(entry["alts"]))
            info = "PGS_IDS=" + ",".join(sorted(entry["pgs_ids"]))
            out.write(f"{chrom}\t{pos}\t.\t{entry['ref']}\t{alt}\t.\tPASS\t{info}\n")

    targets_vcf = targets_root / "pgs_union_sites.vcf.gz"
    run_command(["bcftools", "view", "-Oz", "-o", str(targets_vcf), str(plain_vcf)], check=True)
    run_command(["bcftools", "index", "-t", str(targets_vcf)], check=True)

    payload = {
        "stage": "prs_union_sites",
        "status": "ok",
        "reuse_key": reuse_key,
        "score_hash": score_hash,
        "summary": {
            "union_sites": len(site_map),
            "skipped_rows": skipped,
            "skipped_reasons": skipped_reasons,
            "per_score_counts": per_score_counts,
            "accepted_rows_by_source": accepted_rows_by_source,
            "accepted_multiallelic_rows": multiallelic_rows,
        },
        "artifacts": {
            "targets_vcf": str(targets_vcf),
            "targets_vcf_index": str(targets_vcf) + ".tbi",
            "plain_vcf": str(plain_vcf),
        },
    }
    write_json(union_manifest_path, payload)
    return payload


def select_prs_shards(context: dict[str, Any], targets_vcf: pathlib.Path) -> list[dict[str, Any]]:
    completed = run_command(["tabix", "-l", str(targets_vcf)], check=True)
    target_contigs = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    context_shards = context["shards"]["prs"]
    selected = []
    seen = set(target_contigs)
    for shard in context_shards:
        shard_intervals = [interval for interval in shard["intervals"] if interval in seen]
        if shard_intervals:
            selected.append({**shard, "intervals": shard_intervals})
    return selected


def build_dense_prs_shard(
    *,
    docker_bin: str,
    stage_root: pathlib.Path,
    cram_path: pathlib.Path,
    reference_fasta: pathlib.Path,
    shard: dict[str, Any],
    targets_vcf: pathlib.Path,
    sample_id: str,
    pairhmm_threads: int,
    java_memory: str,
) -> str:
    shard_dir = ensure_dir(stage_root / "dense" / "shards" / shard["name"])
    log_dir = ensure_dir(stage_root / "logs")
    shard_targets = shard_dir / f"{shard['name']}.targets.vcf.gz"
    gvcf = shard_dir / f"{sample_id}.{shard['name']}.g.vcf.gz"
    dense_vcf = shard_dir / f"{sample_id}.{shard['name']}.dense.vcf.gz"
    region_arg = ",".join(shard["intervals"])

    if not complete_artifact(shard_targets):
        remove_if_exists(shard_targets, f"{shard_targets}.tbi")
        run_command(
            ["bcftools", "view", "-r", region_arg, "-Oz", "-o", str(shard_targets), str(targets_vcf)],
            check=True,
            log_path=log_dir / f"{shard['name']}.targets.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(shard_targets)],
            check=True,
            log_path=log_dir / f"{shard['name']}.targets.index.log",
        )

    if not complete_artifact(gvcf):
        remove_if_exists(gvcf, f"{gvcf}.tbi")
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={
                    cram_path.parent: "/input",
                    reference_fasta.parent: "/ref",
                    stage_root: "/out",
                },
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "HaplotypeCaller",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-I", f"/input/{cram_path.name}",
                    "-L", f"/out/dense/shards/{shard['name']}/{shard_targets.name}",
                    "-O", f"/out/dense/shards/{shard['name']}/{gvcf.name}",
                    "-ERC", "GVCF",
                    "--native-pair-hmm-threads", str(pairhmm_threads),
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / f"{shard['name']}.haplotypecaller.log",
        )

    if not complete_artifact(dense_vcf):
        remove_if_exists(dense_vcf, f"{dense_vcf}.tbi")
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={
                    cram_path.parent: "/input",
                    reference_fasta.parent: "/ref",
                    stage_root: "/out",
                },
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "GenotypeGVCFs",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/dense/shards/{shard['name']}/{gvcf.name}",
                    "-L", f"/out/dense/shards/{shard['name']}/{shard_targets.name}",
                    "--include-non-variant-sites",
                    "-O", f"/out/dense/shards/{shard['name']}/{dense_vcf.name}",
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / f"{shard['name']}.genotypegvcfs.log",
        )
    return str(dense_vcf)


def concat_dense_vcf(stage_root: pathlib.Path, reference_fasta: pathlib.Path, sample_id: str, dense_shards: list[str]) -> str:
    concat_vcf = stage_root / "dense" / f"{sample_id}.prs_dense.concat.vcf.gz"
    final_vcf = stage_root / "dense" / f"{sample_id}.prs_dense.vcf.gz"
    log_dir = ensure_dir(stage_root / "logs")
    if not complete_artifact(concat_vcf):
        remove_if_exists(concat_vcf, f"{concat_vcf}.tbi")
        run_command(
            ["bcftools", "concat", "-a", "-Oz", "-o", str(concat_vcf), *dense_shards],
            check=True,
            log_path=log_dir / "dense.concat.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(concat_vcf)],
            check=True,
            log_path=log_dir / "dense.concat.index.log",
        )
    if not complete_artifact(final_vcf):
        remove_if_exists(final_vcf, f"{final_vcf}.tbi")
        run_command(
            ["bcftools", "norm", "-f", str(reference_fasta), "-m", "-both", "-Oz", "-o", str(final_vcf), str(concat_vcf)],
            check=True,
            log_path=log_dir / "dense.norm.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(final_vcf)],
            check=True,
            log_path=log_dir / "dense.norm.index.log",
        )
    return str(final_vcf)


def complete_pfile(prefix: pathlib.Path) -> bool:
    return artifacts_exist(
        [
            str(prefix) + ".pgen",
            str(prefix) + ".psam",
            str(prefix) + ".pvar.zst",
        ]
    )


def build_pgsc_input_pfile(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    dense_vcf: pathlib.Path,
) -> dict[str, Any]:
    docker_bin = context["host_checks"]["docker"] or "docker"
    ensure_docker_ready(docker_bin)

    prep_root = ensure_dir(stage_root / "compatible")
    sampleset = safe_sampleset_name(context["sample_id"])
    output_prefix = prep_root / f"GRCh38_{sampleset}_ALL"
    sex_code = plink_sex_code(context["manifest"]["analysis"]["prs"].get("sample_sex"))
    chrom_filter = "1-22, X, Y, XY" if sex_code else "1-22"
    max_cpus = max(int(context["manifest"]["runtime"]["max_cpus"]), 1)
    thread_count = max(1, min(4, max_cpus // 2 or 1))
    memory_mb = parse_memory_to_mb(context["manifest"]["runtime"]["max_memory"])
    log_path = prep_root / "plink2.make_pfile.log"

    if not complete_pfile(output_prefix):
        remove_if_exists(
            str(output_prefix) + ".pgen",
            str(output_prefix) + ".psam",
            str(output_prefix) + ".pvar.zst",
            str(output_prefix) + ".pvar",
            str(output_prefix) + ".log",
        )

        sex_file = prep_root / "sample_sex.tsv"
        args = [
            "plink2",
            "--threads", str(thread_count),
            "--memory", str(memory_mb),
            "--set-all-var-ids", "@:#:$r:$a",
            "--new-id-max-allele-len", "256",
            "--max-alleles", "2",
        ]
        if sex_code:
            sex_file.write_text(
                "#FID\tIID\tSEX\n"
                f"0\t{context['sample_id']}\t{sex_code}\n",
                encoding="utf-8",
            )
            args.extend(["--update-sex", f"/work/{sex_file.name}"])

        args.extend(
            [
                "--vcf", f"/input/{dense_vcf.name}",
                "--allow-extra-chr",
                "--chr", chrom_filter,
                "--make-pgen", "vzs", "pvar-cols=-xheader,-maybequal,-maybefilter,-maybeinfo,-maybecm",
                "--out", f"/work/{output_prefix.name}",
            ]
        )

        run_command(
            build_docker_run_command(
                image=PLINK2_IMAGE,
                mounts={
                    dense_vcf.parent: "/input",
                    prep_root: "/work",
                },
                inner_args=args,
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_path,
        )

    return {
        "path_prefix": str(output_prefix),
        "sample_sex": context["manifest"]["analysis"]["prs"].get("sample_sex"),
        "sex_chromosomes_included": bool(sex_code),
        "artifacts": {
            "pgen": str(output_prefix) + ".pgen",
            "psam": str(output_prefix) + ".psam",
            "pvar": str(output_prefix) + ".pvar.zst",
            "log": str(log_path),
        },
    }


def run_pgsc_calc(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    dense_vcf: pathlib.Path,
) -> tuple[int, pathlib.Path, dict[str, Any]]:
    ids = context["manifest"]["analysis"]["prs"]["pgs_ids"]
    samplesheet = stage_root / "samplesheet.csv"
    ids_file = stage_root / "pgs_ids.txt"
    outdir = ensure_dir(stage_root / "pgsc_calc")
    workdir = ensure_dir(stage_root / "work")
    genotypes_cache = ensure_dir(stage_root / "genotypes_cache")
    log_path = stage_root / "pgsc_calc.log"
    pfile_info = build_pgsc_input_pfile(context=context, stage_root=stage_root, dense_vcf=dense_vcf)

    samplesheet.write_text(
        "sampleset,path_prefix,chrom,format\n"
        f"{safe_sampleset_name(context['sample_id'])},{pfile_info['path_prefix']},,pfile\n",
        encoding="utf-8",
    )
    ids_file.write_text("\n".join(ids) + "\n", encoding="utf-8")

    args = [
        "nextflow",
        "run",
        "pgscatalog/pgsc_calc",
        "-profile", context["manifest"]["runtime"]["profile"],
        "--input", str(samplesheet),
        "--pgs_id", ",".join(ids),
        "--target_build", "GRCh38",
        "--min_overlap", str(context["manifest"]["analysis"]["prs"]["min_overlap"]),
        "--max_memory", context["manifest"]["runtime"]["max_memory"],
        "--outdir", str(outdir),
        "--genotypes_cache", str(genotypes_cache),
        "-work-dir", str(workdir),
        "-resume",
    ]

    env = {
        "NXF_HOME": context["runtime"]["nextflow_home"],
    }
    completed = run_command(args, cwd=context["repo_root"], env=env, log_path=log_path, check=False)
    return completed.returncode, log_path, pfile_info


def find_first(root: pathlib.Path, pattern: str) -> str | None:
    matches = sorted(root.rglob(pattern))
    return str(matches[0]) if matches else None


def load_existing_stage(stage_manifest_path: pathlib.Path, cache_key: str) -> dict[str, Any] | None:
    if not stage_manifest_path.exists():
        return None
    payload = load_json(stage_manifest_path)
    if payload.get("cache_key") != cache_key or payload.get("status") != "ok":
        return None
    artifacts = payload.get("artifacts") or {}
    if not artifacts_exist([artifacts.get("aggregated_scores")]):
        return None
    payload["reused"] = True
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    cfg = context["manifest"]["analysis"]["prs"]
    stage_root = ensure_dir(pathlib.Path(context["cache"]["prs_root"]))
    live_status_path = stage_root / "stage_status.json"
    stage_manifest_path = stage_root / "stage_manifest.json"
    scoring_root = stage_root / "scoring"
    cache_key = f"v{STAGE_VERSION}:{context['input_fingerprint']['id']}:initializing"

    try:
        progress: dict[str, Any] | None = None
        dense_info: dict[str, Any] | None = None

        def write_live_status(status: str, *, phase: str | None = None, message: str | None = None, diagnostics: dict[str, Any] | None = None, artifacts: dict[str, Any] | None = None) -> None:
            payload = stage_payload(
                "prs",
                status,
                cache_key=cache_key,
                input_fingerprint_id=context["input_fingerprint"]["id"],
                dense_genotypes=dense_info,
                progress=progress,
                phase=phase,
                message=message,
                diagnostics=diagnostics,
                artifacts=artifacts or {"log_dir": str(stage_root / "logs")},
            )
            write_json(live_status_path, {key: value for key, value in payload.items() if value is not None})

        if not cfg["enabled"]:
            payload = stage_payload("prs", "skipped", message="PRS stage disabled.")
            write_json(args.output, payload)
            return

        reference_stage = context["reference_bundle"]
        imported_dense = context["reuse"].get("prs_dense_vcf")

        if imported_dense:
            dense_vcf = pathlib.Path(imported_dense["active"])
            dense_info = stage_payload(
                "prs_dense_genotypes",
                "ok",
                reused=False,
                source_mode="import_prs_dense_vcf",
                import_details={
                    "source": imported_dense["original"],
                    "active": imported_dense["active"],
                    "provenance": imported_dense["provenance"],
                    "trusted_provenance": imported_dense["trusted_provenance"],
                    "stage_manifest": imported_dense.get("stage_manifest"),
                },
                artifacts={
                    "dense_vcf": str(dense_vcf),
                    "dense_vcf_index": imported_dense["sidecars"].get("index", str(dense_vcf) + ".tbi"),
                },
            )
            score_hash = "imported"
            progress = {"phase": "scoring_imported_dense_vcf"}
            write_live_status("running", phase="scoring_imported_dense_vcf")
        else:
            if reference_stage.get("status") != "ready":
                payload = stage_payload(
                    "prs",
                    "failed",
                    message="PRS dense genotyping requires a validated local reference bundle.",
                    dense_genotypes=stage_payload("prs_dense_genotypes", "failed", message="Reference bundle not ready."),
                )
                write_json(stage_manifest_path, payload)
                write_json(args.output, payload)
                return

            reference_fasta = pathlib.Path(reference_stage["fasta"])
            reference_fai = pathlib.Path(reference_stage["fai"])
            union_stage = build_union_sites(context=context, stage_root=stage_root, reference_fasta=reference_fasta, reference_fai=reference_fai)
            score_hash = union_stage["score_hash"]
            dense_root = ensure_dir(stage_root / score_hash)
            dense_vcf = dense_root / "dense" / f"{context['sample_id']}.prs_dense.vcf.gz"

            if complete_artifact(dense_vcf):
                dense_info = stage_payload(
                    "prs_dense_genotypes",
                    "ok",
                    reused=True,
                    score_hash=score_hash,
                    summary=union_stage.get("summary"),
                    artifacts={
                        "dense_vcf": str(dense_vcf),
                        "dense_vcf_index": str(dense_vcf) + ".tbi",
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                    },
                )
                progress = {"phase": "scoring_reused_dense_vcf"}
                write_live_status("running", phase="scoring_reused_dense_vcf")
            else:
                docker_bin = context["host_checks"]["docker"] or "docker"
                ensure_docker_ready(docker_bin)
                shards = select_prs_shards(context, pathlib.Path(union_stage["artifacts"]["targets_vcf"]))
                worker_plan = plan_worker_pool(
                    stage_name="prs_dense_genotypes",
                    shard_count=len(shards),
                    max_cpus=int(context["manifest"]["runtime"]["max_cpus"]),
                    max_memory=context["manifest"]["runtime"]["max_memory"],
                    requested_workers=cfg.get("workers"),
                    requested_pairhmm_threads=cfg.get("pairhmm_threads"),
                    default_worker_cap=4,
                    default_pairhmm_cap=2,
                    min_memory_per_worker_gb=2.5,
                    free_gb=context["disk"].get("free_gb"),
                    min_free_gb=context["disk"].get("min_free_gb"),
                )
                progress = {
                    "phase": "building_dense_genotypes",
                    "score_hash": score_hash,
                    "total_shards": len(shards),
                    "completed_shards": 0,
                    "pending_shards": len(shards),
                    "completed_shard_names": [],
                }
                dense_info = stage_payload(
                    "prs_dense_genotypes",
                    "running",
                    score_hash=score_hash,
                    summary=union_stage.get("summary"),
                    worker_plan=worker_plan,
                    artifacts={
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                    },
                )
                write_live_status("running", phase="building_dense_genotypes")

                tasks = [
                    (
                        lambda shard=shard: build_dense_prs_shard(
                            docker_bin=docker_bin,
                            stage_root=dense_root,
                            cram_path=pathlib.Path(context["inputs"]["active"]["cram"]),
                            reference_fasta=reference_fasta,
                            shard=shard,
                            targets_vcf=pathlib.Path(union_stage["artifacts"]["targets_vcf"]),
                            sample_id=context["sample_id"],
                            pairhmm_threads=worker_plan["pairhmm_threads"],
                            java_memory=worker_plan["java_memory"],
                        )
                    )
                    for shard in shards
                ]
                shard_names = [shard["name"] for shard in shards]

                def on_shard_complete(index: int, _result: Any) -> None:
                    assert progress is not None
                    shard_name = shard_names[index]
                    completed = progress["completed_shard_names"]
                    if shard_name not in completed:
                        completed.append(shard_name)
                    progress["completed_shards"] = len(completed)
                    progress["pending_shards"] = len(shards) - len(completed)
                    progress["last_completed_shard"] = shard_name
                    write_live_status("running", phase="building_dense_genotypes")

                dense_shards = run_tasks_in_parallel(tasks, max_workers=worker_plan["workers"], on_complete=on_shard_complete)
                progress = {
                    **(progress or {}),
                    "phase": "concatenating_dense_vcf",
                }
                write_live_status("running", phase="concatenating_dense_vcf")
                dense_vcf_path = concat_dense_vcf(dense_root, reference_fasta, context["sample_id"], dense_shards)
                dense_info = stage_payload(
                    "prs_dense_genotypes",
                    "ok",
                    reused=False,
                    score_hash=score_hash,
                    summary=union_stage.get("summary"),
                    worker_plan=worker_plan,
                    artifacts={
                        "dense_vcf": dense_vcf_path,
                        "dense_vcf_index": dense_vcf_path + ".tbi",
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                    },
                )
                progress = {
                    **(progress or {}),
                    "phase": "scoring_dense_vcf",
                }
                write_live_status("running", phase="scoring_dense_vcf")

        scoring_root = ensure_dir(stage_root / score_hash / "scoring")
        stage_manifest_path = scoring_root / "stage_manifest.json"
        cache_key = ":".join(
            [
                f"v{STAGE_VERSION}",
                context["input_fingerprint"]["id"],
                score_hash,
                str(cfg["min_overlap"]),
                str(bool(imported_dense)),
                str(cfg.get("sample_sex") or "autosome_only"),
            ]
        )
        reused = load_existing_stage(stage_manifest_path, cache_key)
        if reused and (reused.get("dense_genotypes") or {}).get("artifacts", {}).get("dense_vcf") == dense_info.get("artifacts", {}).get("dense_vcf"):
            write_json(live_status_path, reused)
            write_json(args.output, reused)
            return

        exit_code, log_path, pfile_info = run_pgsc_calc(context=context, stage_root=scoring_root, dense_vcf=pathlib.Path(dense_info["artifacts"]["dense_vcf"]))
        log_text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
        aggregated_scores = find_first(scoring_root / "pgsc_calc", "aggregated_scores.txt.gz")
        report_html = find_first(scoring_root / "pgsc_calc", "report.html")
        match_summary = find_first(scoring_root / "pgsc_calc", "*summary*.csv")
        overlap_failures = parse_overlap_diagnostics(log_text)

        if exit_code != 0 or not aggregated_scores:
            payload = stage_payload(
                "prs",
                "failed",
                cache_key=cache_key,
                input_fingerprint_id=context["input_fingerprint"]["id"],
                message="pgsc_calc did not produce final PRS scores.",
                dense_genotypes=dense_info,
                diagnostics={
                    "exit_code": exit_code,
                    "overlap_failures": overlap_failures,
                    "log_path": str(log_path),
                    "sample_sex": pfile_info["sample_sex"],
                    "sex_chromosomes_included": pfile_info["sex_chromosomes_included"],
                },
                artifacts={
                    "dense_vcf": dense_info["artifacts"]["dense_vcf"],
                    "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"],
                    "targets_vcf": dense_info["artifacts"].get("targets_vcf"),
                    "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index"),
                    "compatible_pgen": pfile_info["artifacts"]["pgen"],
                    "compatible_psam": pfile_info["artifacts"]["psam"],
                    "compatible_pvar": pfile_info["artifacts"]["pvar"],
                    "compatible_plink_log": pfile_info["artifacts"]["log"],
                    "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                    "pgsc_log": str(log_path),
                },
            )
            write_json(stage_manifest_path, payload)
            write_json(live_status_path, payload)
            write_json(args.output, payload)
            return

        preview_rows = []
        with gzip.open(aggregated_scores, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for idx, row in enumerate(reader):
                preview_rows.append({"PGS": row.get("PGS", ""), "IID": row.get("IID", ""), "SUM": row.get("SUM", "")})
                if idx >= 9:
                    break

        payload = stage_payload(
            "prs",
            "ok",
            cache_key=cache_key,
            input_fingerprint_id=context["input_fingerprint"]["id"],
            reused=False,
            dense_genotypes=dense_info,
            summary={
                "num_requested_scores": len(cfg["pgs_ids"]),
                "union_sites": (dense_info.get("summary") or {}).get("union_sites"),
                "preview_rows": preview_rows,
                "overlap_failures": overlap_failures,
                "sample_sex": pfile_info["sample_sex"],
                "sex_chromosomes_included": pfile_info["sex_chromosomes_included"],
            },
            artifacts={
                "dense_vcf": dense_info["artifacts"]["dense_vcf"],
                "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"],
                "targets_vcf": dense_info["artifacts"].get("targets_vcf"),
                "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index"),
                "compatible_pgen": pfile_info["artifacts"]["pgen"],
                "compatible_psam": pfile_info["artifacts"]["psam"],
                "compatible_pvar": pfile_info["artifacts"]["pvar"],
                "compatible_plink_log": pfile_info["artifacts"]["log"],
                "aggregated_scores": aggregated_scores,
                "pgsc_report_html": report_html,
                "match_summary": match_summary,
                "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                "pgsc_log": str(log_path),
            },
        )
        write_json(stage_manifest_path, payload)
        write_json(live_status_path, payload)
        write_json(args.output, payload)
    except Exception as exc:
        payload = stage_payload(
            "prs",
            "failed",
            cache_key=cache_key,
            input_fingerprint_id=context["input_fingerprint"]["id"],
            message=str(exc),
            dense_genotypes=dense_info,
            progress=progress,
            artifacts={"log_dir": str(scoring_root / "logs")},
        )
        write_json(stage_manifest_path, payload)
        write_json(live_status_path, payload)
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
