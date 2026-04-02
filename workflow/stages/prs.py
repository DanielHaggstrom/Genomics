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
    build_docker_run_command,
    complete_artifact,
    complete_artifacts_exist,
    ensure_dir,
    ensure_docker_ready,
    hash_strings,
    init_shard_progress,
    load_json,
    parse_memory_to_mb,
    plan_worker_pool,
    record_shard_completion,
    remove_if_exists,
    run_command,
    run_tasks_in_parallel,
    stage_payload,
    write_json,
)


STAGE_VERSION = 7
UNION_SITE_VERSION = 3
GATK_IMAGE = "broadinstitute/gatk:4.6.1.0"
PLINK2_IMAGE = "ghcr.io/pgscatalog/plink2:2.00a5.10"
OVERLAP_RE = re.compile(r"Score (\S+) fails minimum matching threshold \(([\d.]+)% variants match\)")
BASE_COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


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

    def get_sequence(self, contig: str, pos: int, sequence_length: int) -> str:
        return "".join(self.get_base(contig, pos + offset) for offset in range(sequence_length))


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


def canonical_score_id(score_name: str) -> str:
    return re.sub(r"_hmPOS_GRCh\d+$", "", score_name)


def init_score_audit(pgs_ids: list[str]) -> dict[str, dict[str, int]]:
    audit: dict[str, dict[str, int]] = {}
    for pgs_id in pgs_ids:
        audit[pgs_id] = {
            "accepted_rows": 0,
            "accepted_explicit_other_allele": 0,
            "accepted_hm_inferOtherAllele": 0,
            "accepted_multiallelic_rows": 0,
            "skipped_rows": 0,
        }
    return audit


def increment_score_skip(score_audit: dict[str, dict[str, int]], pgs_id: str, reason: str) -> None:
    bucket = score_audit[pgs_id]
    bucket["skipped_rows"] += 1
    key = f"skipped_{reason}"
    bucket[key] = bucket.get(key, 0) + 1


def record_score_accept(score_audit: dict[str, dict[str, int]], pgs_id: str, resolved: dict[str, Any]) -> None:
    bucket = score_audit[pgs_id]
    bucket["accepted_rows"] += 1
    bucket[f"accepted_{resolved['source']}"] = bucket.get(f"accepted_{resolved['source']}", 0) + 1
    if resolved["multiallelic"]:
        bucket["accepted_multiallelic_rows"] += 1
    if resolved.get("strand_flipped"):
        bucket["accepted_strand_flip_rows"] = bucket.get("accepted_strand_flip_rows", 0) + 1


def write_tsv(path: pathlib.Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_pre_match_audit(
    *,
    targets_root: pathlib.Path,
    score_hash: str,
    pgs_ids: list[str],
    union_summary: dict[str, Any],
    score_audit: dict[str, dict[str, int]],
) -> dict[str, str]:
    audit_root = ensure_dir(targets_root / "diagnostics")
    json_path = audit_root / "pre_match_audit.json"
    tsv_path = audit_root / "pre_match_audit.tsv"

    ordered_rows = []
    for pgs_id in pgs_ids:
        bucket = score_audit[pgs_id]
        row = {
            "pgs_id": pgs_id,
            **bucket,
        }
        ordered_rows.append(row)

    fieldnames = sorted({key for row in ordered_rows for key in row.keys()}, key=lambda item: (item != "pgs_id", item))
    write_tsv(tsv_path, fieldnames, ordered_rows)
    write_json(
        json_path,
        {
            "stage": "prs_pre_match_audit",
            "status": "ok",
            "score_hash": score_hash,
            "summary": union_summary,
            "rows": ordered_rows,
        },
    )
    return {
        "pre_match_audit_json": str(json_path),
        "pre_match_audit_tsv": str(tsv_path),
    }


def write_overlap_failure_artifacts(
    *,
    scoring_root: pathlib.Path,
    min_overlap: float,
    overlap_failures: list[dict[str, str]],
    log_path: pathlib.Path,
    match_summary: str | None,
) -> dict[str, str] | dict[str, str | None]:
    diagnostics_root = ensure_dir(scoring_root / "diagnostics")
    json_path = diagnostics_root / "match_diagnostics.json"
    tsv_path = diagnostics_root / "overlap_failures.tsv"

    ordered_rows = []
    for item in overlap_failures:
        ordered_rows.append(
            {
                "score": item["score"],
                "pgs_id": canonical_score_id(item["score"]),
                "match_percent": item["match_percent"],
                "min_overlap": str(min_overlap),
                "passed_threshold": "false",
            }
        )

    if ordered_rows:
        write_tsv(tsv_path, ["score", "pgs_id", "match_percent", "min_overlap", "passed_threshold"], ordered_rows)
    else:
        write_tsv(tsv_path, ["score", "pgs_id", "match_percent", "min_overlap", "passed_threshold"], [])

    write_json(
        json_path,
        {
            "stage": "prs_match_diagnostics",
            "status": "ok",
            "min_overlap": min_overlap,
            "overlap_failures": ordered_rows,
            "pgsc_log": str(log_path),
            "match_summary": match_summary,
        },
    )
    return {
        "match_diagnostics_json": str(json_path),
        "overlap_failures_tsv": str(tsv_path),
    }


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
    return {token for token in tokens if token and re.fullmatch(r"[ACGT]+", token)}


def complement_base_tokens(tokens: set[str]) -> set[str]:
    complemented = set()
    for token in tokens:
        complement = "".join(BASE_COMPLEMENTS.get(base, "") for base in token)
        if len(complement) != len(token):
            return set()
        complemented.add(complement)
    return complemented


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
    strand_flipped = False
    if ref not in all_alleles:
        complemented_all = complement_base_tokens(all_alleles)
        if ref not in complemented_all:
            return None, "reference_not_in_score_alleles"
        strand_flipped = True
        effect = BASE_COMPLEMENTS[effect]
        other_alleles = complement_base_tokens(other_alleles)
        all_alleles = complemented_all

    alts = all_alleles - {ref}
    if not alts:
        return None, "no_alt_after_reference_match"

    return {
        "source": source,
        "ref_allele": ref,
        "effect_allele": effect,
        "other_alleles": sorted(other_alleles),
        "alts": sorted(alts),
        "multiallelic": len(alts) > 1,
        "strand_flipped": strand_flipped,
    }, None


def resolve_score_row_alleles_for_site(
    row: dict[str, str],
    *,
    fasta: Any,
    contig: str,
    pos: int,
) -> tuple[dict[str, Any] | None, str | None]:
    effect_alleles = parse_base_tokens(row.get("effect_allele"))
    if len(effect_alleles) != 1:
        return None, "unsupported_effect_allele"
    effect = next(iter(effect_alleles))
    complement_effect_set = complement_base_tokens({effect})
    complement_effect = next(iter(complement_effect_set)) if complement_effect_set else ""

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
            source = "effect_allele_only"
            try:
                ref_allele = fasta.get_sequence(contig, pos, len(effect))
            except Exception:
                return None, "reference_lookup_failed"
            strand_flipped = False
            oriented_effect = effect
            if effect == ref_allele:
                alts: set[str] = set()
            elif complement_effect and complement_effect == ref_allele:
                strand_flipped = True
                oriented_effect = complement_effect
                alts = set()
            else:
                alts = {oriented_effect}
            return {
                "source": source,
                "ref_allele": ref_allele,
                "effect_allele": oriented_effect,
                "other_alleles": [],
                "alts": sorted(alts),
                "multiallelic": len(alts) > 1,
                "strand_flipped": strand_flipped,
            }, None

    all_alleles = set(other_alleles)
    all_alleles.add(effect)

    try:
        reference_sequences = {
            len(allele): fasta.get_sequence(contig, pos, len(allele))
            for allele in all_alleles
        }
    except Exception:
        return None, "reference_lookup_failed"

    for strand_flipped, oriented_effect, oriented_other in [
        (False, effect, set(other_alleles)),
        (True, complement_effect, complement_base_tokens(other_alleles)),
    ]:
        if strand_flipped and (not oriented_effect or not oriented_other and other_alleles):
            continue
        oriented_all = set(oriented_other)
        oriented_all.add(oriented_effect)
        matching_refs = [
            allele
            for allele in oriented_all
            if reference_sequences.get(len(allele)) == allele
        ]
        if not matching_refs:
            continue
        ref_allele = sorted(matching_refs, key=len, reverse=True)[0]
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": oriented_effect,
                "other_allele": "/".join(sorted(oriented_other)),
            },
            ref_allele,
        )
        if resolved is None:
            if reason == "reference_not_in_score_alleles":
                continue
            return None, reason
        resolved["source"] = source
        resolved["ref_allele"] = ref_allele
        resolved["strand_flipped"] = strand_flipped
        return resolved, None

    return None, "reference_not_in_score_alleles"


def fetch_scorefile_url(pgs_id: str) -> str:
    with urllib.request.urlopen(f"https://www.pgscatalog.org/rest/score/{pgs_id}", timeout=60) as response:
        payload = json.load(response)
    return payload["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]


def resolve_cached_scorefile(score_cache_root: pathlib.Path, pgs_id: str) -> pathlib.Path | None:
    exact = score_cache_root / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if exact.exists():
        return exact
    matches = sorted(score_cache_root.glob(f"{pgs_id}_*.txt.gz"))
    return matches[0] if matches else None


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
        required_artifacts = [
            artifacts.get("targets_vcf"),
            artifacts.get("pre_match_audit_tsv"),
            artifacts.get("pre_match_audit_json"),
        ]
        if payload.get("reuse_key") == reuse_key and complete_artifacts_exist(required_artifacts):
            return payload

    score_cache_root = ensure_dir(pathlib.Path(context["runtime"]["tool_cache_root"]) / "pgs_scorefiles" / "GRCh38")
    fasta = IndexedFasta(reference_fasta, reference_fai)
    contig_order = {contig: idx for idx, contig in enumerate(fasta.index.keys())}
    site_map: dict[tuple[str, int, str], dict[str, object]] = {}
    per_score_counts: dict[str, int] = {}
    per_score_audit = init_score_audit(pgs_ids)
    accepted_rows_by_source = {
        "explicit_other_allele": 0,
        "hm_inferOtherAllele": 0,
        "effect_allele_only": 0,
    }
    multiallelic_rows = 0
    strand_flip_rows = 0
    skipped = 0
    skipped_reasons: dict[str, int] = {}

    for pgs_id in pgs_ids:
        score_path = resolve_cached_scorefile(score_cache_root, pgs_id)
        if score_path is None:
            url = fetch_scorefile_url(pgs_id)
            score_path = score_cache_root / pathlib.Path(url).name
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
                    increment_score_skip(per_score_audit, pgs_id, "missing_position")
                    continue

                chrom = normalize_contig(chrom_raw, prefix)
                pos = int(pos_raw)
                resolved, reason = resolve_score_row_alleles_for_site(
                    row,
                    fasta=fasta,
                    contig=chrom,
                    pos=pos,
                )
                if resolved is None:
                    skipped += 1
                    skipped_reasons[reason] = skipped_reasons.get(reason, 0) + 1
                    increment_score_skip(per_score_audit, pgs_id, reason)
                    continue

                key = (chrom, pos, resolved["ref_allele"])
                entry = site_map.setdefault(key, {"ref": resolved["ref_allele"], "alts": set(), "pgs_ids": set()})
                entry["alts"].update(resolved["alts"])
                entry["pgs_ids"].add(pgs_id)
                per_score_counts[pgs_id] += 1
                accepted_rows_by_source[resolved["source"]] = accepted_rows_by_source.get(resolved["source"], 0) + 1
                record_score_accept(per_score_audit, pgs_id, resolved)
                if resolved["multiallelic"]:
                    multiallelic_rows += 1
                if resolved.get("strand_flipped"):
                    strand_flip_rows += 1

    plain_vcf = targets_root / "pgs_union_sites.vcf"
    with open(plain_vcf, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, meta in fasta.index.items():
            out.write(f"##contig=<ID={contig},length={meta[0]}>\n")
        out.write("##INFO=<ID=PGS_IDS,Number=.,Type=String,Description=\"PGS IDs contributing to this union site\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref_allele in sorted(site_map.keys(), key=lambda item: (contig_order[item[0]], item[1], item[2])):
            entry = site_map[(chrom, pos, ref_allele)]
            alt = ",".join(sorted(entry["alts"])) or "."
            info = "PGS_IDS=" + ",".join(sorted(entry["pgs_ids"]))
            out.write(f"{chrom}\t{pos}\t.\t{entry['ref']}\t{alt}\t.\tPASS\t{info}\n")

    targets_vcf = targets_root / "pgs_union_sites.vcf.gz"
    run_command(["bcftools", "view", "-Oz", "-o", str(targets_vcf), str(plain_vcf)], check=True)
    run_command(["bcftools", "index", "-t", str(targets_vcf)], check=True)

    summary = {
        "union_sites": len(site_map),
        "skipped_rows": skipped,
        "skipped_reasons": skipped_reasons,
        "per_score_counts": per_score_counts,
        "accepted_rows_by_source": accepted_rows_by_source,
        "accepted_multiallelic_rows": multiallelic_rows,
        "accepted_strand_flip_rows": strand_flip_rows,
    }
    diagnostic_artifacts = write_pre_match_audit(
        targets_root=targets_root,
        score_hash=score_hash,
        pgs_ids=pgs_ids,
        union_summary=summary,
        score_audit=per_score_audit,
    )
    payload = {
        "stage": "prs_union_sites",
        "status": "ok",
        "reuse_key": reuse_key,
        "score_hash": score_hash,
        "summary": summary,
        "artifacts": {
            "targets_vcf": str(targets_vcf),
            "targets_vcf_index": str(targets_vcf) + ".tbi",
            "plain_vcf": str(plain_vcf),
            **diagnostic_artifacts,
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
) -> dict[str, Any]:
    shard_dir = ensure_dir(stage_root / "dense" / "shards" / shard["name"])
    log_dir = ensure_dir(stage_root / "logs")
    shard_targets = shard_dir / f"{shard['name']}.targets.vcf.gz"
    gvcf = shard_dir / f"{sample_id}.{shard['name']}.g.vcf.gz"
    dense_vcf = shard_dir / f"{sample_id}.{shard['name']}.dense.vcf.gz"
    region_arg = ",".join(shard["intervals"])
    shard_targets_reused = complete_artifact(shard_targets)
    gvcf_reused = complete_artifact(gvcf)
    dense_vcf_reused = complete_artifact(dense_vcf)

    if not shard_targets_reused:
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

    if not gvcf_reused:
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

    if not dense_vcf_reused:
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
    return {
        "shard_name": shard["name"],
        "dense_vcf": str(dense_vcf),
        "gvcf": str(gvcf),
        "targets_vcf": str(shard_targets),
        "reused": dense_vcf_reused,
        "reused_targets_vcf": shard_targets_reused,
        "reused_gvcf": gvcf_reused,
        "reused_dense_vcf": dense_vcf_reused,
    }


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
    return complete_artifacts_exist(
        [
            str(prefix) + ".pgen",
            str(prefix) + ".psam",
            str(prefix) + ".pvar.zst",
        ]
    )


def imported_stage_payload(imported: dict[str, Any], *, stage: str, source_mode: str, artifact_key: str) -> dict[str, Any]:
    target_path = pathlib.Path(imported["active"])
    artifacts = {
        artifact_key: str(target_path),
        f"{artifact_key}_index": imported["sidecars"].get("index", str(target_path) + ".tbi"),
    }
    return stage_payload(
        stage,
        "ok",
        reused=False,
        source_mode=source_mode,
        import_details={
            "source": imported["original"],
            "active": imported["active"],
            "provenance": imported["provenance"],
            "trusted_provenance": imported["trusted_provenance"],
            "stage_manifest": imported.get("stage_manifest"),
            "fingerprint": imported.get("fingerprint"),
        },
        artifacts=artifacts,
    )


def build_pgsc_input_pfile(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    input_vcf: pathlib.Path,
    input_kind: str,
    dosage_field: str | None = None,
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
                "--vcf", f"/input/{input_vcf.name}",
                "--allow-extra-chr",
                "--chr", chrom_filter,
                "--make-pgen", "vzs", "pvar-cols=-xheader,-maybequal,-maybefilter,-maybeinfo,-maybecm",
                "--out", f"/work/{output_prefix.name}",
            ]
        )
        if dosage_field:
            args.insert(args.index("--allow-extra-chr"), f"dosage={dosage_field}")

        run_command(
            build_docker_run_command(
                image=PLINK2_IMAGE,
                mounts={
                    input_vcf.parent: "/input",
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
        "input_kind": input_kind,
        "dosage_field": dosage_field,
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
    scoring_vcf: pathlib.Path,
    input_kind: str,
    dosage_field: str | None = None,
) -> tuple[int, pathlib.Path, dict[str, Any]]:
    ids = context["manifest"]["analysis"]["prs"]["pgs_ids"]
    samplesheet = stage_root / "samplesheet.csv"
    ids_file = stage_root / "pgs_ids.txt"
    outdir = ensure_dir(stage_root / "pgsc_calc")
    workdir = ensure_dir(stage_root / "work")
    genotypes_cache = ensure_dir(stage_root / "genotypes_cache")
    log_path = stage_root / "pgsc_calc.log"
    pfile_info = build_pgsc_input_pfile(
        context=context,
        stage_root=stage_root,
        input_vcf=scoring_vcf,
        input_kind=input_kind,
        dosage_field=dosage_field,
    )

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
    required_artifacts = [
        artifacts.get("aggregated_scores"),
        artifacts.get("scoring_input_vcf"),
        artifacts.get("compatible_pgen"),
        artifacts.get("compatible_psam"),
        artifacts.get("compatible_pvar"),
    ]
    if not complete_artifacts_exist(required_artifacts):
        return None
    dense_artifacts = (payload.get("dense_genotypes") or {}).get("artifacts") or {}
    imputed_artifacts = (payload.get("imputed_genotypes") or {}).get("artifacts") or {}
    if not complete_artifacts_exist([dense_artifacts.get("dense_vcf")]) and not complete_artifacts_exist([imputed_artifacts.get("imputed_vcf")]):
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
        imputed_info: dict[str, Any] | None = None

        def write_live_status(status: str, *, phase: str | None = None, message: str | None = None, diagnostics: dict[str, Any] | None = None, artifacts: dict[str, Any] | None = None) -> None:
            payload = stage_payload(
                "prs",
                status,
                cache_key=cache_key,
                input_fingerprint_id=context["input_fingerprint"]["id"],
                dense_genotypes=dense_info,
                imputed_genotypes=imputed_info,
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
        imported_imputed = context["reuse"].get("imputed_prs_vcf")
        imputation_cfg = cfg.get("imputation") or {}

        if imported_dense and imported_imputed:
            payload = stage_payload(
                "prs",
                "failed",
                message="Configure at most one PRS input import: reuse.prs_dense_vcf or reuse.imputed_prs_vcf.",
            )
            write_json(stage_manifest_path, payload)
            write_json(args.output, payload)
            return

        if imported_imputed:
            imputed_info = imported_stage_payload(
                imported_imputed,
                stage="prs_imputed_genotypes",
                source_mode="import_imputed_prs_vcf",
                artifact_key="imputed_vcf",
            )
            score_hash = "imputed-" + imported_imputed["fingerprint"]["id"][:16]
            progress = {"phase": "scoring_imported_imputed_vcf"}
            write_live_status("running", phase="scoring_imported_imputed_vcf")
        elif imported_dense:
            dense_info = imported_stage_payload(
                imported_dense,
                stage="prs_dense_genotypes",
                source_mode="import_prs_dense_vcf",
                artifact_key="dense_vcf",
            )
            score_hash = "dense-" + imported_dense["fingerprint"]["id"][:16]
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
                        "pre_match_audit_json": union_stage["artifacts"].get("pre_match_audit_json"),
                        "pre_match_audit_tsv": union_stage["artifacts"].get("pre_match_audit_tsv"),
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
                progress = init_shard_progress(
                    phase="building_dense_genotypes",
                    total_shards=len(shards),
                    score_hash=score_hash,
                )
                dense_info = stage_payload(
                    "prs_dense_genotypes",
                    "running",
                    score_hash=score_hash,
                    summary=union_stage.get("summary"),
                    worker_plan=worker_plan,
                    artifacts={
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                        "pre_match_audit_json": union_stage["artifacts"].get("pre_match_audit_json"),
                        "pre_match_audit_tsv": union_stage["artifacts"].get("pre_match_audit_tsv"),
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

                def on_shard_complete(index: int, result: Any) -> None:
                    assert progress is not None
                    shard_name = result.get("shard_name") or shard_names[index]
                    record_shard_completion(progress, shard_name=shard_name, reused=bool(result.get("reused")))
                    write_live_status("running", phase="building_dense_genotypes")

                shard_results = run_tasks_in_parallel(tasks, max_workers=worker_plan["workers"], on_complete=on_shard_complete)
                dense_shards = [str(item["dense_vcf"]) for item in shard_results]
                progress = {
                    **(progress or {}),
                    "phase": "concatenating_dense_vcf",
                    "reused_dense_shards": sum(1 for item in shard_results if item.get("reused")),
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
                        "pre_match_audit_json": union_stage["artifacts"].get("pre_match_audit_json"),
                        "pre_match_audit_tsv": union_stage["artifacts"].get("pre_match_audit_tsv"),
                    },
                )
                progress = {
                    **(progress or {}),
                    "phase": "scoring_dense_vcf",
                }
                write_live_status("running", phase="scoring_dense_vcf")

            if imputation_cfg.get("enabled"):
                payload = stage_payload(
                    "prs",
                    "failed",
                    message=(
                        "Local PRS imputation backend is not implemented yet. "
                        "Provide reuse.imputed_prs_vcf to score an imputed PRS VCF, "
                        "or disable analysis.prs.imputation.enabled to score the direct dense PRS VCF."
                    ),
                    dense_genotypes=dense_info,
                    imputed_genotypes=stage_payload(
                        "prs_imputed_genotypes",
                        "failed",
                        message="No imported imputed PRS VCF provided, and local imputation backend is not implemented yet.",
                    ),
                )
                write_json(stage_manifest_path, payload)
                write_json(live_status_path, payload)
                write_json(args.output, payload)
                return

        scoring_root = ensure_dir(stage_root / score_hash / "scoring")
        stage_manifest_path = scoring_root / "stage_manifest.json"
        scoring_input = imputed_info or dense_info
        scoring_input_artifacts = (scoring_input or {}).get("artifacts") or {}
        scoring_input_path = scoring_input_artifacts.get("imputed_vcf") or scoring_input_artifacts.get("dense_vcf")
        scoring_input_index = scoring_input_artifacts.get("imputed_vcf_index") or scoring_input_artifacts.get("dense_vcf_index")
        scoring_input_kind = "imputed_vcf" if imputed_info else "dense_vcf"
        scoring_input_fingerprint = (imported_imputed or imported_dense or {}).get("fingerprint", {}).get("id") or score_hash
        dosage_field = imputation_cfg.get("dosage_field") if imputed_info else None
        cache_key = ":".join(
            [
                f"v{STAGE_VERSION}",
                context["input_fingerprint"]["id"],
                score_hash,
                str(cfg["min_overlap"]),
                scoring_input_kind,
                str(scoring_input_fingerprint),
                str(cfg.get("sample_sex") or "autosome_only"),
            ]
        )
        reused = load_existing_stage(stage_manifest_path, cache_key)
        if reused and (reused.get("artifacts") or {}).get("scoring_input_vcf") == scoring_input_path:
            write_json(live_status_path, reused)
            write_json(args.output, reused)
            return

        exit_code, log_path, pfile_info = run_pgsc_calc(
            context=context,
            stage_root=scoring_root,
            scoring_vcf=pathlib.Path(scoring_input_path),
            input_kind=scoring_input_kind,
            dosage_field=dosage_field,
        )
        log_text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
        aggregated_scores = find_first(scoring_root / "pgsc_calc", "aggregated_scores.txt.gz")
        report_html = find_first(scoring_root / "pgsc_calc", "report.html")
        match_summary = find_first(scoring_root / "pgsc_calc", "*summary*.csv")
        overlap_failures = parse_overlap_diagnostics(log_text)
        match_diagnostic_artifacts = write_overlap_failure_artifacts(
            scoring_root=scoring_root,
            min_overlap=float(cfg["min_overlap"]),
            overlap_failures=overlap_failures,
            log_path=log_path,
            match_summary=match_summary,
        )

        if exit_code != 0 or not aggregated_scores:
            payload = stage_payload(
                "prs",
                "failed",
                cache_key=cache_key,
                input_fingerprint_id=context["input_fingerprint"]["id"],
                message="pgsc_calc did not produce final PRS scores.",
                dense_genotypes=dense_info,
                imputed_genotypes=imputed_info,
                diagnostics={
                    "exit_code": exit_code,
                    "overlap_failures": overlap_failures,
                    "log_path": str(log_path),
                    "sample_sex": pfile_info["sample_sex"],
                    "sex_chromosomes_included": pfile_info["sex_chromosomes_included"],
                    "scoring_input_kind": scoring_input_kind,
                    "dosage_field": pfile_info.get("dosage_field"),
                },
                artifacts={
                    "scoring_input_vcf": scoring_input_path,
                    "scoring_input_index": scoring_input_index,
                    "dense_vcf": dense_info["artifacts"]["dense_vcf"] if dense_info else None,
                    "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"] if dense_info else None,
                    "targets_vcf": dense_info["artifacts"].get("targets_vcf") if dense_info else None,
                    "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index") if dense_info else None,
                    "pre_match_audit_json": dense_info["artifacts"].get("pre_match_audit_json") if dense_info else None,
                    "pre_match_audit_tsv": dense_info["artifacts"].get("pre_match_audit_tsv") if dense_info else None,
                    "imputed_vcf": imputed_info["artifacts"]["imputed_vcf"] if imputed_info else None,
                    "imputed_vcf_index": imputed_info["artifacts"]["imputed_vcf_index"] if imputed_info else None,
                    "compatible_pgen": pfile_info["artifacts"]["pgen"],
                    "compatible_psam": pfile_info["artifacts"]["psam"],
                    "compatible_pvar": pfile_info["artifacts"]["pvar"],
                    "compatible_plink_log": pfile_info["artifacts"]["log"],
                    "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                    "pgsc_log": str(log_path),
                    **match_diagnostic_artifacts,
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
            imputed_genotypes=imputed_info,
            progress=progress,
            summary={
                "num_requested_scores": len(cfg["pgs_ids"]),
                "union_sites": (dense_info.get("summary") or {}).get("union_sites") if dense_info else None,
                "preview_rows": preview_rows,
                "overlap_failures": overlap_failures,
                "sample_sex": pfile_info["sample_sex"],
                "sex_chromosomes_included": pfile_info["sex_chromosomes_included"],
                "scoring_input_kind": scoring_input_kind,
            },
            artifacts={
                "scoring_input_vcf": scoring_input_path,
                "scoring_input_index": scoring_input_index,
                "dense_vcf": dense_info["artifacts"]["dense_vcf"] if dense_info else None,
                "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"] if dense_info else None,
                "targets_vcf": dense_info["artifacts"].get("targets_vcf") if dense_info else None,
                "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index") if dense_info else None,
                "pre_match_audit_json": dense_info["artifacts"].get("pre_match_audit_json") if dense_info else None,
                "pre_match_audit_tsv": dense_info["artifacts"].get("pre_match_audit_tsv") if dense_info else None,
                "imputed_vcf": imputed_info["artifacts"]["imputed_vcf"] if imputed_info else None,
                "imputed_vcf_index": imputed_info["artifacts"]["imputed_vcf_index"] if imputed_info else None,
                "compatible_pgen": pfile_info["artifacts"]["pgen"],
                "compatible_psam": pfile_info["artifacts"]["psam"],
                "compatible_pvar": pfile_info["artifacts"]["pvar"],
                "compatible_plink_log": pfile_info["artifacts"]["log"],
                "aggregated_scores": aggregated_scores,
                "pgsc_report_html": report_html,
                "match_summary": match_summary,
                "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                "pgsc_log": str(log_path),
                **match_diagnostic_artifacts,
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
            imputed_genotypes=imputed_info,
            progress=progress,
            artifacts={"log_dir": str(scoring_root / "logs")},
        )
        write_json(stage_manifest_path, payload)
        write_json(live_status_path, payload)
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
