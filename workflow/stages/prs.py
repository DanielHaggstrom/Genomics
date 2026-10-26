#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pathlib
import re
import shlex
import subprocess
import sys
from typing import Any, Callable

from workflow.core.runtime import (
    REPO_ROOT,
    build_docker_run_command,
    completed_vcf,
    record_vcf_checkpoint,
    fetch_url_text,
    download_url_to_path,
    checkpoint_complete,
    record_checkpoint,
    complete_artifacts_exist,
    ensure_dir,
    ensure_docker_ready,
    germline_gvcf_dir,
    hash_strings,
    init_shard_progress,
    load_json,
    load_yaml,
    normalize_path,
    needs_staging,
    parse_memory_to_mb,
    plan_worker_pool,
    record_shard_completion,
    remove_if_exists,
    run_command,
    run_tasks_in_parallel,
    sample_file_fingerprint,
    stage_file,
    stage_payload,
    write_json,
)


STAGE_VERSION = 9
UNION_SITE_VERSION = 5
IMPUTATION_STAGE_VERSION = 3
PGSC_CALC_REVISION = "v2.3.0"
GATK_IMAGE = "broadinstitute/gatk:4.6.1.0"
PLINK2_IMAGE = "ghcr.io/pgscatalog/plink2:2.00a5.10"
BEAGLE_URL = "https://faculty.washington.edu/browning/beagle/beagle.29Oct24.c8e.jar"
REFERENCE_DB_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources"
OVERLAP_RE = re.compile(r"Score (\S+) fails minimum matching threshold \(([\d.]+)% variants match\)")
CONTIG_HEADER_RE = re.compile(r"^##contig=<ID=([^,>]+)")
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


def is_valid_gzip_file(path: pathlib.Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except (OSError, EOFError):
        return False


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


def base_contig_name(contig: str, prefix: str) -> str:
    if prefix == "chr" and contig.startswith("chr"):
        return contig[3:]
    return contig


def inject_missing_contig_headers(header_text: str, contigs: list[str]) -> str:
    lines = header_text.splitlines()
    declared = set()
    for line in lines:
        found = CONTIG_HEADER_RE.match(line)
        if found:
            declared.add(found.group(1))
    missing = [contig for contig in contigs if contig and contig not in declared]
    if not missing:
        return header_text

    output: list[str] = []
    inserted = False
    for line in lines:
        if not inserted and line.startswith("#CHROM"):
            output.extend(f"##contig=<ID={contig}>" for contig in missing)
            inserted = True
        output.append(line)
    if not inserted:
        output.extend(f"##contig=<ID={contig}>" for contig in missing)
    return "\n".join(output) + "\n"


def is_manifest_path(raw: str) -> bool:
    return any(token in raw for token in ("/", "\\")) or raw.endswith(".yaml")


def allowed_prs_base_chromosomes(context: dict[str, Any]) -> list[str]:
    # The direct materializer currently validates diploid autosomes only. Imported
    # genotypes can retain sex chromosomes once ploidy was established upstream.
    return [str(i) for i in range(1, 23)]


def allowed_prs_contigs(context: dict[str, Any]) -> set[str]:
    prefix = context["contig_prefix"]
    return {normalize_contig(chrom, prefix) for chrom in allowed_prs_base_chromosomes(context)}


def resolve_reference_db_archive_url(reference_db: str) -> str:
    cleaned = reference_db.strip()
    if cleaned.startswith(("http://", "https://")):
        return cleaned
    return f"{REFERENCE_DB_BASE_URL}/{cleaned}.tar.zst"


def locate_reference_db_prefix(extracted_root: pathlib.Path, build: str) -> pathlib.Path:
    matches = sorted(extracted_root.rglob(f"{build}_*_ALL.pvar.zst"))
    if not matches:
        raise RuntimeError(f"Could not locate a {build} reference-db pvar.zst in {extracted_root}")
    pvar_path = matches[0]
    if not pvar_path.name.endswith(".pvar.zst"):
        raise RuntimeError(f"Unexpected reference-db pvar name: {pvar_path.name}")
    prefix_name = pvar_path.name[: -len(".pvar.zst")]
    prefix_path = pvar_path.with_name(prefix_name)
    for suffix in (".pvar.zst", ".pgen", ".psam"):
        candidate = prefix_path.with_name(prefix_path.name + suffix)
        if not candidate.exists():
            raise RuntimeError(f"Reference database is missing expected file: {candidate}")
    return prefix_path


def extract_tar_zst_archive(archive_path: pathlib.Path, extracted_root: pathlib.Path, *, log_path: pathlib.Path) -> None:
    run_command(["tar", "--zstd", "-xf", str(archive_path), "-C", str(extracted_root)],
                check=True, log_path=log_path)


def ensure_reference_db_assets(context: dict[str, Any]) -> dict[str, Any]:
    reference_db = str(context["manifest"]["analysis"]["prs"].get("reference_db") or "").strip()
    if not reference_db:
        raise RuntimeError("analysis.prs.reference_db is required when PRS is enabled.")

    tool_cache_root = pathlib.Path(context["runtime"]["tool_cache_root"])
    reference_db_root = ensure_dir(tool_cache_root / "pgs_reference_db")

    if is_manifest_path(reference_db) or reference_db.endswith(".tar.zst"):
        archive_path = pathlib.Path(normalize_path(reference_db, base_dir=pathlib.Path(context["manifest_path"]).parent))
        reference_db_id = archive_path.name.replace(".tar.zst", "")
        if not archive_path.exists():
            raise RuntimeError(f"Reference database archive not found: {archive_path}")
    else:
        reference_db_id = reference_db
        archive_path = reference_db_root / f"{reference_db_id}.tar.zst"
        if not archive_path.exists():
            download_url_to_path(resolve_reference_db_archive_url(reference_db_id), archive_path)

    archive_fingerprint = sample_file_fingerprint(archive_path)["id"][:16]
    extracted_root = ensure_dir(reference_db_root / f"{reference_db_id}-{archive_fingerprint}")
    ready_marker = extracted_root / ".ready"
    if not ready_marker.exists():
        extract_tar_zst_archive(
            archive_path,
            extracted_root,
            log_path=extracted_root / "extract.log",
        )
        ready_marker.write_text(archive_fingerprint + "\n", encoding="utf-8")

    prefix_path = locate_reference_db_prefix(extracted_root, context["manifest"]["reference"]["build"])
    return {
        "reference_db_id": reference_db_id,
        "archive_path": str(archive_path),
        "archive_fingerprint": archive_fingerprint,
        "extracted_root": str(extracted_root),
        "prefix_path": str(prefix_path),
    }


def export_reference_db_pvar_for_contig(
    *,
    context: dict[str, Any],
    reference_db: dict[str, Any],
    targets_root: pathlib.Path,
    base_chrom: str,
) -> pathlib.Path:
    docker_bin = context["host_checks"]["docker"] or "docker"
    ensure_docker_ready(docker_bin)

    output_root = ensure_dir(targets_root / "reference_db")
    pvar_path = output_root / f"{base_chrom}.pvar"
    receipt = pvar_path.with_suffix(".completed.json")
    if checkpoint_complete(receipt, [pvar_path]):
        return pvar_path
    remove_if_exists(pvar_path, receipt)

    extracted_root = pathlib.Path(reference_db["extracted_root"])
    prefix_path = pathlib.Path(reference_db["prefix_path"])
    relative_prefix = prefix_path.relative_to(extracted_root).as_posix()
    run_command(
        build_docker_run_command(
            image=PLINK2_IMAGE,
            mounts={
                extracted_root: ("/refdb", "ro"),
                output_root: "/out",
            },
            inner_args=[
                "plink2",
                "--pfile", f"/refdb/{relative_prefix}", "vzs",
                "--allow-extra-chr",
                "--chr", base_chrom,
                "--make-just-pvar",
                "--out", f"/out/{base_chrom}",
            ],
            docker_bin=docker_bin,
        ),
        check=True,
        log_path=output_root / f"{base_chrom}.plink2.log",
    )
    if not pvar_path.exists():
        raise RuntimeError(f"Reference database export did not produce {pvar_path}")
    record_checkpoint(receipt, [pvar_path])
    return pvar_path


def parse_pvar_header(raw_header: str) -> dict[str, int]:
    fields = raw_header.lstrip("#").rstrip("\n").split("\t")
    index = {name: idx for idx, name in enumerate(fields)}
    required = ("CHROM", "POS", "REF", "ALT")
    missing = [name for name in required if name not in index]
    if missing:
        raise RuntimeError(f"Reference database pvar header is missing columns: {', '.join(missing)}")
    return index


def iter_reference_db_sites(pvar_path: pathlib.Path, *, target_prefix: str):
    with open(pvar_path, "r", encoding="utf-8") as handle:
        header_index: dict[str, int] | None = None
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#"):
                header_index = parse_pvar_header(raw)
                continue
            if header_index is None:
                raise RuntimeError(f"Reference database pvar header not found in {pvar_path}")
            fields = raw.rstrip("\n").split("\t")
            chrom = normalize_contig(fields[header_index["CHROM"]], target_prefix)
            pos = int(fields[header_index["POS"]])
            ref = fields[header_index["REF"]]
            alts = [token for token in fields[header_index["ALT"]].split(",") if token and token != "."]
            yield {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alts": alts,
            }


def write_scaffold_plain_vcf(
    *,
    destination: pathlib.Path,
    contigs: list[str],
    contig_lengths: dict[str, int],
    rows: list[dict[str, Any]],
) -> pathlib.Path:
    with open(destination, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig in contigs:
            length = contig_lengths.get(contig)
            if length is not None:
                out.write(f"##contig=<ID={contig},length={length}>\n")
            else:
                out.write(f"##contig=<ID={contig}>\n")
        out.write("##INFO=<ID=PGS_IDS,Number=.,Type=String,Description=\"PGS IDs contributing to this scaffold site\">\n")
        out.write("##INFO=<ID=SITE_SOURCE,Number=.,Type=String,Description=\"Origin of the scaffold site (reference_db or score_extra)\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            info_parts = []
            if row.get("pgs_ids"):
                info_parts.append("PGS_IDS=" + ",".join(sorted(row["pgs_ids"])))
            if row.get("site_source"):
                info_parts.append("SITE_SOURCE=" + ",".join(sorted(row["site_source"])))
            out.write(
                f"{row['chrom']}\t{row['pos']}\t.\t{row['ref']}\t{','.join(row['alts']) or '.'}\t.\tPASS\t"
                + (";".join(info_parts) if info_parts else ".")
                + "\n"
            )
    return destination


def compress_vcf(path: pathlib.Path, *, log_path: pathlib.Path) -> pathlib.Path:
    output = path.with_suffix(path.suffix + ".gz")
    if completed_vcf(output):
        return output
    remove_if_exists(output, f"{output}.tbi")
    run_command(["bcftools", "view", "-Oz", "-o", str(output), str(path)], check=True, log_path=log_path)
    run_command(["tabix", "-f", "-p", "vcf", str(output)], check=True, log_path=log_path.with_name(log_path.stem + ".index.log"))
    record_vcf_checkpoint(output)
    return output


def merge_scaffold_contig_vcfs(
    *,
    backbone_vcf: pathlib.Path | None,
    extras_vcf: pathlib.Path | None,
    output_vcf: pathlib.Path,
    log_path: pathlib.Path,
) -> pathlib.Path:
    if output_vcf.exists() and completed_vcf(output_vcf):
        return output_vcf
    remove_if_exists(output_vcf, f"{output_vcf}.tbi")
    if backbone_vcf and extras_vcf:
        command = (
            "set -euo pipefail; "
            f"bcftools concat -a -Ou {shlex.quote(str(backbone_vcf))} {shlex.quote(str(extras_vcf))} "
            f"| bcftools sort -Ou "
            f"| bcftools norm --rm-dup exact -Oz -o {shlex.quote(str(output_vcf))}; "
            f"tabix -f -p vcf {shlex.quote(str(output_vcf))}"
        )
        run_command(["bash", "-c", command], check=True, log_path=log_path)
        record_vcf_checkpoint(output_vcf)
        return output_vcf
    source = backbone_vcf or extras_vcf
    if source is None:
        raise RuntimeError("Cannot merge scaffold contig VCFs without an input source.")
    run_command(["bcftools", "view", "-Oz", "-o", str(output_vcf), str(source)], check=True, log_path=log_path)
    run_command(["tabix", "-f", "-p", "vcf", str(output_vcf)], check=True, log_path=log_path.with_name(log_path.stem + ".index.log"))
    record_vcf_checkpoint(output_vcf)
    return output_vcf


def resolve_imputation_panel_manifest_path(context: dict[str, Any]) -> pathlib.Path:
    manifest_path = pathlib.Path(context["manifest_path"])
    panel_ref = ((context["manifest"]["analysis"]["prs"].get("imputation") or {}).get("reference_panel") or "").strip()
    if not panel_ref:
        raise RuntimeError("analysis.prs.imputation.reference_panel is required when PRS imputation is enabled.")

    if is_manifest_path(panel_ref):
        candidate = pathlib.Path(normalize_path(panel_ref, base_dir=manifest_path.parent))
        if candidate.exists():
            return candidate
        raise RuntimeError(f"Imputation panel manifest not found: {candidate}")

    candidates = [
        manifest_path.parent / "imputation_panels" / f"{panel_ref}.yaml",
        REPO_ROOT / "manifests" / "local" / "imputation_panels" / f"{panel_ref}.yaml",
        REPO_ROOT / "manifests" / "examples" / "imputation_panels" / f"{panel_ref}.yaml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise RuntimeError(f"No imputation panel manifest found for reference_panel={panel_ref}")


def infer_panel_index_path(path: pathlib.Path) -> pathlib.Path | None:
    candidates = [pathlib.Path(str(path) + ".tbi"), pathlib.Path(str(path) + ".csi")]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def ensure_staged_path(source: pathlib.Path, destination_dir: pathlib.Path, linux_cache_root: pathlib.Path) -> pathlib.Path:
    if needs_staging(source, linux_cache_root):
        fingerprint = sample_file_fingerprint(source)
        return pathlib.Path(stage_file(source, destination_dir, fingerprint=fingerprint))
    return source


def load_imputation_panel(context: dict[str, Any]) -> dict[str, Any]:
    manifest_path = resolve_imputation_panel_manifest_path(context)
    panel = load_yaml(manifest_path) or {}
    if not isinstance(panel, dict):
        raise RuntimeError(f"Imputation panel manifest must be a YAML object: {manifest_path}")

    build = str(panel.get("build") or "").strip()
    if build and build != context["manifest"]["reference"]["build"]:
        raise RuntimeError(
            f"Imputation panel build mismatch: panel={build} manifest={context['manifest']['reference']['build']}"
        )

    chrom_vcf_pattern = panel.get("chrom_vcf_pattern")
    if not chrom_vcf_pattern:
        raise RuntimeError(f"Imputation panel manifest is missing chrom_vcf_pattern: {manifest_path}")

    if "contig_prefix" in panel:
        panel_prefix = str(panel.get("contig_prefix") or "")
    else:
        panel_prefix = str(context["contig_prefix"] or "")
    chromosomes = [str(item) for item in (panel.get("chromosomes") or [str(i) for i in range(1, 23)] + ["X"])]
    base_chromosomes = sorted({base_contig_name(item, panel_prefix) for item in chromosomes})

    normalized_pattern = normalize_path(chrom_vcf_pattern, base_dir=manifest_path.parent)
    map_pattern = panel.get("map_pattern")
    normalized_map = normalize_path(map_pattern, base_dir=manifest_path.parent) if map_pattern else None

    return {
        "panel_id": str(panel.get("panel_id") or manifest_path.stem),
        "manifest_path": str(manifest_path),
        "manifest_fingerprint_id": sample_file_fingerprint(manifest_path)["id"],
        "build": context["manifest"]["reference"]["build"],
        "contig_prefix": panel_prefix,
        "chrom_vcf_pattern": normalized_pattern,
        "map_pattern": normalized_map,
        "chromosomes": chromosomes,
        "base_chromosomes": base_chromosomes,
    }


def resolve_panel_resources_for_contig(panel: dict[str, Any], *, contig: str, target_prefix: str) -> dict[str, Any] | None:
    base_chrom = base_contig_name(contig, target_prefix)
    if base_chrom not in set(panel["base_chromosomes"]):
        return None

    panel_chrom = normalize_contig(base_chrom, panel["contig_prefix"])
    format_vars = {
        "chrom": panel_chrom,
        "base_chrom": base_chrom,
        "chrom_no_prefix": base_chrom,
    }
    vcf_path = pathlib.Path(panel["chrom_vcf_pattern"].format(**format_vars))
    map_pattern = panel.get("map_pattern")
    map_path = pathlib.Path(map_pattern.format(**format_vars)) if map_pattern else None
    return {
        "base_chrom": base_chrom,
        "panel_chrom": panel_chrom,
        "vcf": vcf_path,
        "map": map_path,
    }


def ensure_beagle_jar(tool_cache_root: pathlib.Path) -> pathlib.Path:
    beagle_root = ensure_dir(tool_cache_root / "beagle")
    jar_path = beagle_root / pathlib.Path(BEAGLE_URL).name
    if jar_path.exists():
        return jar_path
    download_url_to_path(BEAGLE_URL, jar_path)
    return jar_path


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
    payload = json.loads(fetch_url_text(f"https://www.pgscatalog.org/rest/score/{pgs_id}"))
    return payload["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]


def resolve_cached_scorefile(score_cache_root: pathlib.Path, pgs_id: str) -> pathlib.Path | None:
    exact = score_cache_root / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if exact.exists() and is_valid_gzip_file(exact):
        return exact
    remove_if_exists(exact)
    matches = sorted(score_cache_root.glob(f"{pgs_id}_*.txt.gz"))
    for candidate in matches:
        if is_valid_gzip_file(candidate):
            return candidate
        remove_if_exists(candidate)
    return None


def build_union_sites(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    reference_fasta: pathlib.Path,
    reference_fai: pathlib.Path,
) -> dict[str, Any]:
    pgs_ids = sorted(context["manifest"]["analysis"]["prs"]["pgs_ids"])
    prefix = context["contig_prefix"]
    reference_db = ensure_reference_db_assets(context)
    backbone_base_chromosomes = allowed_prs_base_chromosomes(context)
    score_hash = hash_strings(
        [
            f"scaffold-v{UNION_SITE_VERSION}",
            GATK_IMAGE,
            reference_db["reference_db_id"],
            reference_db["archive_fingerprint"],
            sample_file_fingerprint(reference_fasta)["id"],
            str(context["manifest"]["analysis"]["prs"].get("min_depth", 10)),
            str(context["manifest"]["analysis"]["prs"].get("min_gq", 20)),
            *backbone_base_chromosomes,
            *pgs_ids,
        ]
    )[:16]
    targets_root = ensure_dir(stage_root / score_hash / "targets")
    union_manifest_path = targets_root / "union_manifest.json"
    reuse_key = {
        "union_site_version": UNION_SITE_VERSION,
        "pgs_ids": pgs_ids,
        "reference_db_id": reference_db["reference_db_id"],
        "reference_db_archive_fingerprint": reference_db["archive_fingerprint"],
        "backbone_base_chromosomes": backbone_base_chromosomes,
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
    contig_lengths = {contig: meta[0] for contig, meta in fasta.index.items()}
    allowed_contigs = allowed_prs_contigs(context)
    extras_by_contig: dict[str, dict[tuple[str, int, str], dict[str, object]]] = {}
    per_score_counts: dict[str, int] = {}
    score_metadata: dict[str, dict[str, str]] = {}
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
    backbone_sites = 0
    score_extra_sites = 0
    scaffold_contigs: list[str] = []
    backbone_contig_vcfs: list[pathlib.Path] = []
    extra_contig_vcfs: list[pathlib.Path] = []

    for base_chrom in backbone_base_chromosomes:
        pvar_path = export_reference_db_pvar_for_contig(
            context=context,
            reference_db=reference_db,
            targets_root=targets_root,
            base_chrom=base_chrom,
        )
        target_contig = normalize_contig(base_chrom, prefix)
        rows = []
        for record in iter_reference_db_sites(pvar_path, target_prefix=prefix):
            if record["chrom"] not in allowed_contigs:
                continue
            rows.append(
                {
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "ref": record["ref"],
                    "alts": record["alts"],
                    "pgs_ids": set(),
                    "site_source": {"reference_db"},
                }
            )
        if not rows:
            continue
        backbone_sites += len(rows)
        scaffold_contigs.append(target_contig)
        backbone_plain = targets_root / f"{target_contig}.backbone.vcf"
        write_scaffold_plain_vcf(
            destination=backbone_plain,
            contigs=[target_contig],
            contig_lengths=contig_lengths,
            rows=rows,
        )
        backbone_contig_vcfs.append(
            compress_vcf(backbone_plain, log_path=targets_root / f"{target_contig}.backbone.compress.log")
        )

    for pgs_id in pgs_ids:
        score_path = resolve_cached_scorefile(score_cache_root, pgs_id)
        if score_path is None:
            url = fetch_scorefile_url(pgs_id)
            score_path = score_cache_root / pathlib.Path(url).name
            download_url_to_path(url, score_path, validator=is_valid_gzip_file)

        per_score_counts[pgs_id] = 0
        with gzip.open(score_path, "rt", encoding="utf-8", errors="replace") as handle:
            header = None
            for raw in handle:
                if raw.startswith("#"):
                    key, separator, value = raw[1:].strip().partition("=")
                    if separator and key == "trait_reported":
                        score_metadata.setdefault(pgs_id, {})["trait"] = value
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
                if chrom not in allowed_contigs:
                    skipped += 1
                    skipped_reasons["non_scaffold_contig"] = skipped_reasons.get("non_scaffold_contig", 0) + 1
                    increment_score_skip(per_score_audit, pgs_id, "non_scaffold_contig")
                    continue
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
                if not resolved["alts"]:
                    skipped += 1
                    skipped_reasons["effect_only_without_backbone_site"] = skipped_reasons.get("effect_only_without_backbone_site", 0) + 1
                    increment_score_skip(per_score_audit, pgs_id, "effect_only_without_backbone_site")
                    continue
                contig_map = extras_by_contig.setdefault(chrom, {})
                entry = contig_map.setdefault(key, {"ref": resolved["ref_allele"], "alts": set(), "pgs_ids": set(), "site_source": {"score_extra"}})
                new_site = len(entry["alts"]) == 0 and len(entry["pgs_ids"]) == 0
                entry["alts"].update(resolved["alts"])
                entry["pgs_ids"].add(pgs_id)
                if new_site:
                    score_extra_sites += 1
                per_score_counts[pgs_id] += 1
                accepted_rows_by_source[resolved["source"]] = accepted_rows_by_source.get(resolved["source"], 0) + 1
                record_score_accept(per_score_audit, pgs_id, resolved)
                if resolved["multiallelic"]:
                    multiallelic_rows += 1
                if resolved.get("strand_flipped"):
                    strand_flip_rows += 1

    for chrom, contig_map in extras_by_contig.items():
        rows = [
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref_allele,
                "alts": sorted(entry["alts"]),
                "pgs_ids": entry["pgs_ids"],
                "site_source": entry["site_source"],
            }
            for (_, pos, ref_allele), entry in sorted(contig_map.items(), key=lambda item: (item[0][1], item[0][2]))
        ]
        extras_plain = targets_root / f"{chrom}.score_extra.vcf"
        write_scaffold_plain_vcf(
            destination=extras_plain,
            contigs=[chrom],
            contig_lengths=contig_lengths,
            rows=rows,
        )
        extra_contig_vcfs.append(compress_vcf(extras_plain, log_path=targets_root / f"{chrom}.score_extra.compress.log"))
        if chrom not in scaffold_contigs:
            scaffold_contigs.append(chrom)

    per_contig_scaffold_vcfs: list[str] = []
    for chrom in sorted(set(scaffold_contigs), key=lambda item: contig_order.get(item, 10**9)):
        backbone_vcf = next((path for path in backbone_contig_vcfs if path.name.startswith(f"{chrom}.backbone")), None)
        extras_vcf = next((path for path in extra_contig_vcfs if path.name.startswith(f"{chrom}.score_extra")), None)
        contig_output = targets_root / f"{chrom}.scaffold.vcf.gz"
        merge_scaffold_contig_vcfs(
            backbone_vcf=backbone_vcf,
            extras_vcf=extras_vcf,
            output_vcf=contig_output,
            log_path=targets_root / f"{chrom}.scaffold.merge.log",
        )
        per_contig_scaffold_vcfs.append(str(contig_output))

    targets_vcf = targets_root / "pgs_scaffold_sites.vcf.gz"
    if not completed_vcf(targets_vcf):
        remove_if_exists(targets_vcf, f"{targets_vcf}.tbi")
        run_command(
            ["bcftools", "concat", "-a", "-Oz", "-o", str(targets_vcf), *per_contig_scaffold_vcfs],
            check=True,
            log_path=targets_root / "scaffold.concat.log",
        )
        run_command(["tabix", "-f", "-p", "vcf", str(targets_vcf)], check=True, log_path=targets_root / "scaffold.concat.index.log")
        record_vcf_checkpoint(targets_vcf)

    scaffold_site_count = int((run_command(["bcftools", "index", "-n", str(targets_vcf)], check=True).stdout or "0").strip())
    summary = {
        "reference_db": reference_db["reference_db_id"],
        "scaffold_sites": scaffold_site_count,
        "union_sites": scaffold_site_count,
        "backbone_sites": backbone_sites,
        "backbone_base_chromosomes": backbone_base_chromosomes,
        "score_extra_sites": score_extra_sites,
        "skipped_rows": skipped,
        "skipped_reasons": skipped_reasons,
        "per_score_counts": per_score_counts,
        "score_metadata": score_metadata,
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
            "reference_db_archive": reference_db["archive_path"],
            "reference_db_root": reference_db["extracted_root"],
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
    gvcf_root: pathlib.Path | None = None,
    min_depth: int = 10,
    min_gq: int = 20,
) -> dict[str, Any]:
    shard_dir = ensure_dir(stage_root / "dense" / "shards" / shard["name"])
    log_dir = ensure_dir(stage_root / "logs")
    shard_targets = shard_dir / f"{shard['name']}.scaffold.targets.vcf.gz"
    gvcf_dir = germline_gvcf_dir(
        gvcf_root or stage_root / "germline_gvcf", cram_path=cram_path,
        reference_fasta=reference_fasta, image=GATK_IMAGE, intervals=shard["intervals"],
    )
    gvcf = gvcf_dir / "germline.g.vcf.gz"
    dense_vcf = shard_dir / f"{sample_id}.{shard['name']}.scaffold.vcf.gz"
    qc_path = shard_dir / "genotype_qc.json"
    region_arg = ",".join(shard["intervals"])
    targets_receipt = shard_dir / "targets.completed.json"
    gvcf_receipt = gvcf_dir / "completed.json"
    dense_receipt = shard_dir / "scaffold.completed.json"
    target_artifacts = [shard_targets, pathlib.Path(str(shard_targets) + ".tbi")]
    gvcf_artifacts = [gvcf, pathlib.Path(str(gvcf) + ".tbi")]
    dense_artifacts = [dense_vcf, pathlib.Path(str(dense_vcf) + ".tbi"), qc_path]
    shard_targets_reused = checkpoint_complete(targets_receipt, target_artifacts)
    gvcf_reused = checkpoint_complete(gvcf_receipt, gvcf_artifacts)
    dense_vcf_reused = checkpoint_complete(dense_receipt, dense_artifacts)

    if not dense_vcf_reused:
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
            record_checkpoint(targets_receipt, target_artifacts)

        if not gvcf_reused:
            remove_if_exists(gvcf, f"{gvcf}.tbi")
            run_command(
                build_docker_run_command(
                    image=GATK_IMAGE,
                    mounts=[
                        (cram_path.parent, "/input", "ro"),
                        (reference_fasta.parent, "/ref", "ro"),
                        (gvcf_dir, "/gvcf", "rw"),
                    ],
                    inner_args=[
                        "gatk",
                        "--java-options",
                        f"-Xmx{java_memory}",
                        "HaplotypeCaller",
                        "-R", f"/ref/{reference_fasta.name}",
                        "-I", f"/input/{cram_path.name}",
                        "-O", f"/gvcf/{gvcf.name}",
                        "-ERC", "GVCF",
                        "--native-pair-hmm-threads", str(pairhmm_threads),
                        *[item for interval in shard["intervals"] for item in ("-L", interval)],
                    ],
                    docker_bin=docker_bin,
                ),
                check=True,
                log_path=log_dir / f"{shard['name']}.haplotypecaller.log",
            )
            record_checkpoint(gvcf_receipt, gvcf_artifacts)

        remove_if_exists(dense_vcf, f"{dense_vcf}.tbi")
        finalized = shard_dir / "finalized.vcf.gz"
        normalized = shard_dir / "finalized.norm.vcf.gz"
        normalized_targets = shard_dir / "targets.norm.vcf.gz"
        finalized_receipt = shard_dir / "finalized.completed.json"
        finalized_artifacts = [finalized, pathlib.Path(str(finalized) + ".tbi")]
        if not checkpoint_complete(finalized_receipt, finalized_artifacts):
            remove_if_exists(finalized, f"{finalized}.tbi", finalized_receipt)
            run_command(
                build_docker_run_command(
                    image=GATK_IMAGE,
                    mounts=[(reference_fasta.parent, "/ref", "ro"),
                            (gvcf_dir, "/gvcf", "ro"), (shard_dir, "/work", "rw")],
                    inner_args=["gatk", "--java-options", f"-Xmx{java_memory}", "GenotypeGVCFs",
                                "-R", f"/ref/{reference_fasta.name}", "-V", f"/gvcf/{gvcf.name}",
                                "-L", f"/work/{shard_targets.name}",
                                "--include-non-variant-sites", "true",
                                "-O", f"/work/{finalized.name}"],
                    docker_bin=docker_bin,
                ), check=True, log_path=log_dir / f"{shard['name']}.genotypegvcfs.log",
            )
            record_checkpoint(finalized_receipt, finalized_artifacts)
        for source, destination in ((finalized, normalized), (shard_targets, normalized_targets)):
            receipt = destination.with_name(destination.name + ".completed.json")
            if not checkpoint_complete(receipt, [destination]):
                remove_if_exists(destination, receipt)
                run_command(
                    ["bcftools", "norm", "--check-ref", "e", "-f", str(reference_fasta),
                     "--rm-dup", "exact", "-Oz", "-o", str(destination), str(source)],
                    check=True, log_path=log_dir / f"{shard['name']}.{destination.stem}.log",
                )
                record_checkpoint(receipt, [destination])
        from workflow.stages.prs_genotypes import materialize_scaffold
        plain = shard_dir / "scaffold.materialized.vcf"
        with gzip.open(normalized, "rt") as calls, gzip.open(normalized_targets, "rt") as targets, plain.open("w") as output:
            qc = materialize_scaffold(calls, targets, output, sample_id=sample_id,
                                      min_depth=min_depth, min_gq=min_gq)
        run_command(["bcftools", "sort", "-Oz", "-o", str(dense_vcf), str(plain)],
                    check=True, log_path=log_dir / f"{shard['name']}.scaffold.log")
        run_command(["tabix", "-f", "-p", "vcf", str(dense_vcf)], check=True,
                    log_path=log_dir / f"{shard['name']}.scaffold.index.log")
        write_json(qc_path, qc)
        record_checkpoint(dense_receipt, dense_artifacts)
        remove_if_exists(plain)
    return {
        "shard_name": shard["name"],
        "dense_vcf": str(dense_vcf),
        "gvcf": str(gvcf),
        "targets_vcf": str(shard_targets),
        "reused": dense_vcf_reused,
        "reused_targets_vcf": shard_targets_reused,
        "reused_gvcf": gvcf_reused,
        "reused_dense_vcf": dense_vcf_reused,
        "genotype_qc": load_json(qc_path),
    }


def concat_dense_vcf(stage_root: pathlib.Path, reference_fasta: pathlib.Path, sample_id: str, dense_shards: list[str]) -> str:
    final_vcf = stage_root / "dense" / f"{sample_id}.prs_scaffold.vcf.gz"
    if not completed_vcf(final_vcf):
        remove_if_exists(final_vcf, f"{final_vcf}.tbi")
        pipeline = (
            "set -euo pipefail; "
            f"bcftools concat -a -Ou {shlex.join(dense_shards)} "
            f"| bcftools norm -f {shlex.quote(str(reference_fasta))} --rm-dup exact -Oz -o {shlex.quote(str(final_vcf))}; "
            f"tabix -f -p vcf {shlex.quote(str(final_vcf))}"
        )
        run_command(["bash", "-c", pipeline], check=True, log_path=stage_root / "logs" / "dense.concat.log")
        record_vcf_checkpoint(final_vcf)
    return str(final_vcf)



def filter_imputed_expression(info_threshold: float) -> str:
    return f"(INFO/IMP=1 && INFO/DR2>={info_threshold}) || INFO/IMP=0"


def beagle_failed_due_to_no_common_markers(error: subprocess.CalledProcessError) -> bool:
    combined = "\n".join(part for part in [error.output, error.stderr] if part)
    return "no markers in common" in combined.lower()


def build_imputation_shards(context: dict[str, Any], dense_vcf: pathlib.Path, panel: dict[str, Any]) -> list[dict[str, Any]]:
    dense_contigs = [line.strip() for line in run_command(["tabix", "-l", str(dense_vcf)], check=True).stdout.splitlines() if line.strip()]
    shards = []
    for shard in context["shards"]["prs"]:
        shard_contigs = [contig for contig in shard["intervals"] if contig in dense_contigs]
        if not shard_contigs:
            continue
        if any(base_contig_name(contig, context["contig_prefix"]) not in {str(i) for i in range(1, 23)} for contig in shard_contigs):
            raise ValueError("Local PRS imputation currently supports diploid autosomes only; import validated imputed sex-chromosome genotypes separately")
        panel_resources = [resolve_panel_resources_for_contig(panel, contig=contig, target_prefix=context["contig_prefix"]) for contig in shard_contigs]
        imputation_mode = "impute" if all(panel_resources) else "passthrough"
        if imputation_mode == "impute":
            base_chroms = sorted({resource["base_chrom"] for resource in panel_resources if resource})
            if len(base_chroms) != 1:
                raise RuntimeError(f"Imputation shard {shard['name']} spans multiple panel chromosomes: {base_chroms}")
            panel_resource = panel_resources[0]
        else:
            panel_resource = None
        shards.append({**shard, "intervals": shard_contigs, "imputation_mode": imputation_mode, "panel_resource": panel_resource})
    return shards


def stage_panel_resources(panel_resource: dict[str, Any], *, stage_root: pathlib.Path, linux_cache_root: pathlib.Path) -> dict[str, pathlib.Path]:
    panel_dir = ensure_dir(stage_root / "panel_cache" / panel_resource["base_chrom"])
    panel_vcf = pathlib.Path(panel_resource["vcf"])
    if not panel_vcf.exists():
        raise RuntimeError(f"Imputation panel VCF not found: {panel_vcf}")

    panel_index = infer_panel_index_path(panel_vcf)
    if panel_index is None:
        raise RuntimeError(f"Imputation panel VCF is missing a tabix/csi index: {panel_vcf}")

    staged_vcf = ensure_staged_path(panel_vcf, panel_dir, linux_cache_root)
    staged_index = ensure_staged_path(panel_index, panel_dir, linux_cache_root)

    staged_map = None
    map_path = panel_resource.get("map")
    if map_path is not None:
        map_target = pathlib.Path(map_path)
        if not map_target.exists():
            raise RuntimeError(f"Imputation genetic map not found: {map_target}")
        staged_map = ensure_staged_path(map_target, panel_dir, linux_cache_root)

    return {
        "panel_vcf": staged_vcf,
        "panel_index": staged_index,
        "map": staged_map,
    }


def extract_dense_contig_vcf(
    *,
    dense_vcf: pathlib.Path,
    shard_dir: pathlib.Path,
    shard_name: str,
    region_arg: str,
) -> pathlib.Path:
    target_vcf = shard_dir / f"{shard_name}.target.vcf.gz"
    if completed_vcf(target_vcf):
        return target_vcf

    remove_if_exists(target_vcf, f"{target_vcf}.tbi")
    run_command(
        ["bcftools", "view", "-r", region_arg, "-Oz", "-o", str(target_vcf), str(dense_vcf)],
        check=True,
        log_path=shard_dir / f"{shard_name}.target.extract.log",
    )
    run_command(
        ["tabix", "-f", "-p", "vcf", str(target_vcf)],
        check=True,
        log_path=shard_dir / f"{shard_name}.target.extract.index.log",
    )
    record_vcf_checkpoint(target_vcf)
    return target_vcf


def rename_contig_for_panel(
    *,
    input_vcf: pathlib.Path,
    shard_dir: pathlib.Path,
    shard_name: str,
    source_contig: str,
    target_contig: str,
) -> pathlib.Path:
    if source_contig == target_contig:
        return input_vcf

    renamed_vcf = shard_dir / f"{shard_name}.panel_target.vcf.gz"
    if completed_vcf(renamed_vcf):
        return renamed_vcf

    rename_map = shard_dir / f"{shard_name}.rename.tsv"
    rename_map.write_text(f"{source_contig}\t{target_contig}\n", encoding="utf-8")
    header_result = run_command(
        ["bcftools", "view", "-h", str(input_vcf)],
        check=True,
        log_path=shard_dir / f"{shard_name}.rename.header.log", capture_output=True,
    )
    header_fixed_input = input_vcf
    updated_header = inject_missing_contig_headers(header_result.stdout, [source_contig])
    if updated_header != header_result.stdout:
        header_path = shard_dir / f"{shard_name}.rename.header.txt"
        header_path.write_text(updated_header, encoding="utf-8")
        header_fixed_input = shard_dir / f"{shard_name}.rename.headerfix.vcf.gz"
        remove_if_exists(header_fixed_input, f"{header_fixed_input}.tbi")
        run_command(
            ["bcftools", "reheader", "-h", str(header_path), "-o", str(header_fixed_input), str(input_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.rename.reheader.log",
        )
    remove_if_exists(renamed_vcf, f"{renamed_vcf}.tbi")
    run_command(
        ["bcftools", "annotate", "--rename-chrs", str(rename_map), "-Oz", "-o", str(renamed_vcf), str(header_fixed_input)],
        check=True,
        log_path=shard_dir / f"{shard_name}.rename.log",
    )
    run_command(
        ["tabix", "-f", "-p", "vcf", str(renamed_vcf)],
        check=True,
        log_path=shard_dir / f"{shard_name}.rename.index.log",
    )
    record_vcf_checkpoint(renamed_vcf)
    return renamed_vcf


def build_dense_passthrough_shard(
    *,
    dense_vcf: pathlib.Path,
    shard_dir: pathlib.Path,
    shard_name: str,
    sample_id: str,
    region_arg: str,
) -> tuple[pathlib.Path, bool]:
    passthrough_vcf = shard_dir / f"{sample_id}.{shard_name}.dense_passthrough.vcf.gz"
    passthrough_reused = completed_vcf(passthrough_vcf)
    if not passthrough_reused:
        remove_if_exists(passthrough_vcf, f"{passthrough_vcf}.tbi")
        run_command(
            ["bcftools", "view", "-r", region_arg, "-Oz", "-o", str(passthrough_vcf), str(dense_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.passthrough.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(passthrough_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.passthrough.index.log",
        )
        record_vcf_checkpoint(passthrough_vcf)
    return passthrough_vcf, passthrough_reused


def target_has_imputable_markers(target_vcf: pathlib.Path) -> bool:
    opener = gzip.open if target_vcf.name.endswith((".gz", ".bgz")) else open
    with opener(target_vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) >= 5 and any(re.fullmatch(r"[ACGT]+", alt) for alt in fields[4].split(",")):
                return True
    return False


def finalize_imputed_vcf(
    *,
    input_vcf: pathlib.Path,
    shard_dir: pathlib.Path,
    shard_name: str,
    reference_fasta: pathlib.Path,
    info_threshold: float,
    output_name: str,
) -> pathlib.Path:
    norm_vcf = shard_dir / f"{shard_name}.imputed.norm.vcf.gz"
    if not completed_vcf(norm_vcf):
        remove_if_exists(norm_vcf, f"{norm_vcf}.tbi")
        run_command(
            ["bcftools", "norm", "-f", str(reference_fasta), "-m", "-both", "-Oz", "-o", str(norm_vcf), str(input_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.norm.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(norm_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.norm.index.log",
        )
        record_vcf_checkpoint(norm_vcf)

    final_vcf = shard_dir / output_name
    if not completed_vcf(final_vcf):
        remove_if_exists(final_vcf, f"{final_vcf}.tbi")
        run_command(
            ["bcftools", "view", "-i", filter_imputed_expression(info_threshold), "-Oz", "-o", str(final_vcf), str(norm_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.filter.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(final_vcf)],
            check=True,
            log_path=shard_dir / f"{shard_name}.filter.index.log",
        )
        record_vcf_checkpoint(final_vcf)
    return final_vcf


def build_imputed_prs_shard(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    dense_vcf: pathlib.Path,
    reference_fasta: pathlib.Path,
    shard: dict[str, Any],
    panel: dict[str, Any],
    beagle_jar: pathlib.Path,
    java_memory: str,
    beagle_threads: int,
    info_threshold: float,
) -> dict[str, Any]:
    shard_dir = ensure_dir(stage_root / "shards" / shard["name"])
    region_arg = ",".join(shard["intervals"])

    if shard["imputation_mode"] == "passthrough":
        passthrough_vcf, passthrough_reused = build_dense_passthrough_shard(
            dense_vcf=dense_vcf,
            shard_dir=shard_dir,
            shard_name=shard["name"],
            sample_id=context["sample_id"],
            region_arg=region_arg,
        )
        return {
            "shard_name": shard["name"],
            "mode": "passthrough",
            "output_vcf": str(passthrough_vcf),
            "reused": passthrough_reused,
            "panel_chrom": None,
            "passthrough_reason": "panel_unavailable",
        }

    final_vcf = shard_dir / f"{context['sample_id']}.{shard['name']}.imputed.filtered.vcf.gz"
    if completed_vcf(final_vcf):
        return {
            "shard_name": shard["name"],
            "mode": "imputed",
            "output_vcf": str(final_vcf),
            "reused": True,
            "panel_chrom": shard["panel_resource"]["panel_chrom"],
        }

    target_vcf = extract_dense_contig_vcf(
        dense_vcf=dense_vcf,
        shard_dir=shard_dir,
        shard_name=shard["name"],
        region_arg=region_arg,
    )
    target_input = rename_contig_for_panel(
        input_vcf=target_vcf,
        shard_dir=shard_dir,
        shard_name=shard["name"],
        source_contig=shard["intervals"][0],
        target_contig=shard["panel_resource"]["panel_chrom"],
    )
    if not target_has_imputable_markers(target_input):
        passthrough_vcf, passthrough_reused = build_dense_passthrough_shard(
            dense_vcf=dense_vcf,
            shard_dir=shard_dir,
            shard_name=shard["name"],
            sample_id=context["sample_id"],
            region_arg=region_arg,
        )
        return {
            "shard_name": shard["name"],
            "mode": "passthrough",
            "output_vcf": str(passthrough_vcf),
            "reused": passthrough_reused,
            "panel_chrom": shard["panel_resource"]["panel_chrom"],
            "passthrough_reason": "no_imputable_markers",
        }
    panel_resources = stage_panel_resources(
        shard["panel_resource"],
        stage_root=stage_root,
        linux_cache_root=pathlib.Path(context["runtime"]["linux_cache_root"]),
    )
    beagle_out_prefix = shard_dir / f"{context['sample_id']}.{shard['name']}.beagle"
    beagle_vcf = pathlib.Path(str(beagle_out_prefix) + ".vcf.gz")
    beagle_log = shard_dir / f"{shard['name']}.beagle.log"
    beagle_receipt = shard_dir / "beagle.completed.json"
    if not checkpoint_complete(beagle_receipt, [beagle_vcf]):
        remove_if_exists(beagle_vcf, beagle_log, beagle_receipt)
        args = [
            "java",
            f"-Xmx{java_memory}",
            "-jar",
            str(beagle_jar),
            f"gt={target_input}",
            f"ref={panel_resources['panel_vcf']}",
            f"out={beagle_out_prefix}",
            f"chrom={shard['panel_resource']['panel_chrom']}",
            f"nthreads={beagle_threads}",
        ]
        if panel_resources["map"] is not None:
            args.append(f"map={panel_resources['map']}")
        try:
            run_command(
                args,
                check=True,
                log_path=beagle_log,
            )
            record_checkpoint(beagle_receipt, [beagle_vcf])
        except subprocess.CalledProcessError as error:
            if not beagle_failed_due_to_no_common_markers(error):
                raise
            remove_if_exists(beagle_vcf)
            passthrough_vcf, passthrough_reused = build_dense_passthrough_shard(
                dense_vcf=dense_vcf,
                shard_dir=shard_dir,
                shard_name=shard["name"],
                sample_id=context["sample_id"],
                region_arg=region_arg,
            )
            return {
                "shard_name": shard["name"],
                "mode": "passthrough",
                "output_vcf": str(passthrough_vcf),
                "reused": passthrough_reused,
                "panel_chrom": shard["panel_resource"]["panel_chrom"],
                "passthrough_reason": "no_panel_overlap_window",
            }

    renamed_back = rename_contig_for_panel(
        input_vcf=beagle_vcf,
        shard_dir=shard_dir,
        shard_name=f"{shard['name']}.restore",
        source_contig=shard["panel_resource"]["panel_chrom"],
        target_contig=shard["intervals"][0],
    )
    final_vcf = finalize_imputed_vcf(
        input_vcf=renamed_back,
        shard_dir=shard_dir,
        shard_name=shard["name"],
        reference_fasta=reference_fasta,
        info_threshold=info_threshold,
        output_name=f"{context['sample_id']}.{shard['name']}.imputed.filtered.vcf.gz",
    )
    return {
        "shard_name": shard["name"],
        "mode": "imputed",
        "output_vcf": str(final_vcf),
        "reused": False,
        "panel_chrom": shard["panel_resource"]["panel_chrom"],
    }


def concat_imputed_vcf(stage_root: pathlib.Path, sample_id: str, shard_vcfs: list[str]) -> str:
    final_vcf = stage_root / f"{sample_id}.imputed_prs.vcf.gz"
    if not completed_vcf(final_vcf):
        remove_if_exists(final_vcf, f"{final_vcf}.tbi")
        pipeline = (
            "set -euo pipefail; "
            f"bcftools concat -a -Ou {shlex.join(shard_vcfs)} "
            f"| bcftools sort -Oz -o {shlex.quote(str(final_vcf))}; "
            f"tabix -f -p vcf {shlex.quote(str(final_vcf))}"
        )
        run_command(["bash", "-c", pipeline], check=True, log_path=stage_root / "concat.log")
        record_vcf_checkpoint(final_vcf)
    return str(final_vcf)



def load_existing_imputation_stage(stage_manifest_path: pathlib.Path, cache_key: str) -> dict[str, Any] | None:
    if not stage_manifest_path.exists():
        return None
    payload = load_json(stage_manifest_path)
    if payload.get("cache_key") != cache_key or payload.get("status") != "ok":
        return None
    artifacts = payload.get("artifacts") or {}
    if not artifacts.get("imputed_vcf") or not completed_vcf(pathlib.Path(artifacts["imputed_vcf"])):
        return None
    payload["reused"] = True
    return payload


def build_or_reuse_imputed_vcf(
    *,
    context: dict[str, Any],
    stage_root: pathlib.Path,
    dense_info: dict[str, Any],
    score_hash: str,
    reference_fasta: pathlib.Path,
    update_status: Callable[[dict[str, Any], dict[str, Any]], None] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    imputation_cfg = context["manifest"]["analysis"]["prs"].get("imputation") or {}
    panel = load_imputation_panel(context)
    dense_vcf = pathlib.Path(dense_info["artifacts"]["dense_vcf"])
    shards = build_imputation_shards(context, dense_vcf, panel)
    resource_fingerprints = []
    for shard in shards:
        resource = shard.get("panel_resource") or {}
        for key in ("vcf", "map"):
            if resource.get(key):
                resource_fingerprints.append(sample_file_fingerprint(resource[key])["id"])
    cache_key = hash_strings(
        [
            f"v{STAGE_VERSION}",
            f"imputation-v{IMPUTATION_STAGE_VERSION}",
            context["input_fingerprint"]["id"],
            score_hash,
            panel["panel_id"],
            panel["manifest_path"],
            panel["manifest_fingerprint_id"],
            sample_file_fingerprint(dense_vcf)["id"],
            sample_file_fingerprint(reference_fasta)["id"],
            BEAGLE_URL,
            *resource_fingerprints,
            str(imputation_cfg.get("info_threshold")),
            str(imputation_cfg.get("backend")),
        ]
    )
    imputation_root = ensure_dir(stage_root / score_hash / "imputation" / cache_key)
    stage_manifest_path = imputation_root / "stage_manifest.json"
    existing = load_existing_imputation_stage(stage_manifest_path, cache_key)
    if existing is not None:
        return existing, {"phase": "scoring_reused_imputed_vcf", "score_hash": score_hash}

    beagle_jar = ensure_beagle_jar(pathlib.Path(context["runtime"]["tool_cache_root"]))
    worker_plan = plan_worker_pool(
        stage_name="prs_imputation",
        shard_count=len(shards),
        max_cpus=int(context["manifest"]["runtime"]["max_cpus"]),
        max_memory=context["manifest"]["runtime"]["max_memory"],
        requested_workers=context["manifest"]["analysis"]["prs"].get("workers"),
        requested_pairhmm_threads=context["manifest"]["analysis"]["prs"].get("pairhmm_threads"),
        default_worker_cap=6,
        default_pairhmm_cap=3,
        min_memory_per_worker_gb=2.0,
        free_gb=context["disk"].get("free_gb"),
        min_free_gb=context["disk"].get("min_free_gb"),
    )
    progress = init_shard_progress(
        phase="imputing_dense_vcf",
        total_shards=len(shards),
        score_hash=score_hash,
    )
    running_payload = stage_payload(
        "prs_imputed_genotypes",
        "running",
        cache_key=cache_key,
        score_hash=score_hash,
        worker_plan=worker_plan,
        summary={
            "backend": imputation_cfg.get("backend"),
            "panel_id": panel["panel_id"],
            "info_threshold": imputation_cfg.get("info_threshold"),
            "target_shards": len(shards),
        },
        artifacts={
            "panel_manifest": panel["manifest_path"],
            "dense_vcf": dense_info["artifacts"]["dense_vcf"],
            "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"],
        },
    )
    if update_status is not None:
        update_status(progress, running_payload)

    tasks = [
        (
            lambda shard=shard: build_imputed_prs_shard(
                context=context,
                stage_root=imputation_root,
                dense_vcf=dense_vcf,
                reference_fasta=reference_fasta,
                shard=shard,
                panel=panel,
                beagle_jar=beagle_jar,
                java_memory=worker_plan["java_memory"],
                beagle_threads=worker_plan["pairhmm_threads"],
                info_threshold=float(imputation_cfg.get("info_threshold", 0.8)),
            )
        )
        for shard in shards
    ]
    shard_names = [shard["name"] for shard in shards]

    def on_complete(index: int, result: Any) -> None:
        record_shard_completion(progress, shard_name=result.get("shard_name") or shard_names[index], reused=bool(result.get("reused")))
        if update_status is not None:
            update_status(progress, running_payload)

    shard_results = run_tasks_in_parallel(tasks, max_workers=worker_plan["workers"], on_complete=on_complete)
    imputed_outputs = [result["output_vcf"] for result in shard_results]
    progress = {
        **progress,
        "phase": "concatenating_imputed_vcf",
        "reused_imputed_shards": sum(1 for item in shard_results if item.get("reused")),
        "imputed_shards": sum(1 for item in shard_results if item.get("mode") == "imputed"),
        "passthrough_shards": sum(1 for item in shard_results if item.get("mode") == "passthrough"),
    }
    if update_status is not None:
        update_status(progress, running_payload)

    panel_vcf = concat_imputed_vcf(imputation_root, context["sample_id"], imputed_outputs)
    # Reference-panel imputation must not erase direct-only score sites or replace
    # confidently observed genotypes. Supply DS for retained direct genotypes too.
    from workflow.stages.prs_genotypes import augment_scaffold
    augmented_plain = imputation_root / "augmented.vcf"
    augmented_vcf = imputation_root / f"{context['sample_id']}.augmented_prs.vcf.gz"
    with gzip.open(dense_vcf, "rt") as direct, gzip.open(panel_vcf, "rt") as imputed, augmented_plain.open("w") as output:
        genotype_qc = augment_scaffold(direct, imputed, output, sample_id=context["sample_id"])
    run_command(["bcftools", "sort", "-Oz", "-o", str(augmented_vcf), str(augmented_plain)],
                check=True, log_path=imputation_root / "augment.log")
    run_command(["tabix", "-f", "-p", "vcf", str(augmented_vcf)], check=True,
                log_path=imputation_root / "augment.index.log")
    record_vcf_checkpoint(augmented_vcf)
    remove_if_exists(augmented_plain)
    final_vcf = str(augmented_vcf)
    payload = stage_payload(
        "prs_imputed_genotypes",
        "ok",
        cache_key=cache_key,
        reused=False,
        score_hash=score_hash,
        worker_plan=worker_plan,
        summary={
            "backend": imputation_cfg.get("backend"),
            "panel_id": panel["panel_id"],
            "info_threshold": imputation_cfg.get("info_threshold"),
            "target_shards": len(shards),
            "imputed_shards": progress["imputed_shards"],
            "passthrough_shards": progress["passthrough_shards"],
            "genotype_qc": genotype_qc,
        },
        artifacts={
            "panel_manifest": panel["manifest_path"],
            "dense_vcf": dense_info["artifacts"]["dense_vcf"],
            "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"],
            "imputed_vcf": final_vcf,
            "imputed_vcf_index": f"{final_vcf}.tbi",
            "beagle_jar": str(beagle_jar),
        },
    )
    write_json(stage_manifest_path, payload)
    progress = {**progress, "phase": "scoring_imputed_vcf"}
    if update_status is not None:
        update_status(progress, payload)
    return payload, progress


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


def scoring_input_identity(context: dict[str, Any], input_vcf: pathlib.Path, input_kind: str, dosage_field: str | None) -> dict[str, Any]:
    return {
        "version": STAGE_VERSION,
        "input": sample_file_fingerprint(input_vcf)["id"],
        "input_kind": input_kind,
        "sample_id": context["sample_id"],
        "sample_sex": context["manifest"]["analysis"]["prs"].get("sample_sex"),
        "dosage_field": dosage_field,
        "plink_image": PLINK2_IMAGE,
    }


def build_pgsc_input_pfile(
    *, context: dict[str, Any], stage_root: pathlib.Path, input_vcf: pathlib.Path,
    input_kind: str, dosage_field: str | None = None,
) -> dict[str, Any]:
    identity = scoring_input_identity(context, input_vcf, input_kind, dosage_field)
    cache_key = hash_strings([json.dumps(identity, sort_keys=True)])
    prep_root = ensure_dir(stage_root / "compatible" / cache_key)
    output_prefix = prep_root / f"GRCh38_{safe_sampleset_name(context['sample_id'])}_ALL"
    provenance = prep_root / "conversion.json"
    sex_code = plink_sex_code(identity["sample_sex"])
    chrom_filter = "1-22,X,Y,XY" if sex_code else "1-22"
    max_cpus = max(int(context["manifest"]["runtime"]["max_cpus"]), 1)
    log_path = prep_root / "plink2.make_pfile.log"
    normalized_vcf = prep_root / "callable.vcf.gz"
    reused = complete_pfile(output_prefix) and provenance.is_file() and load_json(provenance) == identity
    if not reused:
        docker_bin = context["host_checks"]["docker"] or "docker"
        ensure_docker_ready(docker_bin)
        remove_if_exists(provenance)
        sample_names = run_command(["bcftools", "query", "-l", str(input_vcf)], check=True).stdout.splitlines()
        if len(sample_names) != 1:
            raise ValueError("PRS scoring requires exactly one sample in the input VCF")
        # A missing target row must not count toward site overlap or be mean-imputed
        # by PLINK. DS inputs require an actual dosage; direct inputs require GT.
        if dosage_field and not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", dosage_field):
            raise ValueError("Invalid VCF dosage field")
        expression = f'FMT/{dosage_field}!="."' if dosage_field else 'GT!="mis"'
        filtered = prep_root / "callable.original_name.vcf.gz"
        run_command(["bcftools", "view", "-i", expression, "-f", ".,PASS", "-Oz", "-o", str(filtered), str(input_vcf)],
                    check=True, log_path=prep_root / "filter_callable.log")
        sample_file = prep_root / "sample_name.txt"
        sample_file.write_text(context["sample_id"] + "\n", encoding="utf-8")
        run_command(["bcftools", "reheader", "-s", str(sample_file), "-o", str(normalized_vcf), str(filtered)],
                    check=True, log_path=prep_root / "reheader.log")
        args = ["plink2", "--threads", str(max(1, min(4, max_cpus))),
                "--memory", str(parse_memory_to_mb(context["manifest"]["runtime"]["max_memory"])),
                "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "512",
                "--max-alleles", "2", "--vcf", f"/work/{normalized_vcf.name}"]
        if dosage_field:
            args.append(f"dosage={dosage_field}")
        if sex_code:
            sex_file = prep_root / "sample_sex.tsv"
            sex_file.write_text(f"#FID\tIID\tSEX\n0\t{context['sample_id']}\t{sex_code}\n", encoding="utf-8")
            args.extend(["--update-sex", f"/work/{sex_file.name}"])
        args.extend(["--allow-extra-chr", "--chr", chrom_filter,
                     "--make-pgen", "vzs", "pvar-cols=-xheader,-maybequal,-maybefilter,-maybeinfo,-maybecm",
                     "--out", f"/work/{output_prefix.name}"])
        run_command(build_docker_run_command(image=PLINK2_IMAGE, mounts={prep_root: ("/work", "rw")},
                                             inner_args=args, docker_bin=docker_bin),
                    check=True, log_path=log_path)
        if not complete_pfile(output_prefix):
            raise RuntimeError("PLINK conversion did not produce a complete scoring input")
        write_json(provenance, identity)
        remove_if_exists(filtered)
    return {
        "format": "pfile", "path_prefix": str(output_prefix), "cache_key": cache_key,
        "input_kind": input_kind, "dosage_field": dosage_field, "reused": reused,
        "sample_sex": identity["sample_sex"], "sex_chromosomes_allowed": bool(sex_code),
        "renamed_chr_prefix": None,
        "artifacts": {"pgen": str(output_prefix) + ".pgen", "psam": str(output_prefix) + ".psam",
                      "pvar": str(output_prefix) + ".pvar.zst", "log": str(log_path),
                      "provenance": str(provenance), "normalized_vcf": str(normalized_vcf)},
    }


def build_pgsc_calc_samplesheet(
    *,
    sampleset: str,
    scoring_input_info: dict[str, Any],
    dosage_field: str | None = None,
) -> str:
    header = ["sampleset", "path_prefix", "chrom", "format"]
    row = [sampleset, scoring_input_info["path_prefix"], "", scoring_input_info["format"]]
    if scoring_input_info["format"] == "vcf" and dosage_field:
        header.append("vcf_genotype_field")
        row.append(dosage_field)
    return ",".join(header) + "\n" + ",".join(row) + "\n"


def build_pgsc_calc_resource_config(*, max_memory: str, max_cpus: int = 1, match_cpus: int = 1, score_cpus: int = 1) -> str:
    return "\n".join(
        [
            "executor {",
            f"  cpus = {max_cpus}",
            f"  memory = '{max_memory}'",
            "}",
            "process {",
            "  withName: 'PGSCATALOG_PGSCCALC:PGSCCALC:MATCH:MATCH_VARIANTS' {",
            f"    cpus = {match_cpus}",
            f"    memory = '{max_memory}'",
            "    time = '24.h'",
            "  }",
            "  withName: 'PGSCATALOG_PGSCCALC:PGSCCALC:APPLY_SCORE:PLINK2_SCORE' {",
            f"    cpus = {score_cpus}",
            f"    memory = '{max_memory}'",
            "    time = '24.h'",
            "  }",
            "  withName: 'PGSCATALOG_PGSCCALC:PGSCCALC:APPLY_SCORE:SCORE_AGGREGATE' {",
            f"    cpus = {score_cpus}",
            f"    memory = '{max_memory}'",
            "    time = '24.h'",
            "  }",
            "}",
            "",
        ]
    )


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
    outdir = stage_root / "pgsc_calc"
    workdir = stage_root / "work"
    genotypes_cache = ensure_dir(stage_root / "genotypes_cache")
    log_path = stage_root / "pgsc_calc.log"
    nested_config = stage_root / "pgsc_calc.resources.config"
    # Stable per-cache ownership survives resume without changing task configuration.
    docker_label = "genomics.scoring=" + hash_strings([str(stage_root.resolve())])[:32]
    scoring_input_info = build_pgsc_input_pfile(
        context=context,
        stage_root=stage_root,
        input_vcf=scoring_vcf,
        input_kind=input_kind,
        dosage_field=dosage_field,
    )

    samplesheet.write_text(
        build_pgsc_calc_samplesheet(
            sampleset=safe_sampleset_name(context["sample_id"]),
            scoring_input_info=scoring_input_info,
            dosage_field=dosage_field,
        ),
        encoding="utf-8",
    )

    ensure_dir(outdir)
    ensure_dir(workdir)
    remove_if_exists(log_path)
    nested_config.write_text(
        build_pgsc_calc_resource_config(max_memory=context["manifest"]["runtime"]["max_memory"],
                                       max_cpus=int(context["manifest"]["runtime"]["max_cpus"]))
        + f"\ndocker.runOptions = '-u $(id -u):$(id -g) --label {docker_label}'\n",
        encoding="utf-8",
    )

    args = [
        "nextflow",
        "run",
        "pgscatalog/pgsc_calc",
        "-r", PGSC_CALC_REVISION,
        "-resume",
        "-profile", context["manifest"]["runtime"]["profile"],
        "-c", str(nested_config),
        "--input", str(samplesheet),
        "--pgs_id", ",".join(ids),
        "--target_build", "GRCh38",
        "--min_overlap", str(context["manifest"]["analysis"]["prs"]["min_overlap"]),
        "--max_memory", context["manifest"]["runtime"]["max_memory"],
        "--max_cpus", str(context["manifest"]["runtime"]["max_cpus"]),
        "--outdir", str(outdir),
        "--genotypes_cache", str(genotypes_cache),
        "-work-dir", str(workdir),
    ]
    if context["manifest"]["analysis"]["prs"].get("run_ancestry"):
        reference_db = ensure_reference_db_assets(context)
        args.extend(["--run_ancestry", reference_db["archive_path"]])

    env = {
        "NXF_HOME": context["runtime"]["nextflow_home"],
    }
    try:
        completed = run_command(args, cwd=stage_root, env=env, log_path=log_path, check=False)
    finally:
        # Nextflow starts daemon-owned containers outside its process group. Reap
        # only this scoring cache's containers, including after forced JVM shutdown.
        docker_bin = (context.get("host_checks") or {}).get("docker") or "docker"
        try:
            containers = subprocess.run([docker_bin, "ps", "-aq", "--filter", f"label={docker_label}"],
                                        capture_output=True, text=True, timeout=5, check=True).stdout.split()
            if containers:
                subprocess.run([docker_bin, "rm", "-f", *containers], capture_output=True, timeout=10, check=True)
        except (OSError, subprocess.SubprocessError) as exc:
            print(f"Could not verify scoring container cleanup for {docker_label}: {exc}", file=sys.stderr)
    return completed.returncode, log_path, scoring_input_info


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
    has_scoring_surface = complete_artifacts_exist(
        [
            artifacts.get("normalized_scoring_input_vcf"),
        ]
    ) or complete_artifacts_exist(
        [
            artifacts.get("scoring_input_pgen"),
            artifacts.get("scoring_input_psam"),
            artifacts.get("scoring_input_pvar"),
        ]
    )
    if not has_scoring_surface:
        return None
    required_artifacts = [
        artifacts.get("aggregated_scores"),
        artifacts.get("scoring_input_vcf"),
    ]
    if not complete_artifacts_exist(required_artifacts):
        return None
    dense_artifacts = (payload.get("dense_genotypes") or {}).get("artifacts") or {}
    imputed_artifacts = (payload.get("imputed_genotypes") or {}).get("artifacts") or {}
    if not complete_artifacts_exist([dense_artifacts.get("dense_vcf")]) and not complete_artifacts_exist([imputed_artifacts.get("imputed_vcf")]):
        return None
    payload["reused"] = True
    return payload


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

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
            dense_vcf = dense_root / "dense" / f"{context['sample_id']}.prs_scaffold.vcf.gz"
            dense_qc_path = dense_root / "dense" / "genotype_qc.json"

            if completed_vcf(dense_vcf) and dense_qc_path.is_file():
                dense_info = stage_payload(
                    "prs_dense_genotypes",
                    "ok",
                    reused=True,
                    score_hash=score_hash,
                    summary={**(union_stage.get("summary") or {}), "genotype_qc": load_json(dense_qc_path)},
                    artifacts={
                        "dense_vcf": str(dense_vcf),
                        "genotype_qc_json": str(dense_qc_path),
                        "dense_vcf_index": str(dense_vcf) + ".tbi",
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                        "reference_db_archive": union_stage["artifacts"].get("reference_db_archive"),
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
                    default_worker_cap=6,
                    default_pairhmm_cap=3,
                    min_memory_per_worker_gb=2.0,
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
                        "reference_db_archive": union_stage["artifacts"].get("reference_db_archive"),
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
                            gvcf_root=stage_root.parent / "germline_gvcf",
                            min_depth=cfg.get("min_depth", 10),
                            min_gq=cfg.get("min_gq", 20),
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
                genotype_qc = {key: sum(item["genotype_qc"][key] for item in shard_results)
                               for key in ("total_sites", "callable_sites", "missing_sites")}
                genotype_qc["reasons"] = {}
                for item in shard_results:
                    for reason, count in item["genotype_qc"]["reasons"].items():
                        genotype_qc["reasons"][reason] = genotype_qc["reasons"].get(reason, 0) + count
                write_json(dense_qc_path, genotype_qc)
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
                    summary={**(union_stage.get("summary") or {}), "genotype_qc": genotype_qc},
                    worker_plan=worker_plan,
                    artifacts={
                        "dense_vcf": dense_vcf_path,
                        "genotype_qc_json": str(dense_qc_path),
                        "dense_vcf_index": dense_vcf_path + ".tbi",
                        "targets_vcf": union_stage["artifacts"]["targets_vcf"],
                        "targets_vcf_index": union_stage["artifacts"]["targets_vcf_index"],
                        "reference_db_archive": union_stage["artifacts"].get("reference_db_archive"),
                        "pre_match_audit_json": union_stage["artifacts"].get("pre_match_audit_json"),
                        "pre_match_audit_tsv": union_stage["artifacts"].get("pre_match_audit_tsv"),
                    },
                )
                progress = {
                    **(progress or {}),
                    "phase": "scoring_dense_vcf",
                }
                write_live_status("running", phase="scoring_dense_vcf")

        if dense_info and imputation_cfg.get("enabled"):
            if reference_stage.get("status") != "ready":
                raise ValueError("PRS imputation requires a validated reference, including for imported scaffolds")
            reference_fasta = pathlib.Path(reference_stage["fasta"])
            def update_imputation_status(current_progress: dict[str, Any], current_payload: dict[str, Any]) -> None:
                nonlocal progress, imputed_info
                progress = current_progress
                imputed_info = current_payload
                write_live_status("running", phase=current_progress.get("phase"))

            imputed_info, progress = build_or_reuse_imputed_vcf(
                context=context,
                stage_root=stage_root,
                dense_info=dense_info,
                score_hash=score_hash,
                reference_fasta=reference_fasta,
                update_status=update_imputation_status,
            )
            write_live_status("running", phase=(progress or {}).get("phase"))

        scoring_input = imputed_info or dense_info
        scoring_input_artifacts = (scoring_input or {}).get("artifacts") or {}
        scoring_input_path = scoring_input_artifacts.get("imputed_vcf") or scoring_input_artifacts.get("dense_vcf")
        scoring_input_index = scoring_input_artifacts.get("imputed_vcf_index") or scoring_input_artifacts.get("dense_vcf_index")
        scoring_input_kind = "imputed_vcf" if imputed_info else "dense_vcf"
        dosage_field = imputation_cfg.get("dosage_field") if imputed_info else None
        scoring_identity = scoring_input_identity(context, pathlib.Path(scoring_input_path), scoring_input_kind, dosage_field)
        ancestry_fingerprint = (ensure_reference_db_assets(context)["archive_fingerprint"]
                                if cfg.get("run_ancestry") else None)
        cache_key = hash_strings([json.dumps({
            "scoring_input": scoring_identity, "pgs_ids": sorted(cfg["pgs_ids"]),
            "min_overlap": cfg["min_overlap"], "run_ancestry": bool(cfg.get("run_ancestry")),
            "ancestry_reference": ancestry_fingerprint, "revision": PGSC_CALC_REVISION,
        }, sort_keys=True)])
        scoring_root = ensure_dir(stage_root / score_hash / "scoring" / cache_key)
        stage_manifest_path = scoring_root / "stage_manifest.json"
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
        scoring_vcf_info = pfile_info
        scoring_input_artifacts = {
            "scoring_input_vcf": scoring_input_path,
            "scoring_input_index": scoring_input_index,
            "scoring_input_format": scoring_vcf_info.get("format"),
            "scoring_input_pgen": scoring_vcf_info["artifacts"].get("pgen"),
            "scoring_input_psam": scoring_vcf_info["artifacts"].get("psam"),
            "scoring_input_pvar": scoring_vcf_info["artifacts"].get("pvar"),
            "scoring_input_pfile_log": scoring_vcf_info["artifacts"].get("log"),
            "normalized_scoring_input_vcf": scoring_vcf_info["artifacts"].get("normalized_vcf"),
            "normalized_scoring_input_vcf_index": scoring_vcf_info["artifacts"].get("normalized_vcf_index"),
            "scoring_input_rename_map": scoring_vcf_info["artifacts"].get("rename_map"),
        }
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
                    "sample_sex": cfg.get("sample_sex"),
                    "sex_chromosomes_allowed": bool(cfg.get("sample_sex")),
                    "scoring_input_kind": scoring_input_kind,
                    "pgsc_input_format": scoring_vcf_info.get("format"),
                    "dosage_field": dosage_field,
                    "renamed_chr_prefix": scoring_vcf_info["renamed_chr_prefix"],
                },
                artifacts={
                    "dense_vcf": dense_info["artifacts"]["dense_vcf"] if dense_info else None,
                    "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"] if dense_info else None,
                    "targets_vcf": dense_info["artifacts"].get("targets_vcf") if dense_info else None,
                    "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index") if dense_info else None,
                    "pre_match_audit_json": dense_info["artifacts"].get("pre_match_audit_json") if dense_info else None,
                    "pre_match_audit_tsv": dense_info["artifacts"].get("pre_match_audit_tsv") if dense_info else None,
                    "imputed_vcf": imputed_info["artifacts"]["imputed_vcf"] if imputed_info else None,
                    "imputed_vcf_index": imputed_info["artifacts"]["imputed_vcf_index"] if imputed_info else None,
                    "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                    "pgsc_log": str(log_path),
                    **scoring_input_artifacts,
                    **match_diagnostic_artifacts,
                },
            )
            write_json(stage_manifest_path, payload)
            write_json(live_status_path, payload)
            write_json(args.output, payload)
            return

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
                "pgsc_calc_revision": PGSC_CALC_REVISION,
                "genotype_qc": (scoring_input.get("summary") or {}).get("genotype_qc"),
                "score_metadata": ((dense_info or {}).get("summary") or {}).get("score_metadata", {}),
                "union_sites": (dense_info.get("summary") or {}).get("union_sites") if dense_info else None,
                "overlap_failures": overlap_failures,
                "sample_sex": cfg.get("sample_sex"),
                "sex_chromosomes_allowed": bool(cfg.get("sample_sex")),
                "scoring_input_kind": scoring_input_kind,
                "pgsc_input_format": scoring_vcf_info.get("format"),
                "renamed_chr_prefix": scoring_vcf_info["renamed_chr_prefix"],
            },
            artifacts={
                "dense_vcf": dense_info["artifacts"]["dense_vcf"] if dense_info else None,
                "dense_vcf_index": dense_info["artifacts"]["dense_vcf_index"] if dense_info else None,
                "targets_vcf": dense_info["artifacts"].get("targets_vcf") if dense_info else None,
                "targets_vcf_index": dense_info["artifacts"].get("targets_vcf_index") if dense_info else None,
                "pre_match_audit_json": dense_info["artifacts"].get("pre_match_audit_json") if dense_info else None,
                "pre_match_audit_tsv": dense_info["artifacts"].get("pre_match_audit_tsv") if dense_info else None,
                "imputed_vcf": imputed_info["artifacts"]["imputed_vcf"] if imputed_info else None,
                "imputed_vcf_index": imputed_info["artifacts"]["imputed_vcf_index"] if imputed_info else None,
                "aggregated_scores": aggregated_scores,
                "pgsc_report_html": report_html,
                "match_summary": match_summary,
                "pgsc_outdir": str(scoring_root / "pgsc_calc"),
                "pgsc_log": str(log_path),
                **scoring_input_artifacts,
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
