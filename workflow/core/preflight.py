#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import re
import shutil

from workflow.core.pipeline_lib import (
    REPO_ROOT,
    detect_executable,
    ensure_dir,
    hash_strings,
    load_json,
    normalize_path,
    run_command,
    sample_file_fingerprint,
    write_json,
)


def detect_contig_prefix(existing_vcf: str | None, cram: str) -> str:
    if existing_vcf and pathlib.Path(existing_vcf).exists():
        completed = run_command(["bcftools", "view", "-h", existing_vcf], check=True)
        for line in completed.stdout.splitlines():
            if line.startswith("##contig=<ID="):
                contig = line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">")
                return "chr" if contig.startswith("chr") else ""

    completed = run_command(["samtools", "view", "-H", cram], check=True)
    for line in completed.stdout.splitlines():
        if not line.startswith("@SQ") or "SN:" not in line:
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                contig = field[3:]
                return "chr" if contig.startswith("chr") else ""
    return ""


def detect_upstream_sample_name(existing_vcf: str | None, cram: str) -> str | None:
    if existing_vcf and pathlib.Path(existing_vcf).exists():
        completed = run_command(["bcftools", "query", "-l", existing_vcf], check=True)
        names = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
        if names:
            return names[0]

    completed = run_command(["samtools", "view", "-H", cram], check=True)
    for line in completed.stdout.splitlines():
        if not line.startswith("@RG"):
            continue
        match = re.search(r"\bSM:([^\t]+)", line)
        if match:
            return match.group(1)
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validated-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    payload = load_json(args.validated_manifest)
    manifest = payload["manifest"]
    manifest_path = pathlib.Path(payload["manifest_path"])

    normalized_inputs = {
        "cram": normalize_path(manifest["inputs"]["cram"], base_dir=manifest_path.parent),
        "crai": normalize_path(manifest["inputs"]["crai"], base_dir=manifest_path.parent),
    }
    normalized_reuse = {
        key: normalize_path(value, base_dir=manifest_path.parent)
        for key, value in (manifest.get("reuse") or {}).items()
        if value
    }

    cram_path = pathlib.Path(normalized_inputs["cram"])
    crai_path = pathlib.Path(normalized_inputs["crai"])
    if not cram_path.exists():
        raise SystemExit(f"Missing CRAM: {cram_path}")
    if not crai_path.exists():
        raise SystemExit(f"Missing CRAI: {crai_path}")

    for optional_file in normalized_reuse.values():
        if optional_file and not pathlib.Path(optional_file).exists():
            raise SystemExit(f"Configured reuse artifact not found: {optional_file}")

    base_dir = pathlib.Path(normalize_path(manifest["outputs"]["base_dir"], base_dir=REPO_ROOT))
    cache_dir = pathlib.Path(normalize_path(manifest["outputs"]["cache_dir"], base_dir=REPO_ROOT))
    run_id = manifest["outputs"]["run_id"]

    ensure_dir(base_dir)
    ensure_dir(cache_dir)

    cram_fp = sample_file_fingerprint(cram_path)
    crai_fp = sample_file_fingerprint(crai_path)
    input_fingerprint = {
        "id": hash_strings([cram_fp["id"], crai_fp["id"]]),
        "cram": cram_fp,
        "crai": crai_fp,
    }

    sample_cache_dir = ensure_dir(cache_dir / "samples" / input_fingerprint["id"])
    monogenic_cache_dir = ensure_dir(sample_cache_dir / "monogenic_callset")
    annotation_cache_dir = ensure_dir(sample_cache_dir / "annotation")
    prs_cache_root = ensure_dir(sample_cache_dir / "prs_genotypes")
    results_run_dir = ensure_dir(base_dir / manifest["sample_id"] / run_id)
    report_dir = ensure_dir(results_run_dir / "report")
    manifests_dir = ensure_dir(results_run_dir / "manifests")

    usage_target = base_dir if base_dir.exists() else base_dir.parent
    disk_usage = shutil.disk_usage(usage_target)
    free_gb = disk_usage.free / (1024 ** 3)
    min_free_gb = manifest["runtime"]["min_free_gb"]
    if free_gb < min_free_gb:
        raise SystemExit(
            f"Only {free_gb:.1f} GiB free on {usage_target}, but manifest requires at least {min_free_gb} GiB"
        )

    contig_prefix = detect_contig_prefix(normalized_reuse.get("existing_vcf"), normalized_inputs["cram"])
    upstream_sample = detect_upstream_sample_name(normalized_reuse.get("existing_vcf"), normalized_inputs["cram"])

    existing_annotated = normalized_reuse.get("existing_annotated_vcf")
    existing_vcf = normalized_reuse.get("existing_vcf")
    monogenic_reference_required = bool(
        manifest["analysis"]["monogenic"]["enabled"] and not (existing_annotated or existing_vcf)
    )
    prs_reference_required = bool(manifest["analysis"]["prs"]["enabled"])

    host_checks = {
        "docker": detect_executable("docker"),
        "nextflow": detect_executable("nextflow"),
        "java": detect_executable("java"),
        "python3": detect_executable("python3"),
        "bcftools": detect_executable("bcftools"),
        "samtools": detect_executable("samtools"),
        "tabix": detect_executable("tabix"),
    }
    missing = [name for name, path in host_checks.items() if path is None]
    if missing:
        raise SystemExit("Missing required executables: " + ", ".join(sorted(missing)))

    context = {
        "context_version": 1,
        "repo_root": str(REPO_ROOT),
        "manifest_path": str(manifest_path),
        "manifest": manifest,
        "sample_id": manifest["sample_id"],
        "normalized_inputs": normalized_inputs,
        "normalized_reuse": normalized_reuse,
        "input_fingerprint": input_fingerprint,
        "cache_dir": str(cache_dir),
        "sample_cache_dir": str(sample_cache_dir),
        "monogenic_cache_dir": str(monogenic_cache_dir),
        "annotation_cache_dir": str(annotation_cache_dir),
        "prs_cache_root": str(prs_cache_root),
        "results_run_dir": str(results_run_dir),
        "report_dir": str(report_dir),
        "manifests_dir": str(manifests_dir),
        "contig_prefix": contig_prefix,
        "upstream_sample_name": upstream_sample,
        "host_checks": host_checks,
        "disk": {
            "usage_target": str(usage_target),
            "free_gb": round(free_gb, 2),
            "min_free_gb": min_free_gb,
        },
        "stage_requirements": {
            "monogenic_reference_required": monogenic_reference_required,
            "prs_reference_required": prs_reference_required,
            "any_reference_required": monogenic_reference_required or prs_reference_required,
        },
    }
    write_json(args.output, context)


if __name__ == "__main__":
    main()
