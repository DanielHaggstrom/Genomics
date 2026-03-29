#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import shutil
from typing import Any

from workflow.core.pipeline_lib import (
    artifacts_exist,
    ensure_dir,
    load_json,
    read_text_if_exists,
    run_command,
    stage_payload,
    write_json,
)


def parse_summary(path: str) -> dict[str, str]:
    out: dict[str, str] = {}
    text = read_text_if_exists(path)
    if not text:
        return out
    for line in text.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        out[key] = value
    return out


STAGE_VERSION = 2


def maybe_reuse(stage_manifest_path: pathlib.Path, expected_key: str) -> dict[str, Any] | None:
    if not stage_manifest_path.exists():
        return None
    payload = load_json(stage_manifest_path)
    if payload.get("cache_key") != expected_key:
        return None
    artifacts = dict(payload.get("artifacts") or {})
    artifacts.pop("monogenic_report", None)
    payload["artifacts"] = artifacts
    required = [
        artifacts.get("priority_variants_tsv"),
        artifacts.get("clinvar_unique_tsv"),
        artifacts.get("lof_unique_tsv"),
    ]
    if not artifacts_exist(required):
        return None
    payload["reused"] = True
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--reference-stage", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    reference_stage = load_json(args.reference_stage)
    cfg = context["manifest"]["analysis"]["monogenic"]
    stage_dir = ensure_dir(pathlib.Path(context["annotation_cache_dir"]) / "monogenic_triage")
    stage_manifest_path = stage_dir / "stage_manifest.json"
    cache_key = f"v{STAGE_VERSION}:{context['input_fingerprint']['id']}:{cfg['mode']}:{cfg['max_af']}"

    try:
        if not cfg["enabled"]:
            payload = stage_payload("monogenic", "skipped", message="Monogenic stage disabled.")
            write_json(args.output, payload)
            return

        reused = maybe_reuse(stage_manifest_path, cache_key)
        if reused:
            write_json(args.output, reused)
            return

        reuse_cfg = context["normalized_reuse"]
        annotated_vcf = reuse_cfg.get("existing_annotated_vcf")
        source_mode = None
        callset_dir = ensure_dir(pathlib.Path(context["monogenic_cache_dir"]))

        if annotated_vcf:
            source_mode = "reuse_existing_annotated_vcf"
        else:
            source_vcf = reuse_cfg.get("existing_vcf")
            if not source_vcf:
                if reference_stage.get("status") != "ready":
                    payload = stage_payload(
                        "monogenic",
                        "failed",
                        cache_key=cache_key,
                        message="No annotated VCF or filtered VCF to reuse, and reference bundle is not ready for CRAM-based calling.",
                    )
                    write_json(stage_manifest_path, payload)
                    write_json(args.output, payload)
                    return

                source_vcf = str(callset_dir / f"{context['sample_id']}.hard-filtered.vcf.gz")
                if not pathlib.Path(source_vcf).exists():
                    run_command(
                        [
                            "bash",
                            "workflow/core/call_variants_from_cram.sh",
                            "--cram",
                            context["normalized_inputs"]["cram"],
                            "--reference-fasta",
                            reference_stage["fasta"],
                            "--outdir",
                            str(callset_dir),
                            "--sample-id",
                            context["sample_id"],
                        ],
                        cwd=context["repo_root"],
                        log_path=stage_dir / "call_variants.log",
                    )
                source_mode = "called_from_cram"
            else:
                cached_source_vcf = callset_dir / f"{context['sample_id']}.reused.hard-filtered.vcf.gz"
                cached_source_tbi = pathlib.Path(str(cached_source_vcf) + ".tbi")
                if not cached_source_vcf.exists():
                    shutil.copy2(source_vcf, cached_source_vcf)
                    source_tbi = pathlib.Path(reuse_cfg.get("existing_vcf_index") or f"{source_vcf}.tbi")
                    if source_tbi.exists():
                        shutil.copy2(source_tbi, cached_source_tbi)
                source_vcf = str(cached_source_vcf)
                source_mode = "reuse_existing_filtered_vcf"

            annotated_vcf = str(pathlib.Path(context["annotation_cache_dir"]) / f"{context['sample_id']}.vep.clinvar_enriched.vcf.gz")
            if not pathlib.Path(annotated_vcf).exists():
                run_command(
                    [
                        "bash",
                        "workflow/monogenic/annotate_with_vep.sh",
                        "--input-vcf",
                        source_vcf,
                        "--output-vcf",
                        annotated_vcf,
                        "--forks",
                        str(cfg["vep_forks"]),
                    ],
                    cwd=context["repo_root"],
                    log_path=stage_dir / "annotate_with_vep.log",
                )

        priority_tsv = stage_dir / "priority_variants.tsv"
        priority_summary = stage_dir / "priority_variants.summary.txt"
        run_command(
            [
                "python3",
                "workflow/monogenic/extract_priority_variants.py",
                "--input-vcf",
                annotated_vcf,
                "--output-tsv",
                str(priority_tsv),
                "--summary-txt",
                str(priority_summary),
            ],
            cwd=context["repo_root"],
            log_path=stage_dir / "extract_priority.log",
        )

        run_command(
            [
                "python3",
                "workflow/monogenic/build_monogenic_shortlists.py",
                "--input-tsv",
                str(priority_tsv),
                "--outdir",
                str(stage_dir),
            ],
            cwd=context["repo_root"],
            log_path=stage_dir / "build_shortlists.log",
        )

        review_queue = stage_dir / "rare_clinvar_review_queue.tsv"
        review_unique = stage_dir / "rare_clinvar_review_unique_variants.tsv"
        run_command(
            [
                "python3",
                "workflow/monogenic/triage_clinvar_enriched.py",
                "--input-vcf",
                annotated_vcf,
                "--output-tsv",
                str(review_queue),
                "--mode",
                cfg["mode"],
                "--max-af",
                str(cfg["max_af"]),
            ],
            cwd=context["repo_root"],
            log_path=stage_dir / "triage_clinvar.log",
        )

        run_command(
            [
                "python3",
                "workflow/monogenic/dedup_shortlist_by_variant.py",
                "--input-tsv",
                str(review_queue),
                "--output-tsv",
                str(review_unique),
            ],
            cwd=context["repo_root"],
            log_path=stage_dir / "dedup_review.log",
        )

        summary = parse_summary(str(priority_summary))
        payload = stage_payload(
            "monogenic",
            "ok",
            cache_key=cache_key,
            source_mode=source_mode,
            reused=False,
            artifacts={
                "annotated_vcf": annotated_vcf,
                "priority_variants_tsv": str(priority_tsv),
                "priority_summary": str(priority_summary),
                "clinvar_unique_tsv": str(stage_dir / "clinvar_plp_unique_variants.tsv"),
                "lof_unique_tsv": str(stage_dir / "lof_unique_variants.tsv"),
                "review_unique_tsv": str(review_unique),
            },
            summary=summary,
        )
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)
    except Exception as exc:
        payload = stage_payload("monogenic", "failed", cache_key=cache_key, message=str(exc))
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
