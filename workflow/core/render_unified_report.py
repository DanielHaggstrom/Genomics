#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import shutil

from workflow.core.pipeline_lib import ensure_dir, load_json, write_json


def copy_if_exists(src: str | None, dst_dir: pathlib.Path) -> str | None:
    if not src:
        return None
    src_path = pathlib.Path(src)
    if not src_path.exists():
        return None
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src_path.name
    shutil.copy2(src_path, dst)
    return str(dst)


def copy_file_if_needed(src: str, dst: pathlib.Path) -> None:
    src_path = pathlib.Path(src)
    try:
        if src_path.resolve() == dst.resolve():
            return
    except FileNotFoundError:
        pass
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_path, dst)


def display_root(raw_base: str, *parts: str) -> str:
    base = pathlib.PurePosixPath(str(raw_base).replace("\\", "/"))
    return str(base.joinpath(*parts))


def display_known_path(path: str | None, context: dict, *, run_only: bool = False) -> str | None:
    if not path:
        return None

    target = pathlib.Path(path)
    run_dir = pathlib.Path(context["results_run_dir"])
    sample_cache_dir = pathlib.Path(context["cache_dir"]) / "samples" / context["input_fingerprint"]["id"]
    outputs_cfg = context["manifest"]["outputs"]
    run_display_root = display_root(outputs_cfg["base_dir"], context["sample_id"], outputs_cfg["run_id"])
    cache_display_root = display_root(outputs_cfg["cache_dir"], "samples", context["input_fingerprint"]["id"])

    try:
        rel = target.relative_to(run_dir)
        return str(pathlib.PurePosixPath(run_display_root) / pathlib.PurePosixPath(rel.as_posix()))
    except ValueError:
        pass

    if run_only:
        return path

    try:
        rel = target.relative_to(sample_cache_dir)
        return str(pathlib.PurePosixPath(cache_display_root) / pathlib.PurePosixPath(rel.as_posix()))
    except ValueError:
        return path


def read_prs_rows(path: str | None, limit: int = 20) -> list[dict[str, str]]:
    if not path or not pathlib.Path(path).exists():
        return []
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            rows.append({"PGS": row.get("PGS", ""), "IID": row.get("IID", ""), "SUM": row.get("SUM", "")})
            if idx + 1 >= limit:
                break
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--reference-probe", required=True)
    parser.add_argument("--reference-stage", required=True)
    parser.add_argument("--monogenic-stage", required=True)
    parser.add_argument("--prs-stage", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    reference_probe = load_json(args.reference_probe)
    reference_stage = load_json(args.reference_stage)
    monogenic_stage = load_json(args.monogenic_stage)
    prs_stage = load_json(args.prs_stage)

    run_dir = pathlib.Path(context["results_run_dir"])
    report_dir = ensure_dir(run_dir / "report")
    manifests_dir = ensure_dir(run_dir / "manifests")
    monogenic_dir = ensure_dir(run_dir / "monogenic")
    prs_dir = ensure_dir(run_dir / "prs")

    copy_file_if_needed(args.context, manifests_dir / "sample_context.json")
    copy_file_if_needed(args.reference_probe, manifests_dir / "reference_probe.json")
    copy_file_if_needed(args.reference_stage, manifests_dir / "reference_stage.json")
    copy_file_if_needed(args.monogenic_stage, manifests_dir / "monogenic_stage.json")
    copy_file_if_needed(args.prs_stage, manifests_dir / "prs_stage.json")

    artifact_manifest = {
        "context": context,
        "reference_probe": reference_probe,
        "reference_stage": reference_stage,
        "monogenic_stage": monogenic_stage,
        "prs_stage": prs_stage,
    }
    artifact_manifest_path = report_dir / "artifact_manifest.json"
    write_json(artifact_manifest_path, artifact_manifest)

    copy_if_exists((monogenic_stage.get("artifacts") or {}).get("priority_variants_tsv"), monogenic_dir)
    copy_if_exists((monogenic_stage.get("artifacts") or {}).get("clinvar_unique_tsv"), monogenic_dir)
    copy_if_exists((monogenic_stage.get("artifacts") or {}).get("lof_unique_tsv"), monogenic_dir)
    copy_if_exists((monogenic_stage.get("artifacts") or {}).get("review_unique_tsv"), monogenic_dir)

    copied_scores = copy_if_exists((prs_stage.get("artifacts") or {}).get("aggregated_scores"), prs_dir)
    copied_prs_report = copy_if_exists((prs_stage.get("artifacts") or {}).get("pgsc_report_html"), prs_dir)
    copy_if_exists((prs_stage.get("artifacts") or {}).get("match_summary"), prs_dir)

    prs_rows = read_prs_rows(copied_scores)
    display_run_dir = display_root(
        context["manifest"]["outputs"]["base_dir"],
        context["sample_id"],
        context["manifest"]["outputs"]["run_id"],
    )
    report_path = report_dir / "report.md"
    with open(report_path, "w", encoding="utf-8") as out:
        out.write(f"# Genome Report: {context['sample_id']}\n\n")
        out.write("This report is a research-use summary generated by the CRAM-first rebuild. ")
        out.write("It combines provenance, monogenic triage, and PRS execution status in one place.\n\n")

        out.write("## Summary\n")
        out.write(f"- Sample ID: `{context['sample_id']}`\n")
        out.write(f"- Canonical input fingerprint: `{context['input_fingerprint']['id']}`\n")
        out.write(f"- CRAM reference probe: `{reference_probe.get('status', 'unknown')}`\n")
        out.write(f"- Monogenic stage: `{monogenic_stage.get('status', 'unknown')}`\n")
        out.write(f"- PRS stage: `{prs_stage.get('status', 'unknown')}`\n")
        out.write(f"- Reference bundle stage: `{reference_stage.get('status', 'unknown')}`\n\n")

        out.write("## Inputs And Provenance\n")
        out.write(f"- CRAM: `{context['manifest']['inputs']['cram']}`\n")
        out.write(f"- CRAI: `{context['manifest']['inputs']['crai']}`\n")
        if context.get("upstream_sample_name"):
            out.write(f"- Upstream sample name seen in source files: `{context['upstream_sample_name']}`\n")
        out.write(f"- Contig prefix detected: `{context.get('contig_prefix') or '(none)'}`\n")
        out.write(f"- Results directory: `{display_run_dir}`\n")
        out.write(f"- Artifact manifest: `{display_known_path(str(artifact_manifest_path), context, run_only=True)}`\n\n")

        out.write("## Reference Resolution\n")
        if reference_probe.get("status") == "ok":
            probe_summary = reference_probe.get("summary") or {}
            hints = reference_probe.get("provenance_hints") or {}
            remote_probe = reference_probe.get("remote_decode_probe") or {}
            out.write("- CRAM header probe: completed successfully.\n")
            out.write(f"- Contigs in CRAM header: `{probe_summary.get('contig_count', 'n/a')}`\n")
            out.write(f"- First contigs: `{', '.join(probe_summary.get('first_contigs') or []) or 'n/a'}`\n")
            out.write(f"- All header contigs have MD5: `{probe_summary.get('all_contigs_have_md5', False)}`\n")
            if hints.get("aligner"):
                out.write(f"- Aligner hint: `{hints['aligner']}`\n")
            if hints.get("ht_reference"):
                out.write(f"- DRAGEN `--ht-reference` hint: `{hints['ht_reference']}`\n")
            if hints.get("ref_dir"):
                out.write(f"- DRAGEN `--ref-dir` hint: `{hints['ref_dir']}`\n")
            if remote_probe.get("status") == "ok":
                out.write("- Remote MD5 decode probe: succeeded.\n")
            elif remote_probe.get("status") == "failed":
                out.write(f"- Remote MD5 decode probe: failed. {remote_probe.get('message', '')}\n")
            else:
                out.write(f"- Remote MD5 decode probe: not attempted. {remote_probe.get('message', '')}\n")
        else:
            out.write(f"- CRAM header probe failed. {reference_probe.get('message', '')}\n")

        if reference_stage.get("status") == "ready":
            comparison = reference_stage.get("comparison") or {}
            out.write("- Local reference bundle: validated for all CRAM contigs.\n")
            out.write(f"- Bundle manifest: `{reference_stage.get('bundle_manifest', 'n/a')}`\n")
            out.write(f"- Exact dictionary match: `{comparison.get('exact_dictionary_match', False)}`\n")
            out.write(f"- Extra contigs in candidate bundle: `{len(comparison.get('extra_contigs') or [])}`\n\n")
        elif reference_stage.get("status") == "skipped":
            out.write(f"- Local reference bundle: not ready. {reference_stage.get('message', '')}\n\n")
        else:
            out.write(f"- Local reference bundle validation failed. {reference_stage.get('message', '')}\n")
            comparison = reference_stage.get("comparison") or {}
            if comparison:
                out.write(f"- Missing contigs: `{len(comparison.get('missing_contigs') or [])}`\n")
                out.write(f"- Length mismatches: `{len(comparison.get('length_mismatches') or [])}`\n")
                out.write(f"- MD5 mismatches: `{len(comparison.get('md5_mismatches') or [])}`\n")
                out.write(f"- Candidate contigs without MD5: `{len(comparison.get('md5_missing_in_candidate') or [])}`\n")
            out.write("\n")

        out.write("## Monogenic\n")
        if monogenic_stage.get("status") == "ok":
            out.write("- Status: completed successfully.\n")
            summary = monogenic_stage.get("summary") or {}
            if summary:
                out.write(f"- Priority transcript rows after dedup: `{summary.get('priority_rows_after_dedup', 'n/a')}`\n")
                out.write(f"- VCF records scanned: `{summary.get('total_variant_records_seen', 'n/a')}`\n")
            out.write("\n")
        elif monogenic_stage.get("status") == "skipped":
            out.write(f"- Status: skipped. {monogenic_stage.get('message', '')}\n\n")
        else:
            out.write(f"- Status: failed. {monogenic_stage.get('message', '')}\n\n")

        out.write("## PRS\n")
        if prs_stage.get("status") == "ok":
            out.write("- Status: completed successfully.\n")
            if copied_scores:
                out.write(f"- Aggregated scores: `{copied_scores}`\n")
            if copied_prs_report:
                out.write(f"- pgsc_calc HTML report: `{copied_prs_report}`\n")
            if prs_rows:
                out.write("\n| PGS | IID | SUM |\n")
                out.write("| --- | --- | --- |\n")
                for row in prs_rows:
                    out.write(f"| {row['PGS']} | {row['IID']} | {row['SUM']} |\n")
                out.write("\n")
        elif prs_stage.get("status") == "failed":
            out.write(f"- Status: failed. {prs_stage.get('message', '')}\n")
            diagnostics = prs_stage.get("diagnostics") or {}
            overlap_failures = diagnostics.get("overlap_failures") or []
            if overlap_failures:
                out.write("- Matching threshold diagnostics:\n")
                for item in overlap_failures:
                    out.write(f"  - `{item['score']}` matched `{item['match_percent']}%`\n")
            if diagnostics.get("log_path"):
                out.write(f"- PRS log: `{display_known_path(diagnostics['log_path'], context)}`\n")
            out.write("\n")
        else:
            out.write(f"- Status: skipped. {prs_stage.get('message', '')}\n\n")

        out.write("## Limitations\n")
        out.write("- This pipeline is for research use and does not produce clinical diagnoses.\n")
        out.write("- PRS on WGS remains experimental here: the rebuild force-genotypes score loci from the CRAM instead of trusting a sparse variant-only VCF.\n")
        out.write("- CRAM-based stages require the exact reference bundle that matches the CRAM. A generic GRCh38 FASTA is not always sufficient.\n\n")

        out.write("## Artifacts\n")
        out.write(f"- Stage manifests live in `{display_known_path(str(manifests_dir), context, run_only=True)}`\n")
        out.write(f"- Monogenic outputs copied to `{display_known_path(str(monogenic_dir), context, run_only=True)}`\n")
        out.write(f"- PRS outputs copied to `{display_known_path(str(prs_dir), context, run_only=True)}`\n")
        out.write(f"- Full artifact manifest: `{display_known_path(str(artifact_manifest_path), context, run_only=True)}`\n")

    enabled_stages = {
        "monogenic": context["manifest"]["analysis"]["monogenic"]["enabled"],
        "prs": context["manifest"]["analysis"]["prs"]["enabled"],
    }
    failed_enabled = [
        stage_name
        for stage_name, enabled in enabled_stages.items()
        if enabled and {"monogenic": monogenic_stage, "prs": prs_stage}[stage_name].get("status") == "failed"
    ]
    summary = {
        "sample_id": context["sample_id"],
        "report_path": str(report_path),
        "artifact_manifest_path": str(artifact_manifest_path),
        "enabled_stages": enabled_stages,
        "stage_statuses": {
            "reference_probe": reference_probe.get("status"),
            "reference_bundle": reference_stage.get("status"),
            "monogenic": monogenic_stage.get("status"),
            "prs": prs_stage.get("status"),
        },
        "failed_enabled_stages": failed_enabled,
        "overall_status": "failed" if failed_enabled else "ok",
    }
    write_json(args.output, summary)


if __name__ == "__main__":
    main()
