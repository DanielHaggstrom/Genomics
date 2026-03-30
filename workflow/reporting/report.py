#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import shutil
from typing import Any

from workflow.core.runtime import ensure_dir, load_json, write_json


def copy_file_if_needed(src: str, dst: pathlib.Path) -> None:
    src_path = pathlib.Path(src)
    try:
        if src_path.resolve() == dst.resolve():
            return
    except FileNotFoundError:
        pass
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_path, dst)


def read_prs_rows(path: str | None, limit: int = 20) -> list[dict[str, str]]:
    if not path or not pathlib.Path(path).exists():
        return []
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            rows.append({"PGS": row.get("PGS", ""), "IID": row.get("IID", ""), "SUM": row.get("SUM", "")})
            if idx + 1 >= limit:
                break
    return rows


def render_report(
    context: dict[str, Any],
    monogenic_stage: dict[str, Any],
    prs_stage: dict[str, Any],
    *,
    context_source: pathlib.Path | None,
    monogenic_source: pathlib.Path | None,
    prs_source: pathlib.Path | None,
    output_summary: pathlib.Path,
) -> None:
    run_dir = pathlib.Path(context["results"]["run_dir"])
    report_dir = ensure_dir(pathlib.Path(context["results"]["report_dir"]))
    manifests_dir = ensure_dir(pathlib.Path(context["results"]["manifests_dir"]))

    if context_source is not None:
        copy_file_if_needed(str(context_source), manifests_dir / "prepared_context.json")
    if monogenic_source is not None:
        copy_file_if_needed(str(monogenic_source), manifests_dir / "monogenic_stage.json")
    if prs_source is not None:
        copy_file_if_needed(str(prs_source), manifests_dir / "prs_stage.json")

    summary_artifacts = {
        "context": context,
        "monogenic_stage": monogenic_stage,
        "prs_stage": prs_stage,
    }
    artifact_manifest_path = report_dir / "artifact_manifest.json"
    write_json(artifact_manifest_path, summary_artifacts)

    report_path = report_dir / "report.md"
    reference_probe = context["reference_probe"]
    reference_stage = context["reference_bundle"]
    prs_dense = prs_stage.get("dense_genotypes") or {}
    prs_rows = read_prs_rows((prs_stage.get("artifacts") or {}).get("aggregated_scores"))

    with open(report_path, "w", encoding="utf-8") as out:
        out.write(f"# Genome Report: {context['sample_id']}\n\n")
        out.write("Research-use summary for the canonical CRAM-first pipeline.\n\n")

        out.write("## Summary\n")
        out.write(f"- Sample ID: `{context['sample_id']}`\n")
        out.write(f"- Input fingerprint: `{context['input_fingerprint']['id']}`\n")
        out.write(f"- Reference probe: `{reference_probe.get('status', 'unknown')}`\n")
        out.write(f"- Reference bundle: `{reference_stage.get('status', 'unknown')}`\n")
        out.write(f"- Monogenic stage: `{monogenic_stage.get('status', 'unknown')}`\n")
        out.write(f"- PRS dense genotypes: `{prs_dense.get('status', 'unknown')}`\n")
        out.write(f"- PRS stage: `{prs_stage.get('status', 'unknown')}`\n\n")

        out.write("## Inputs And Runtime\n")
        out.write(f"- CRAM input: `{context['inputs']['original']['cram']}`\n")
        out.write(f"- Active CRAM path: `{context['inputs']['active']['cram']}`\n")
        out.write(f"- CRAM staged into Linux cache: `{context['inputs']['staged']['cram']}`\n")
        out.write(f"- Active reference FASTA: `{reference_stage.get('fasta', 'n/a')}`\n")
        out.write(f"- Linux cache root: `{context['runtime']['linux_cache_root']}`\n")
        out.write(f"- Results directory: `{run_dir}`\n")
        if context.get("upstream_sample_name"):
            out.write(f"- Upstream sample name: `{context['upstream_sample_name']}`\n")
        out.write("\n")

        out.write("## Reference Resolution\n")
        if reference_probe.get("status") == "ok":
            probe_summary = reference_probe.get("summary") or {}
            out.write(f"- Contigs in CRAM header: `{probe_summary.get('contig_count', 'n/a')}`\n")
            out.write(f"- First contigs: `{', '.join(probe_summary.get('first_contigs') or []) or 'n/a'}`\n")
            out.write(f"- Contig prefix: `{probe_summary.get('contig_prefix', '')}`\n")
            out.write(f"- All header contigs have MD5: `{probe_summary.get('all_contigs_have_md5', False)}`\n")
            hints = reference_probe.get("provenance_hints") or {}
            if hints.get("aligner"):
                out.write(f"- Aligner hint: `{hints['aligner']}`\n")
            if hints.get("ref_dir"):
                out.write(f"- Reference directory hint: `{hints['ref_dir']}`\n")
            remote_probe = reference_probe.get("remote_decode_probe") or {}
            out.write(f"- Remote MD5 decode probe: `{remote_probe.get('status', 'not_attempted')}`\n")
        else:
            out.write(f"- Reference probe failed: {reference_probe.get('message', '')}\n")
        if reference_stage.get("comparison"):
            comparison = reference_stage["comparison"]
            out.write(f"- Exact dictionary match: `{comparison.get('exact_dictionary_match', False)}`\n")
            out.write(f"- Extra contigs in candidate bundle: `{len(comparison.get('extra_contigs') or [])}`\n")
        if reference_stage.get("message"):
            out.write(f"- Bundle note: {reference_stage['message']}\n")
        out.write("\n")

        out.write("## Monogenic\n")
        if monogenic_stage.get("status") == "ok":
            out.write(f"- Source mode: `{monogenic_stage.get('source_mode', 'unknown')}`\n")
            import_details = monogenic_stage.get("import_details")
            if import_details:
                out.write(f"- Import provenance: `{import_details.get('provenance')}`\n")
            summary = monogenic_stage.get("summary") or {}
            out.write(f"- Priority rows after dedup: `{summary.get('priority_rows_after_dedup', 'n/a')}`\n")
            out.write(f"- Variant records scanned: `{summary.get('total_variant_records_seen', 'n/a')}`\n")
            artifacts = monogenic_stage.get("artifacts") or {}
            out.write(f"- Annotated VCF: `{artifacts.get('annotated_vcf', 'n/a')}`\n")
            out.write(f"- ClinVar unique TSV: `{artifacts.get('clinvar_unique_tsv', 'n/a')}`\n")
            out.write(f"- LoF unique TSV: `{artifacts.get('lof_unique_tsv', 'n/a')}`\n")
            out.write(f"- Review queue unique TSV: `{artifacts.get('review_unique_tsv', 'n/a')}`\n")
        else:
            out.write(f"- Status: `{monogenic_stage.get('status', 'unknown')}`\n")
            if monogenic_stage.get("message"):
                out.write(f"- Detail: {monogenic_stage['message']}\n")
        out.write("\n")

        out.write("## PRS\n")
        if prs_dense:
            out.write(f"- Dense genotype stage status: `{prs_dense.get('status', 'unknown')}`\n")
            import_details = prs_dense.get("import_details")
            if import_details:
                out.write(f"- Dense genotype import provenance: `{import_details.get('provenance')}`\n")
            dense_summary = prs_dense.get("summary") or {}
            if dense_summary.get("union_sites") is not None:
                out.write(f"- Union sites requested: `{dense_summary['union_sites']}`\n")
            dense_artifacts = prs_dense.get("artifacts") or {}
            if dense_artifacts.get("dense_vcf"):
                out.write(f"- Dense PRS VCF: `{dense_artifacts['dense_vcf']}`\n")
            if dense_artifacts.get("targets_vcf"):
                out.write(f"- Union targets VCF: `{dense_artifacts['targets_vcf']}`\n")

        out.write(f"- PRS scoring status: `{prs_stage.get('status', 'unknown')}`\n")
        if prs_stage.get("status") == "ok":
            artifacts = prs_stage.get("artifacts") or {}
            out.write(f"- Aggregated scores: `{artifacts.get('aggregated_scores', 'n/a')}`\n")
            if artifacts.get("pgsc_report_html"):
                out.write(f"- pgsc_calc HTML report: `{artifacts['pgsc_report_html']}`\n")
            if prs_rows:
                out.write("\n| PGS | IID | SUM |\n")
                out.write("| --- | --- | --- |\n")
                for row in prs_rows:
                    out.write(f"| {row['PGS']} | {row['IID']} | {row['SUM']} |\n")
                out.write("\n")
        else:
            diagnostics = prs_stage.get("diagnostics") or {}
            if prs_stage.get("message"):
                out.write(f"- Detail: {prs_stage['message']}\n")
            overlap_failures = diagnostics.get("overlap_failures") or []
            if overlap_failures:
                out.write("- Matching threshold diagnostics:\n")
                for item in overlap_failures:
                    out.write(f"  - `{item['score']}` matched `{item['match_percent']}%`\n")
            if diagnostics.get("log_path"):
                out.write(f"- Log: `{diagnostics['log_path']}`\n")
        out.write("\n")

        out.write("## Limitations\n")
        out.write("- Research use only; this pipeline does not produce clinical diagnoses.\n")
        out.write("- Monogenic and PRS stages reuse imported artifacts only at explicit stage boundaries.\n")
        out.write("- PRS is generated from a dense genotype VCF derived from the CRAM, not from a sparse hard-filtered WGS VCF.\n")
        out.write("- CRAM-based stages require the exact matching reference bundle.\n\n")

        out.write("## Artifacts\n")
        out.write(f"- Context + stage manifests: `{manifests_dir}`\n")
        out.write(f"- Artifact manifest: `{artifact_manifest_path}`\n")

    summary = {
        "sample_id": context["sample_id"],
        "report_path": str(report_path),
        "artifact_manifest_path": str(artifact_manifest_path),
        "enabled_stages": {
            "monogenic": context["manifest"]["analysis"]["monogenic"]["enabled"],
            "prs": context["manifest"]["analysis"]["prs"]["enabled"],
        },
        "stage_statuses": {
            "reference_probe": reference_probe.get("status"),
            "reference_bundle": reference_stage.get("status"),
            "monogenic": monogenic_stage.get("status"),
            "prs_dense_genotypes": prs_dense.get("status"),
            "prs": prs_stage.get("status"),
        },
    }
    summary["failed_enabled_stages"] = [
        name
        for name, enabled in summary["enabled_stages"].items()
        if enabled and summary["stage_statuses"][name] == "failed"
    ]
    summary["overall_status"] = "failed" if summary["failed_enabled_stages"] else "ok"
    write_json(output_summary, summary)


def enforce_outcome(summary: dict[str, Any]) -> int:
    if summary.get("failed_enabled_stages"):
        raise SystemExit("Enabled stages failed: " + ", ".join(summary["failed_enabled_stages"]))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    render_parser = subparsers.add_parser("render")
    render_parser.add_argument("--context", required=True)
    render_parser.add_argument("--monogenic-stage", required=True)
    render_parser.add_argument("--prs-stage", required=True)
    render_parser.add_argument("--output", required=True)

    enforce_parser = subparsers.add_parser("enforce")
    enforce_parser.add_argument("--summary", required=True)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "render":
        render_report(
            load_json(args.context),
            load_json(args.monogenic_stage),
            load_json(args.prs_stage),
            context_source=pathlib.Path(args.context),
            monogenic_source=pathlib.Path(args.monogenic_stage),
            prs_source=pathlib.Path(args.prs_stage),
            output_summary=pathlib.Path(args.output),
        )
        return
    if args.command == "enforce":
        raise SystemExit(enforce_outcome(load_json(args.summary)))


if __name__ == "__main__":
    main()
