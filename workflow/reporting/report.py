#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import shutil
from typing import Any

from workflow.core.manifest import derive_execution_mode
from workflow.core.runtime import ensure_dir, load_json, write_json


def copy_file_if_needed(src: str, dst: pathlib.Path) -> None:
    src_path = pathlib.Path(src)
    if not src_path.exists():
        return
    try:
        if src_path.resolve() == dst.resolve():
            return
    except FileNotFoundError:
        pass
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_path, dst)


def publish_stage_artifacts(
    *,
    run_dir: pathlib.Path,
    monogenic_stage: dict[str, Any],
    prs_stage: dict[str, Any],
) -> dict[str, dict[str, str]]:
    published: dict[str, dict[str, str]] = {"monogenic": {}, "prs": {}}

    monogenic_dir = ensure_dir(run_dir / "monogenic")
    monogenic_artifacts = monogenic_stage.get("artifacts") or {}
    for source_key, filename in [
        ("annotated_vcf", "annotated.vcf.gz"),
        ("priority_summary", "priority_variants.summary.txt"),
        ("priority_variants_tsv", "priority_variants.tsv"),
        ("clinvar_unique_tsv", "clinvar_plp_unique_variants.tsv"),
        ("lof_unique_tsv", "lof_unique_variants.tsv"),
        ("review_unique_tsv", "rare_clinvar_review_unique_variants.tsv"),
    ]:
        source = monogenic_artifacts.get(source_key)
        if not source:
            continue
        destination = monogenic_dir / filename
        copy_file_if_needed(source, destination)
        published["monogenic"][source_key] = str(destination)
        index_path = source + ".tbi" if source.endswith(".vcf.gz") else None
        if index_path and pathlib.Path(index_path).exists():
            copy_file_if_needed(index_path, pathlib.Path(str(destination) + ".tbi"))

    prs_dir = ensure_dir(run_dir / "prs")
    prs_artifacts = prs_stage.get("artifacts") or {}
    prs_imputed_artifacts = (prs_stage.get("imputed_genotypes") or {}).get("artifacts") or {}
    dense_artifacts = (prs_stage.get("dense_genotypes") or {}).get("artifacts") or {}
    for source_key, filename in [
        ("dense_vcf", "prs_scaffold.vcf.gz"),
        ("aggregated_scores", "aggregated_scores.txt.gz"),
        ("pgsc_report_html", "pgsc_report.html"),
        ("match_summary", "match_summary.csv"),
        ("match_diagnostics_json", "match_diagnostics.json"),
        ("overlap_failures_tsv", "overlap_failures.tsv"),
        ("imputed_vcf", "imputed_prs.vcf.gz"),
    ]:
        source = dense_artifacts.get(source_key) or prs_artifacts.get(source_key) or prs_imputed_artifacts.get(source_key)
        if not source:
            continue
        destination = prs_dir / filename
        copy_file_if_needed(source, destination)
        published["prs"][source_key] = str(destination)
        index_path = source + ".tbi" if source.endswith(".vcf.gz") else None
        if index_path and pathlib.Path(index_path).exists():
            copy_file_if_needed(index_path, pathlib.Path(str(destination) + ".tbi"))

    for source_key, filename in [
        ("pre_match_audit_tsv", "pre_match_audit.tsv"),
        ("pre_match_audit_json", "pre_match_audit.json"),
    ]:
        source = dense_artifacts.get(source_key) or prs_artifacts.get(source_key)
        if not source:
            continue
        destination = prs_dir / filename
        copy_file_if_needed(source, destination)
        published["prs"][source_key] = str(destination)

    return published


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


def render_import_details(out: Any, import_details: dict[str, Any] | None, *, prefix: str = "Import") -> None:
    if not import_details:
        return
    out.write(f"- {prefix} provenance: `{import_details.get('provenance', 'unknown')}`\n")
    if import_details.get("source"):
        out.write(f"- {prefix} source path: `{import_details['source']}`\n")
    if import_details.get("active") and import_details.get("active") != import_details.get("source"):
        out.write(f"- {prefix} active path: `{import_details['active']}`\n")
    out.write(f"- {prefix} trusted provenance: `{bool(import_details.get('trusted_provenance'))}`\n")
    if import_details.get("stage_manifest"):
        out.write(f"- {prefix} stage manifest: `{import_details['stage_manifest']}`\n")
    elif not import_details.get("trusted_provenance"):
        out.write(f"- {prefix} note: external import without repo-generated stage manifest\n")


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
    published_outputs = publish_stage_artifacts(run_dir=run_dir, monogenic_stage=monogenic_stage, prs_stage=prs_stage)
    write_json(artifact_manifest_path, summary_artifacts)

    report_path = report_dir / "report.md"
    execution_mode = context.get("execution_mode") or derive_execution_mode(context["manifest"])
    reference_probe = context["reference_probe"]
    reference_stage = context["reference_bundle"]
    prs_dense = prs_stage.get("dense_genotypes") or {}
    prs_imputed = prs_stage.get("imputed_genotypes") or {}
    prs_rows = read_prs_rows((prs_stage.get("artifacts") or {}).get("aggregated_scores"))
    monogenic_enabled = bool(context["manifest"]["analysis"]["monogenic"]["enabled"])
    prs_enabled = bool(context["manifest"]["analysis"]["prs"]["enabled"])
    if execution_mode == "monogenic-only":
        report_title = "Monogenic Report"
    elif execution_mode == "prs-only":
        report_title = "PRS Report"
    else:
        report_title = "Genome Report"

    with open(report_path, "w", encoding="utf-8") as out:
        out.write(f"# {report_title}: {context['sample_id']}\n\n")
        out.write("Research-use summary for the canonical CRAM-first pipeline.\n\n")

        out.write("## Summary\n")
        out.write(f"- Sample ID: `{context['sample_id']}`\n")
        out.write(f"- Execution mode: `{execution_mode}`\n")
        out.write(f"- Input fingerprint: `{context['input_fingerprint']['id']}`\n")
        out.write(f"- Reference probe: `{reference_probe.get('status', 'unknown')}`\n")
        out.write(f"- Reference bundle: `{reference_stage.get('status', 'unknown')}`\n")
        out.write(f"- Monogenic stage: `{monogenic_stage.get('status', 'unknown')}`\n")
        out.write(f"- PRS scaffold stage: `{prs_dense.get('status', 'unknown')}`\n")
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

        if monogenic_enabled:
            out.write("## Monogenic\n")
            if monogenic_stage.get("status") == "ok":
                out.write(f"- Source mode: `{monogenic_stage.get('source_mode', 'unknown')}`\n")
                import_details = monogenic_stage.get("import_details")
                render_import_details(out, import_details)
                summary = monogenic_stage.get("summary") or {}
                out.write(f"- Priority rows after dedup: `{summary.get('priority_rows_after_dedup', 'n/a')}`\n")
                out.write(f"- Variant records scanned: `{summary.get('total_variant_records_seen', 'n/a')}`\n")
                artifacts = monogenic_stage.get("artifacts") or {}
                out.write(f"- Annotated VCF: `{artifacts.get('annotated_vcf', 'n/a')}`\n")
                out.write(f"- ClinVar unique TSV: `{artifacts.get('clinvar_unique_tsv', 'n/a')}`\n")
                out.write(f"- LoF unique TSV: `{artifacts.get('lof_unique_tsv', 'n/a')}`\n")
                out.write(f"- Review queue unique TSV: `{artifacts.get('review_unique_tsv', 'n/a')}`\n")
                published_monogenic = published_outputs.get("monogenic") or {}
                if published_monogenic:
                    out.write(f"- Published monogenic outputs: `{run_dir / 'monogenic'}`\n")
            else:
                out.write(f"- Status: `{monogenic_stage.get('status', 'unknown')}`\n")
                if monogenic_stage.get("message"):
                    out.write(f"- Detail: {monogenic_stage['message']}\n")
            out.write("\n")

        if prs_enabled:
            out.write("## PRS\n")
            if prs_dense:
                out.write(f"- Scaffold stage status: `{prs_dense.get('status', 'unknown')}`\n")
                import_details = prs_dense.get("import_details")
                render_import_details(out, import_details, prefix="Scaffold import")
                dense_summary = prs_dense.get("summary") or {}
                if dense_summary.get("reference_db"):
                    out.write(f"- Reference DB backbone: `{dense_summary['reference_db']}`\n")
                if dense_summary.get("scaffold_sites") is not None:
                    out.write(f"- Scaffold sites requested: `{dense_summary['scaffold_sites']}`\n")
                if dense_summary.get("backbone_sites") is not None:
                    out.write(f"- Backbone sites requested: `{dense_summary['backbone_sites']}`\n")
                if dense_summary.get("score_extra_sites") is not None:
                    out.write(f"- Score-extra sites requested: `{dense_summary['score_extra_sites']}`\n")
                if dense_summary.get("accepted_strand_flip_rows") is not None:
                    out.write(f"- Accepted rows rescued by strand complement: `{dense_summary['accepted_strand_flip_rows']}`\n")
                accepted_rows_by_source = dense_summary.get("accepted_rows_by_source") or {}
                if accepted_rows_by_source:
                    out.write(
                        "- Accepted score rows by source: "
                        + ", ".join(f"`{key}`={value}" for key, value in accepted_rows_by_source.items())
                        + "\n"
                    )
                skipped_reasons = dense_summary.get("skipped_reasons") or {}
                if skipped_reasons:
                    out.write(
                        "- Skipped score rows by reason: "
                        + ", ".join(f"`{key}`={value}" for key, value in skipped_reasons.items())
                        + "\n"
                    )
                dense_artifacts = prs_dense.get("artifacts") or {}
                if dense_artifacts.get("dense_vcf"):
                    out.write(f"- PRS scaffold VCF: `{dense_artifacts['dense_vcf']}`\n")
                if dense_artifacts.get("targets_vcf"):
                    out.write(f"- Scaffold targets VCF: `{dense_artifacts['targets_vcf']}`\n")
                if dense_artifacts.get("reference_db_archive"):
                    out.write(f"- Reference DB archive: `{dense_artifacts['reference_db_archive']}`\n")
                if dense_artifacts.get("pre_match_audit_tsv"):
                    out.write(f"- Pre-match audit TSV: `{dense_artifacts['pre_match_audit_tsv']}`\n")
                if dense_artifacts.get("pre_match_audit_json"):
                    out.write(f"- Pre-match audit JSON: `{dense_artifacts['pre_match_audit_json']}`\n")
                published_prs = published_outputs.get("prs") or {}
                if published_prs:
                    out.write(f"- Published PRS outputs: `{run_dir / 'prs'}`\n")

            if prs_imputed:
                out.write(f"- Imputed genotype stage status: `{prs_imputed.get('status', 'unknown')}`\n")
                render_import_details(out, prs_imputed.get("import_details"), prefix="Imputed genotype import")
                imputed_summary = prs_imputed.get("summary") or {}
                if imputed_summary.get("backend"):
                    out.write(f"- Imputation backend: `{imputed_summary['backend']}`\n")
                if imputed_summary.get("panel_id"):
                    out.write(f"- Imputation panel: `{imputed_summary['panel_id']}`\n")
                if imputed_summary.get("info_threshold") is not None:
                    out.write(f"- Imputation INFO threshold: `{imputed_summary['info_threshold']}`\n")
                imputed_artifacts = prs_imputed.get("artifacts") or {}
                if imputed_artifacts.get("imputed_vcf"):
                    out.write(f"- Imputed PRS VCF: `{imputed_artifacts['imputed_vcf']}`\n")

            out.write(f"- PRS scoring status: `{prs_stage.get('status', 'unknown')}`\n")
            if prs_stage.get("status") == "ok":
                artifacts = prs_stage.get("artifacts") or {}
                summary = prs_stage.get("summary") or {}
                out.write(f"- Scoring input kind: `{summary.get('scoring_input_kind', 'unknown')}`\n")
                if summary.get("pgsc_input_format"):
                    out.write(f"- pgsc_calc input format: `{summary.get('pgsc_input_format')}`\n")
                if summary.get("renamed_chr_prefix") is not None:
                    out.write(f"- Scoring contigs normalized for pgsc_calc: `{summary.get('renamed_chr_prefix')}`\n")
                out.write(f"- Aggregated scores: `{artifacts.get('aggregated_scores', 'n/a')}`\n")
                if artifacts.get("pgsc_report_html"):
                    out.write(f"- pgsc_calc HTML report: `{artifacts['pgsc_report_html']}`\n")
                if artifacts.get("match_diagnostics_json"):
                    out.write(f"- Match diagnostics JSON: `{artifacts['match_diagnostics_json']}`\n")
                if artifacts.get("scoring_input_pgen"):
                    out.write(f"- Scoring input PGEN: `{artifacts['scoring_input_pgen']}`\n")
                if artifacts.get("scoring_input_psam"):
                    out.write(f"- Scoring input PSAM: `{artifacts['scoring_input_psam']}`\n")
                if artifacts.get("scoring_input_pvar"):
                    out.write(f"- Scoring input PVAR: `{artifacts['scoring_input_pvar']}`\n")
                if artifacts.get("normalized_scoring_input_vcf"):
                    out.write(f"- Normalized scoring VCF: `{artifacts['normalized_scoring_input_vcf']}`\n")
                if artifacts.get("overlap_failures_tsv"):
                    out.write(f"- Match threshold summary TSV: `{artifacts['overlap_failures_tsv']}`\n")
                if prs_rows:
                    out.write("\n| PGS | IID | SUM |\n")
                    out.write("| --- | --- | --- |\n")
                    for row in prs_rows:
                        out.write(f"| {row['PGS']} | {row['IID']} | {row['SUM']} |\n")
                    out.write("\n")
            else:
                diagnostics = prs_stage.get("diagnostics") or {}
                artifacts = prs_stage.get("artifacts") or {}
                if prs_stage.get("message"):
                    out.write(f"- Detail: {prs_stage['message']}\n")
                if diagnostics.get("scoring_input_kind"):
                    out.write(f"- Scoring input kind: `{diagnostics['scoring_input_kind']}`\n")
                if diagnostics.get("pgsc_input_format"):
                    out.write(f"- pgsc_calc input format: `{diagnostics['pgsc_input_format']}`\n")
                if diagnostics.get("renamed_chr_prefix") is not None:
                    out.write(f"- Scoring contigs normalized for pgsc_calc: `{diagnostics['renamed_chr_prefix']}`\n")
                if artifacts.get("match_diagnostics_json"):
                    out.write(f"- Match diagnostics JSON: `{artifacts['match_diagnostics_json']}`\n")
                if artifacts.get("scoring_input_pgen"):
                    out.write(f"- Scoring input PGEN: `{artifacts['scoring_input_pgen']}`\n")
                if artifacts.get("scoring_input_psam"):
                    out.write(f"- Scoring input PSAM: `{artifacts['scoring_input_psam']}`\n")
                if artifacts.get("scoring_input_pvar"):
                    out.write(f"- Scoring input PVAR: `{artifacts['scoring_input_pvar']}`\n")
                if artifacts.get("normalized_scoring_input_vcf"):
                    out.write(f"- Normalized scoring VCF: `{artifacts['normalized_scoring_input_vcf']}`\n")
                if artifacts.get("overlap_failures_tsv"):
                    out.write(f"- Match threshold summary TSV: `{artifacts['overlap_failures_tsv']}`\n")
                if artifacts.get("match_summary"):
                    out.write(f"- Match summary CSV: `{artifacts['match_summary']}`\n")
                overlap_failures = diagnostics.get("overlap_failures") or []
                if overlap_failures:
                    out.write("\n| Score | Match % |\n")
                    out.write("| --- | ---: |\n")
                    for item in overlap_failures:
                        out.write(f"| {item['score']} | {item['match_percent']} |\n")
                    out.write("\n")
                if diagnostics.get("log_path"):
                    out.write(f"- Log: `{diagnostics['log_path']}`\n")
            out.write("\n")

        out.write("## Limitations\n")
        out.write("- Research use only; this pipeline does not produce clinical diagnoses.\n")
        out.write("- Monogenic and PRS stages reuse imported artifacts only at explicit stage boundaries.\n")
        out.write("- PRS is generated from a broad scaffold VCF derived from the CRAM, with optional scaffold-based imputation before scoring.\n")
        out.write("- CRAM-based stages require the exact matching reference bundle.\n\n")

        out.write("## Artifacts\n")
        out.write(f"- Context + stage manifests: `{manifests_dir}`\n")
        out.write(f"- Artifact manifest: `{artifact_manifest_path}`\n")

    summary = {
        "sample_id": context["sample_id"],
        "execution_mode": execution_mode,
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
            "prs_imputed_genotypes": prs_imputed.get("status"),
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
