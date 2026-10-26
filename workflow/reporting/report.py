#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import tempfile
from typing import Any

from workflow.core.manifest import derive_execution_mode
from workflow.core.runtime import atomic_copy, ensure_dir, load_json, write_json


def copy_file_if_needed(src: str, dst: pathlib.Path) -> bool:
    source = pathlib.Path(src)
    if not source.is_file():
        return False
    if source.resolve() == dst.resolve():
        return True
    if dst.is_file() and (source.stat().st_size, source.stat().st_mtime_ns) == (dst.stat().st_size, dst.stat().st_mtime_ns):
        return True
    atomic_copy(source, dst)
    return True


def publish_stage_artifacts(
    *, run_dir: pathlib.Path, monogenic_stage: dict[str, Any], prs_stage: dict[str, Any],
) -> dict[str, dict[str, str]]:
    # Genotype VCFs are resumable cache artifacts, not reader-facing duplicates.
    sources = {
        "monogenic": monogenic_stage.get("artifacts") or {},
        "prs": {**((prs_stage.get("dense_genotypes") or {}).get("artifacts") or {}),
                **((prs_stage.get("imputed_genotypes") or {}).get("artifacts") or {}),
                **(prs_stage.get("artifacts") or {})},
    }
    names = {
        "monogenic": {
            "priority_summary": "priority_variants.summary.txt",
            "priority_variants_tsv": "priority_variants.tsv",
            "clinvar_unique_tsv": "clinvar_plp_unique_variants.tsv",
            "lof_unique_tsv": "lof_unique_variants.tsv",
            "review_unique_tsv": "rare_clinvar_review_unique_variants.tsv",
        },
        "prs": {
            "genotype_qc_json": "genotype_qc.json",
            "aggregated_scores": "aggregated_scores.txt.gz",
            "pgsc_report_html": "pgsc_report.html",
            "match_summary": "match_summary.csv",
            "match_diagnostics_json": "match_diagnostics.json",
            "overlap_failures_tsv": "overlap_failures.tsv",
            "pre_match_audit_tsv": "pre_match_audit.tsv",
            "pre_match_audit_json": "pre_match_audit.json",
        },
    }
    published = {"monogenic": {}, "prs": {}}
    for stage, files in names.items():
        for key, filename in files.items():
            source = sources[stage].get(key)
            destination = run_dir / stage / filename
            if source and copy_file_if_needed(source, destination):
                published[stage][key] = str(destination)
    return published


def read_prs_rows(path: str | None) -> list[dict[str, str]]:
    if not path or not pathlib.Path(path).exists():
        return []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_match_counts(path: str | None) -> dict[tuple[str, str], dict[str, int]]:
    """Aggregate pgsc_calc summary categories, never its per-candidate match log."""
    counts: dict[tuple[str, str], dict[str, int]] = {}
    if not path or not pathlib.Path(path).exists():
        return counts
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"dataset", "accession", "match_status", "count"}.issubset(reader.fieldnames or []):
            return counts
        for row in reader:
            status = row["match_status"]
            if status not in {"matched", "excluded", "unmatched"}:
                continue
            count = int(row["count"])
            key = (row["dataset"], row["accession"].split("_hmPOS_", 1)[0])
            bucket = counts.setdefault(key, {"matched": 0, "total": 0})
            bucket["total"] += count
            if status == "matched":
                bucket["matched"] += count
    return counts


def markdown_cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def render_prs_results(out: Any, prs_stage: dict[str, Any], prs_rows: list[dict[str, str]], *, run_ancestry: bool) -> None:
    out.write("## PRS Results\n\n")
    summary = prs_stage.get("summary") or {}
    artifacts = prs_stage.get("artifacts") or {}
    qc = summary.get("genotype_qc") or ((prs_stage.get("dense_genotypes") or {}).get("summary") or {}).get("genotype_qc")
    if qc:
        out.write(f"Scaffold genotype QC: {qc.get('callable_sites', 'not reported')} callable, "
                  f"{qc.get('missing_sites', 'not reported')} missing out of {qc.get('total_sites', 'not reported')} sites. "
                  "These are overall scaffold counts, not per-score coverage. Missing observations are excluded before scoring.\n\n")
    reference_rows = sum(row.get("sampleset") == "reference" for row in prs_rows)
    prs_rows = [row for row in prs_rows if row.get("sampleset") != "reference"]
    if prs_stage.get("status") == "ok" and prs_rows:
        out.write(f"Calculated {len(prs_rows)} score rows")
        if summary.get("num_requested_scores") is not None:
            out.write(f" from {summary['num_requested_scores']} requested scores")
        out.write(". All target score rows are shown below.\n\n")
        if reference_rows:
            out.write(f"The source table also contains {reference_rows} reference-panel rows; these are not this sample's results.\n\n")
        matching = read_match_counts(artifacts.get("match_summary"))
        columns = [("PGS", "PGS"), ("trait", "Trait"), ("IID", "Sample"), ("SUM", "Raw SUM")]
        columns += [(key, key) for key in ("DENOM", "AVG") if any(row.get(key) for row in prs_rows)]
        columns += [("matched", "Matched variants"), ("match_percent", "Match %")]
        columns += [(key, label) for key, label in (
            ("usable_genotypes", "Usable genotypes"), ("missing_genotypes", "Missing genotypes")
        ) if any(row.get(key) not in (None, "", "NA", "nan") for row in prs_rows)]
        out.write("| " + " | ".join(label for _, label in columns) + " |\n")
        out.write("| " + " | ".join("---" for _ in columns) + " |\n")
        for row in prs_rows:
            values = dict(row)
            score_id = row.get("PGS", "").split("_hmPOS_", 1)[0]
            values["trait"] = row.get("trait") or row.get("trait_reported") or (summary.get("score_metadata", {}).get(score_id) or {}).get("trait") or "not supplied"
            match = matching.get((row.get("sampleset", ""), score_id))
            if match and match["total"]:
                values["matched"] = f"{match['matched']} / {match['total']}"
                values["match_percent"] = f"{100 * match['matched'] / match['total']:.2f}"
            out.write("| " + " | ".join(markdown_cell(values.get(key) or "not reported") for key, _ in columns) + " |\n")
        out.write("\nSUM is the raw weighted sum of effect-allele dosages. Raw sums are not comparable across different scores "
                  "and are not absolute disease probabilities.\n\n")
        if any(row.get("DENOM") for row in prs_rows):
            out.write("DENOM and AVG are retained as reported by pgsc_calc; AVG is SUM divided by DENOM. "
                      "Neither is a percentile or an ancestry adjustment.\n\n")
        out.write("Match % describes allele matching to the scoring file, not the proportion of this sample's genotypes "
                  "that passed quality checks. Passing the matching threshold alone does not validate a score.\n\n")
        if not all(row.get("missing_genotypes") not in (None, "", "NA", "nan") for row in prs_rows):
            out.write("Per-score genotype missingness was not reported for every score; do not interpret an absent count as zero.\n\n")
        if run_ancestry:
            out.write("Ancestry normalization was requested, but this table contains raw scores only. "
                      "Consult the pgsc_calc HTML report for any reference-calibrated results and their reference population.\n\n")
        else:
            out.write("No reference calibration was requested. A reference-relative percentile or Z-score cannot be inferred "
                      "from these raw sums.\n\n")
    else:
        out.write("No final PRS scores are available.\n\n")
        if prs_stage.get("message"):
            out.write(f"{prs_stage['message']}\n\n")
    failures = summary.get("overlap_failures") or (prs_stage.get("diagnostics") or {}).get("overlap_failures") or []
    if failures:
        out.write("Scores excluded by the matching threshold:\n\n| Score | Match % |\n| --- | ---: |\n")
        for item in failures:
            out.write(f"| {markdown_cell(item['score'])} | {markdown_cell(item['match_percent'])} |\n")
        out.write("\n")


def render_report(
    context: dict[str, Any], monogenic_stage: dict[str, Any], prs_stage: dict[str, Any],
    *, output_summary: pathlib.Path,
) -> None:
    run_dir = pathlib.Path(context["results"]["run_dir"])
    report_dir = ensure_dir(pathlib.Path(context["results"]["report_dir"]))
    published = publish_stage_artifacts(run_dir=run_dir, monogenic_stage=monogenic_stage, prs_stage=prs_stage)
    artifact_manifest_path = report_dir / "artifact_manifest.json"
    write_json(artifact_manifest_path, {"context": context, "monogenic_stage": monogenic_stage, "prs_stage": prs_stage})
    mode = context.get("execution_mode") or derive_execution_mode(context["manifest"])
    title = {"monogenic-only": "Monogenic Report", "prs-only": "PRS Report", "both": "Genome Report"}[mode]
    enabled = {stage: bool(context["manifest"]["analysis"][stage]["enabled"]) for stage in ("monogenic", "prs")}
    stages = {"monogenic": monogenic_stage, "prs": prs_stage}
    report_path = report_dir / "report.md"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=report_dir, delete=False) as out:
            temporary = pathlib.Path(out.name)
            out.write(f"# {title}: {context['sample_id']}\n\nResearch-use results; this report does not provide clinical diagnoses.\n\n")
            if enabled["prs"]:
                render_prs_results(out, prs_stage, read_prs_rows((prs_stage.get("artifacts") or {}).get("aggregated_scores")),
                                   run_ancestry=bool(context["manifest"]["analysis"]["prs"].get("run_ancestry")))
            if enabled["monogenic"]:
                out.write("## Monogenic\n\n")
                summary = monogenic_stage.get("summary") or {}
                if monogenic_stage.get("status") == "ok":
                    out.write(f"Priority rows after deduplication: {summary.get('priority_rows_after_dedup', 'not reported')}. "
                              f"Variant records scanned: {summary.get('total_variant_records_seen', 'not reported')}.\n\n")
                elif monogenic_stage.get("message"):
                    out.write(markdown_cell(monogenic_stage["message"]) + "\n\n")
            out.write("## Execution Details\n\n")
            out.write(f"- Execution mode: `{mode}`\n")
            reference = context.get("reference_bundle") or {}
            out.write(f"- Reference bundle: `{reference.get('status', 'unknown')}`\n")
            for stage, active in enabled.items():
                if active:
                    out.write(f"- {stage.upper() if stage == 'prs' else 'Monogenic'} stage: `{stages[stage].get('status', 'unknown')}`\n")
            out.write("\n## Files\n\n")
            for stage, files in published.items():
                for key, path in files.items():
                    out.write(f"- [{stage}: {key.replace('_', ' ')}]({path})\n")
            out.write(f"- [Full provenance and canonical genotype/checkpoint paths]({artifact_manifest_path})\n\n")
            out.write("Genotype VCFs remain in the resumable cache. The provenance file records their locations and all stage diagnostics.\n")
        temporary.replace(report_path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    stage_statuses = {stage: payload.get("status") for stage, payload in stages.items()}
    failed = [stage for stage, active in enabled.items() if active and stage_statuses[stage] != "ok"]
    write_json(output_summary, {
        "sample_id": context["sample_id"], "execution_mode": mode,
        "report_path": str(report_path), "artifact_manifest_path": str(artifact_manifest_path),
        "enabled_stages": enabled, "stage_statuses": stage_statuses,
        "failed_enabled_stages": failed, "overall_status": "failed" if failed else "ok",
    })


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
            output_summary=pathlib.Path(args.output),
        )
        return
    if args.command == "enforce":
        raise SystemExit(enforce_outcome(load_json(args.summary)))


if __name__ == "__main__":
    main()
