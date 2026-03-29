#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import re

from workflow.core.pipeline_lib import ensure_dir, load_json, run_command, stage_payload, write_json


OVERLAP_RE = re.compile(r"Score (\S+) fails minimum matching threshold \(([\d.]+)% variants match\)")


def parse_overlap_diagnostics(log_text: str) -> list[dict[str, str]]:
    matches = []
    for line in log_text.splitlines():
        found = OVERLAP_RE.search(line)
        if not found:
            continue
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


def find_first(root: pathlib.Path, pattern: str) -> str | None:
    matches = sorted(root.rglob(pattern))
    return str(matches[0]) if matches else None


def maybe_reuse(stage_manifest_path: pathlib.Path, cache_key: str) -> dict | None:
    if not stage_manifest_path.exists():
        return None
    payload = load_json(stage_manifest_path)
    if payload.get("cache_key") != cache_key or payload.get("status") != "ok":
        return None
    aggregated = (payload.get("artifacts") or {}).get("aggregated_scores")
    if not aggregated or not pathlib.Path(aggregated).exists():
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
    prs_cfg = context["manifest"]["analysis"]["prs"]
    stage_root = ensure_dir(pathlib.Path(context["prs_cache_root"]))
    ids = prs_cfg.get("pgs_ids") or []
    cache_key = f"{context['input_fingerprint']['id']}:{'|'.join(ids)}:{prs_cfg['min_overlap']}"
    stage_manifest_path = stage_root / "stage_manifest.json"

    try:
        if not prs_cfg["enabled"]:
            payload = stage_payload("prs", "skipped", message="PRS stage disabled.")
            write_json(args.output, payload)
            return

        reused = maybe_reuse(stage_manifest_path, cache_key)
        if reused:
            write_json(args.output, reused)
            return

        if reference_stage.get("status") != "ready":
            payload = stage_payload(
                "prs",
                "failed",
                cache_key=cache_key,
                message="PRS targeted genotyping requires a ready reference bundle that matches the CRAM.",
            )
            write_json(stage_manifest_path, payload)
            write_json(args.output, payload)
            return

        union_json = stage_root / "union_stage.json"
        run_command(
            [
                "python3",
                "-m",
                "workflow.prs.build_pgs_union_sites",
                "--context",
                args.context,
                "--reference-stage",
                args.reference_stage,
                "--output",
                str(union_json),
            ],
            cwd=context["repo_root"],
            log_path=stage_root / "build_union.log",
        )
        union_stage = load_json(union_json)

        target_dir = ensure_dir(stage_root / union_stage["score_hash"] / "targeted_genotypes")
        targeted_vcf = target_dir / f"{context['sample_id']}.prs_targeted.vcf.gz"
        if not targeted_vcf.exists():
            run_command(
                [
                    "bash",
                    "workflow/prs/targeted_call_from_cram.sh",
                    "--cram",
                    context["normalized_inputs"]["cram"],
                    "--reference-fasta",
                    reference_stage["fasta"],
                    "--targets-vcf",
                    union_stage["artifacts"]["targets_vcf"],
                    "--outdir",
                    str(target_dir),
                    "--sample-id",
                    context["sample_id"],
                ],
                cwd=context["repo_root"],
                log_path=target_dir / "targeted_call.log",
            )

        samplesheet = stage_root / "samplesheet.csv"
        pgs_ids_file = stage_root / "pgs_ids.txt"
        pgsc_outdir = stage_root / union_stage["score_hash"] / "pgsc_calc"
        pgsc_workdir = stage_root / union_stage["score_hash"] / "work"
        genotypes_cache = stage_root / union_stage["score_hash"] / "genotypes_cache"
        ensure_dir(pgsc_outdir)
        ensure_dir(pgsc_workdir)
        ensure_dir(genotypes_cache)

        samplesheet.write_text(
            "sampleset,path_prefix,chrom,format,vcf_genotype_field\n"
            f"{safe_sampleset_name(context['sample_id'])},{str(targeted_vcf)[:-7]},,vcf,GT\n",
            encoding="utf-8",
        )
        pgs_ids_file.write_text("\n".join(ids) + "\n", encoding="utf-8")

        log_path = stage_root / "run_pgsc_calc.log"
        completed = run_command(
            [
                "bash",
                "workflow/prs/run_pgsc_calc.sh",
                "--samplesheet",
                str(samplesheet),
                "--pgs-ids-file",
                str(pgs_ids_file),
                "--outdir",
                str(pgsc_outdir),
                "--workdir",
                str(pgsc_workdir),
                "--profile",
                context["manifest"]["runtime"]["profile"],
                "--max-memory",
                context["manifest"]["runtime"]["max_memory"],
                "--min-overlap",
                str(prs_cfg["min_overlap"]),
                "--resume",
                "1",
                "--genotypes-cache",
                str(genotypes_cache),
            ],
            cwd=context["repo_root"],
            log_path=log_path,
            check=False,
        )

        log_text = pathlib.Path(log_path).read_text(encoding="utf-8")
        aggregated_scores = find_first(pgsc_outdir, "aggregated_scores.txt.gz")
        pgsc_report = find_first(pgsc_outdir, "report.html")
        match_summary = find_first(pgsc_outdir, "*summary*.csv")
        overlap_diagnostics = parse_overlap_diagnostics(log_text)

        if completed.returncode != 0 or not aggregated_scores:
            payload = stage_payload(
                "prs",
                "failed",
                cache_key=cache_key,
                message="pgsc_calc did not produce final PRS scores.",
                diagnostics={
                    "exit_code": completed.returncode,
                    "overlap_failures": overlap_diagnostics,
                    "log_path": str(log_path),
                },
                artifacts={
                    "targets_vcf": str(targeted_vcf),
                    "pgsc_outdir": str(pgsc_outdir),
                    "log_path": str(log_path),
                },
            )
            write_json(stage_manifest_path, payload)
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
            reused=False,
            summary={
                "num_requested_scores": len(ids),
                "preview_rows": preview_rows,
                "overlap_failures": overlap_diagnostics,
            },
            artifacts={
                "targets_vcf": str(targeted_vcf),
                "targets_vcf_index": str(targeted_vcf) + ".tbi",
                "aggregated_scores": aggregated_scores,
                "pgsc_report_html": pgsc_report,
                "match_summary": match_summary,
                "pgsc_outdir": str(pgsc_outdir),
                "pgsc_log": str(log_path),
            },
        )
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)
    except Exception as exc:
        payload = stage_payload("prs", "failed", cache_key=cache_key, message=str(exc))
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
