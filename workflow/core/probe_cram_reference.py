#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib

from workflow.core.pipeline_lib import ensure_dir, load_json, run_command, stage_payload, write_json
from workflow.core.reference_lib import (
    extract_reference_hints,
    parse_program_entries,
    parse_read_group_entries,
    parse_sq_entries,
    summarize_contigs,
)


def run_remote_decode_probe(context: dict, first_contig: str, log_path: pathlib.Path) -> dict:
    cache_root = pathlib.Path(context["cache_dir"]) / "references" / "remote_md5_cache"
    ensure_dir(cache_root)
    region = f"{first_contig}:1-100000"
    env = {
        "REF_CACHE": str(cache_root / "%2s" / "%2s" / "%s"),
        "REF_PATH": f"{cache_root}/%2s/%2s/%s:https://www.ebi.ac.uk/ena/cram/md5/%s",
    }
    completed = run_command(
        [
            "bash",
            "-lc",
            'set -euo pipefail; samtools view "$1" "$2" | sed -n "1p" >/dev/null',
            "_",
            context["normalized_inputs"]["cram"],
            region,
        ],
        cwd=context["repo_root"],
        env=env,
        log_path=log_path,
        check=False,
    )
    if completed.returncode == 0:
        return {
            "attempted": True,
            "status": "ok",
            "region": region,
            "cache_dir": str(cache_root),
            "message": "samtools successfully decoded CRAM slices via REF_PATH/ENA MD5 resolution.",
            "log_path": str(log_path),
        }
    return {
        "attempted": True,
        "status": "failed",
        "region": region,
        "cache_dir": str(cache_root),
        "message": "samtools could not decode CRAM slices via REF_PATH/ENA MD5 resolution.",
        "log_path": str(log_path),
        "exit_code": completed.returncode,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    probe_dir = ensure_dir(pathlib.Path(context["sample_cache_dir"]) / "reference_probe")

    try:
        header_completed = run_command(
            ["samtools", "view", "-H", context["normalized_inputs"]["cram"]],
            cwd=context["repo_root"],
            log_path=probe_dir / "cram_header.log",
        )
        header_text = header_completed.stdout
        sq_entries = parse_sq_entries(header_text)
        pg_entries = parse_program_entries(header_text)
        rg_entries = parse_read_group_entries(header_text)
        contig_summary = summarize_contigs(sq_entries)
        reference_hints = extract_reference_hints(pg_entries)

        remote_probe = {
            "attempted": False,
            "status": "not_attempted",
            "message": "Remote MD5 decode probe disabled in manifest.",
        }
        if context["manifest"]["reference"].get("remote_decode_probe"):
            first_contig = contig_summary["first_contigs"][0] if contig_summary["first_contigs"] else None
            if first_contig:
                remote_probe = run_remote_decode_probe(context, first_contig, probe_dir / "remote_decode_probe.log")
            else:
                remote_probe = {
                    "attempted": False,
                    "status": "not_attempted",
                    "message": "No CRAM contigs were parsed from the header, so remote decode probe was skipped.",
                }

        payload = stage_payload(
            "reference_probe",
            "ok",
            summary={
                "contig_count": contig_summary["count"],
                "first_contigs": contig_summary["first_contigs"],
                "contig_prefix": contig_summary["contig_prefix"],
                "all_contigs_have_md5": contig_summary["all_contigs_have_md5"],
                "program_count": len(pg_entries),
                "read_group_count": len(rg_entries),
            },
            cram_header={
                "sq_entries": sq_entries,
                "pg_entries": pg_entries,
                "rg_entries": rg_entries,
            },
            provenance_hints=reference_hints,
            remote_decode_probe=remote_probe,
            artifacts={
                "header_log": str(probe_dir / "cram_header.log"),
                "remote_decode_log": remote_probe.get("log_path"),
            },
        )
        write_json(args.output, payload)
    except Exception as exc:
        payload = stage_payload("reference_probe", "failed", message=str(exc))
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
