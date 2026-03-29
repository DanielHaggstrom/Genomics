#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import re
from typing import Any


def parse_tagged_header_line(line: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for token in line.rstrip().split("\t")[1:]:
        if ":" not in token:
            continue
        key, value = token.split(":", 1)
        fields[key] = value.strip()
    return fields


def parse_sq_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "name": fields.get("SN"),
                "length": int(fields["LN"]) if fields.get("LN") else None,
                "md5": fields.get("M5"),
                "uri": fields.get("UR"),
            }
        )
    return entries


def parse_program_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@PG"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "id": fields.get("ID"),
                "name": fields.get("PN"),
                "version": fields.get("VN"),
                "command_line": fields.get("CL"),
                "previous": fields.get("PP"),
                "description": fields.get("DS"),
            }
        )
    return entries


def parse_read_group_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@RG"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "id": fields.get("ID"),
                "sample": fields.get("SM"),
                "platform_unit": fields.get("PU"),
                "library": fields.get("LB"),
            }
        )
    return entries


def summarize_contigs(sq_entries: list[dict[str, Any]]) -> dict[str, Any]:
    names = [entry["name"] for entry in sq_entries if entry.get("name")]
    primary = [name for name in names if name in {*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY", "chrM", *(str(i) for i in range(1, 23)), "X", "Y", "MT", "M"}]
    has_md5_for_all = all(entry.get("md5") for entry in sq_entries if entry.get("name"))
    contig_prefix = ""
    if names:
        contig_prefix = "chr" if names[0].startswith("chr") else ""
    return {
        "count": len(names),
        "primary_contigs_present": primary,
        "first_contigs": names[:10],
        "contig_prefix": contig_prefix,
        "all_contigs_have_md5": has_md5_for_all,
    }


def extract_reference_hints(pg_entries: list[dict[str, Any]]) -> dict[str, Any]:
    hints: dict[str, Any] = {
        "aligner": None,
        "ht_reference": None,
        "ref_dir": None,
        "variant_annotation_assembly": None,
    }
    command_lines = [entry.get("command_line") or "" for entry in pg_entries]
    combined = "\n".join(command_lines)
    lowered = combined.lower()

    if "dragen" in lowered:
        hints["aligner"] = "DRAGEN"

    patterns = {
        "ht_reference": [
            r"--ht-reference=(\S+)",
            r"--ht-reference\s+(\S+)",
        ],
        "ref_dir": [
            r"--ref-dir=(\S+)",
            r"--ref-dir\s+(\S+)",
        ],
        "variant_annotation_assembly": [
            r"--variant-annotation-assembly=(\S+)",
            r"--variant-annotation-assembly\s+(\S+)",
        ],
    }

    for key, regexes in patterns.items():
        for regex in regexes:
            match = re.search(regex, combined)
            if match:
                hints[key] = match.group(1)
                break

    return hints


def parse_sequence_dict(path: str | pathlib.Path) -> list[dict[str, Any]]:
    text = pathlib.Path(path).read_text(encoding="utf-8")
    return parse_sq_entries(text)


def index_sq_entries(entries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {entry["name"]: entry for entry in entries if entry.get("name")}


def compare_sequence_catalogs(cram_sq_entries: list[dict[str, Any]], candidate_sq_entries: list[dict[str, Any]]) -> dict[str, Any]:
    cram_map = index_sq_entries(cram_sq_entries)
    candidate_map = index_sq_entries(candidate_sq_entries)

    missing_contigs: list[str] = []
    length_mismatches: list[dict[str, Any]] = []
    md5_mismatches: list[dict[str, Any]] = []
    md5_missing_in_candidate: list[str] = []

    for name, cram_entry in cram_map.items():
        candidate_entry = candidate_map.get(name)
        if candidate_entry is None:
            missing_contigs.append(name)
            continue

        if candidate_entry.get("length") != cram_entry.get("length"):
            length_mismatches.append(
                {
                    "contig": name,
                    "cram_length": cram_entry.get("length"),
                    "candidate_length": candidate_entry.get("length"),
                }
            )

        cram_md5 = cram_entry.get("md5")
        candidate_md5 = candidate_entry.get("md5")
        if cram_md5 and not candidate_md5:
            md5_missing_in_candidate.append(name)
        elif cram_md5 and candidate_md5 and cram_md5 != candidate_md5:
            md5_mismatches.append(
                {
                    "contig": name,
                    "cram_md5": cram_md5,
                    "candidate_md5": candidate_md5,
                }
            )

    extra_contigs = sorted(name for name in candidate_map.keys() if name not in cram_map)
    compared_contigs = len(cram_map)
    exact_for_cram = not missing_contigs and not length_mismatches and not md5_mismatches and not md5_missing_in_candidate
    exact_dictionary_match = exact_for_cram and not extra_contigs and len(candidate_map) == len(cram_map)

    return {
        "compared_contigs": compared_contigs,
        "candidate_contigs": len(candidate_map),
        "exact_for_cram": exact_for_cram,
        "exact_dictionary_match": exact_dictionary_match,
        "missing_contigs": missing_contigs,
        "length_mismatches": length_mismatches,
        "md5_mismatches": md5_mismatches,
        "md5_missing_in_candidate": md5_missing_in_candidate,
        "extra_contigs": extra_contigs,
    }
