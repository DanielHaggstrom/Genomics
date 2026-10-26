"""Materialize allele-defined PRS sites from finalized, normalized single-sample calls.

This is not a variant caller or a gVCF decoder. GenotypeGVCFs must finalize the
input first. Only explicit, quality-controlled calls are retained. In particular,
an absent record is never interpreted as reference, and a reference SNP call does
not establish the absence of an indel. Both inputs must use the same reference,
contig names, normalization and coordinate order.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from itertools import groupby
import json
import math
from pathlib import Path
import re
from typing import TextIO

from workflow.core.runtime import check_cancelled


MATERIALIZER_VERSION = 1
_BASES = re.compile(r"[ACGT]+\Z")
_CONTIG = re.compile(r"##contig=<ID=([^,>]+)")


def _header(handle: TextIO, *, sample: bool) -> tuple[list[str], list[str]]:
    metadata = []
    for line in handle:
        if line.startswith("##"):
            metadata.append(line.rstrip())
        elif line.startswith("#CHROM\t"):
            columns = line.rstrip().split("\t")
            if len(columns) != (10 if sample else 8):
                raise ValueError("Expected one sample in calls and a sites-only target VCF")
            return metadata, columns
        else:
            raise ValueError("Invalid VCF header")
    raise ValueError("Missing VCF column header")


def _records(handle: TextIO, ranks: dict[str, int], *, sample: bool):
    previous = None
    for index, line in enumerate(handle):
        if index % 1024 == 0:
            check_cancelled()
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != (10 if sample else 8):
            raise ValueError("Invalid VCF record column count")
        chrom, position = fields[0], int(fields[1])
        if chrom not in ranks or position < 1:
            raise ValueError("VCF record has an undeclared contig or invalid position")
        key = (ranks[chrom], position)
        if previous is not None and key < previous:
            raise ValueError("VCF records must be sorted in contig-header order")
        previous = key
        if sample:
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            if any(alt in {"<NON_REF>", "<*>"} for alt in fields[4].split(",")) or int(info.get("END", position)) > position:
                raise ValueError("Calls must be finalized all-sites VCF, not gVCF blocks or NON_REF alleles")
        yield key, fields


def _output_header(output: TextIO, contigs: dict[str, str], sample_id: str, *, thresholds: tuple[int, int] | None = None) -> None:
    if not sample_id or re.search(r"\s", sample_id):
        raise ValueError("Sample ID must be nonempty and contain no whitespace")
    output.write("##fileformat=VCFv4.2\n")
    for line in contigs.values():
        output.write(line + "\n")
    output.write(f"##source=prs_genotypes_v{MATERIALIZER_VERSION}\n")
    if thresholds:
        output.write(f"##prs_min_depth={thresholds[0]}\n##prs_min_gq={thresholds[1]}\n")
    for name, kind, description in [
        ("GT", "String", "Genotype on the target REF/ALT alleles"),
        ("DS", "Float", "Alternate allele dosage; direct calls take precedence"),
        ("DP", "Integer", "Original per-sample read depth"),
        ("GQ", "Integer", "Original genotype quality"),
        ("RGQ", "Integer", "Original reference genotype quality"),
        ("FT", "String", "PASS or reason this target genotype is missing"),
    ]:
        output.write(f'##FORMAT=<ID={name},Number=1,Type={kind},Description="{description}">\n')
    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")


def _call(target: list[str], calls: list[list[str]], min_depth: int, min_gq: int) -> tuple[str, str, dict[str, str]]:
    ref, alt = target[3:5]
    if not _BASES.fullmatch(ref) or not _BASES.fullmatch(alt) or ref == alt:
        return "./.", "unsupported_target", {}
    if not calls:
        return "./.", "no_record", {}
    if len(calls) != 1:
        return "./.", "duplicate_call_position", {}
    call = calls[0]
    values = dict(zip(call[8].split(":"), call[9].split(":")))
    evidence = {key: values.get(key, ".") for key in ("DP", "GQ", "RGQ")}
    if call[6] not in {".", "PASS"} or values.get("FT", ".") not in {".", "PASS"}:
        return "./.", "filtered_call", evidence
    if call[3] != ref:
        return "./.", "reference_mismatch", evidence
    gt = values.get("GT", ".")
    indices = re.split(r"[/|]", gt)
    if any(index == "." for index in indices):
        return "./.", "no_call", evidence
    if len(indices) != 2:
        return "./.", "unsupported_ploidy", evidence
    alleles = [call[3], *call[4].split(",")]
    try:
        if any(not index.isdigit() or int(index) >= len(alleles) for index in indices):
            raise ValueError
        genotype = [alleles[int(index)] for index in indices]
    except ValueError as exc:
        raise ValueError("Invalid genotype allele index") from exc
    try:
        depth = int(evidence["DP"])
        quality_key = "RGQ" if indices == ["0", "0"] and evidence["RGQ"] != "." else "GQ"
        quality = int(evidence[quality_key])
    except ValueError:
        return "./.", "missing_quality", evidence
    if depth < min_depth:
        return "./.", "low_depth", evidence
    if quality < min_gq:
        return "./.", "low_genotype_quality", evidence
    if any(allele not in {ref, alt} for allele in genotype):
        return "./.", "other_allele", evidence
    # A concrete variant call establishes this indel; an uninformative anchor
    # base alone does not. SNP reference calls can safely retain scaffold ALT.
    if alt not in alleles and (len(ref) != 1 or len(alt) != 1):
        return "./.", "unresolved_indel", evidence
    mapped = ["0" if allele == ref else "1" for allele in genotype]
    separator = "|" if "|" in gt else "/"
    return separator.join(mapped), "PASS", evidence


def materialize_scaffold(
    calls: TextIO,
    targets: TextIO,
    output: TextIO,
    *,
    sample_id: str | None = None,
    min_depth: int = 10,
    min_gq: int = 20,
) -> dict:
    """Stream targets and calls into a VCF with explicit alleles and callability.

    Uncallable rows keep their target alleles and missing GT, plus FORMAT/FT.
    Consumers must exclude these rows before scoring; site overlap is not a
    substitute for sample callability. Thresholds are research QC defaults.
    """
    if min_depth < 1 or min_gq < 1:
        raise ValueError("Depth and genotype-quality thresholds must be positive")
    call_metadata, call_columns = _header(calls, sample=True)
    target_metadata, _ = _header(targets, sample=False)
    contigs = {}
    for line in [*call_metadata, *target_metadata]:
        match = _CONTIG.match(line)
        if match:
            contigs.setdefault(match[1], line)
    ranks = {chrom: index for index, chrom in enumerate(contigs)}
    _output_header(output, contigs, sample_id or call_columns[9], thresholds=(min_depth, min_gq))
    grouped_calls = groupby(_records(calls, ranks, sample=True), key=lambda item: item[0])
    current = next(grouped_calls, None)
    current_rows = [row for _, row in current[1]] if current else []
    counts: Counter[str] = Counter()
    total = 0
    for key, group in groupby(_records(targets, ranks, sample=False), key=lambda item: item[0]):
        while current is not None and current[0] < key:
            current = next(grouped_calls, None)
            current_rows = [row for _, row in current[1]] if current else []
        observed = current_rows if current is not None and current[0] == key else []
        seen = set()
        for _, target in group:
            for alt in target[4].split(","):
                allele_key = (target[3], alt)
                if allele_key in seen:
                    continue
                seen.add(allele_key)
                row = [*target[:4], alt, *target[5:]]
                gt, reason, evidence = _call(row, observed, min_depth, min_gq)
                counts[reason] += 1
                total += 1
                sample = ":".join([gt, *(evidence.get(field, ".") for field in ("DP", "GQ", "RGQ")), reason])
                output.write("\t".join([*row[:5], ".", ".", ".", "GT:DP:GQ:RGQ:FT", sample]) + "\n")
    # Validate the tail too, so malformed calls cannot be hidden beyond the last
    # requested target. groupby advances lazily; consume each group explicitly.
    for _, group in grouped_calls:
        for _ in group:
            pass
    return {
        "materializer_version": MATERIALIZER_VERSION,
        "total_sites": total,
        "callable_sites": counts["PASS"],
        "missing_sites": total - counts["PASS"],
        "reasons": dict(counts),
        "min_depth": min_depth,
        "min_gq": min_gq,
    }


def augment_scaffold(direct: TextIO, imputed: TextIO, output: TextIO, *, sample_id: str | None = None) -> dict:
    """Augment QC-filtered direct calls with normalized biallelic imputed calls.

    This preserves direct-only sites and exact direct genotypes. Callers must
    quality-filter imputed records (e.g. DR2) first; a DS value is not itself a
    quality criterion. Unsupported/uncalled records remain missing for auditing.
    """
    direct_metadata, direct_columns = _header(direct, sample=True)
    imputed_metadata, imputed_columns = _header(imputed, sample=True)
    if direct_columns[9] != imputed_columns[9]:
        raise ValueError("Direct and imputed VCFs must contain the same sample")
    contigs = {}
    for line in [*direct_metadata, *imputed_metadata]:
        match = _CONTIG.match(line)
        if match:
            contigs.setdefault(match[1], line)
    ranks = {chrom: index for index, chrom in enumerate(contigs)}
    _output_header(output, contigs, sample_id or direct_columns[9])
    def groups(handle: TextIO):
        for key, group in groupby(_records(handle, ranks, sample=True), key=lambda item: item[0]):
            rows = {}
            for _, row in group:
                alleles = (row[3], row[4])
                if not all(_BASES.fullmatch(allele) for allele in alleles) or row[3] == row[4]:
                    raise ValueError("Merge inputs must contain concrete biallelic REF/ALT records")
                if alleles in rows:
                    raise ValueError("Duplicate allele-defined record in merge input")
                rows[alleles] = row
            yield key, rows
    iterators = [iter(groups(direct)), iter(groups(imputed))]
    current = [next(iterator, None) for iterator in iterators]
    counts: Counter[str] = Counter()
    while any(item is not None for item in current):
        key = min(item[0] for item in current if item is not None)
        sources = [item[1] if item is not None and item[0] == key else {} for item in current]
        observed_position = False
        for row in sources[0].values():
            values = dict(zip(row[8].split(":"), row[9].split(":")))
            if row[6] in {".", "PASS"} and values.get("FT", ".") in {".", "PASS"} and re.fullmatch(r"[01][/|][01]", values.get("GT", ".")):
                observed_position = True
                break
        for alleles in sorted(set(sources[0]) | set(sources[1])):
            records = [source.get(alleles) for source in sources]
            chosen = next(row for row in records if row is not None)
            selected = {}
            source_name = "missing"
            for source_index, row in enumerate(records):
                # A different allele pair at an already observed position cannot
                # be filled with a contradictory population-panel genotype.
                if source_index == 1 and observed_position:
                    continue
                if row is None or row[6] not in {".", "PASS"}:
                    continue
                values = dict(zip(row[8].split(":"), row[9].split(":")))
                if values.get("FT", ".") not in {".", "PASS"}:
                    continue
                gt = values.get("GT", ".")
                if not re.fullmatch(r"[01][/|][01]", gt):
                    continue
                # Exact direct genotypes win even when imputation disagrees.
                dosage = float(gt[0] == "1") + float(gt[2] == "1")
                if source_index == 1 and values.get("DS", ".") != ".":
                    try:
                        dosage = float(values["DS"])
                    except ValueError:
                        continue
                    if not math.isfinite(dosage) or not 0 <= dosage <= 2:
                        continue
                selected = {**values, "GT": gt, "DS": f"{dosage:g}", "FT": "PASS"}
                chosen = row
                source_name = "direct" if source_index == 0 else "imputed"
                break
            counts[source_name] += 1
            if not selected:
                selected = dict(zip(chosen[8].split(":"), chosen[9].split(":")))
                selected.update(GT="./.", DS=".", FT="no_usable_call")
            fields = ("GT", "DS", "DP", "GQ", "RGQ", "FT")
            output.write("\t".join([*chosen[:5], ".", ".", ".", ":".join(fields), ":".join(selected.get(field, ".") for field in fields)]) + "\n")
        for index, item in enumerate(current):
            if item is not None and item[0] == key:
                current[index] = next(iterators[index], None)
    return {"total_sites": sum(counts.values()), "callable_sites": counts["direct"] + counts["imputed"], "missing_sites": counts["missing"], "sources": dict(counts)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", nargs="?", choices=("materialize", "merge"), default="materialize")
    parser.add_argument("--calls")
    parser.add_argument("--targets")
    parser.add_argument("--direct")
    parser.add_argument("--imputed")
    parser.add_argument("--output", required=True, help="Plain-text output VCF")
    parser.add_argument("--sample-id")
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--min-gq", type=int, default=20)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    def opener(path: str):
        return gzip.open(path, "rt", encoding="utf-8") if path.endswith((".gz", ".bgz")) else open(path, encoding="utf-8")
    required = (args.direct, args.imputed) if args.operation == "merge" else (args.calls, args.targets)
    if not all(required):
        parser.error("merge needs --direct/--imputed; materialize needs --calls/--targets")
    inputs = [Path(path) for path in required]
    destinations = [Path(args.output), Path(args.summary)]
    if destinations[0].resolve() == destinations[1].resolve():
        parser.error("VCF and summary outputs must be different files")
    for destination in destinations:
        if any(destination.resolve() == source.resolve() or
               (destination.exists() and destination.samefile(source)) for source in inputs):
            parser.error("Refusing to overwrite an input file or its alias")
    with opener(required[0]) as first, opener(required[1]) as second, open(args.output, "w", encoding="utf-8") as output:
        if args.operation == "merge":
            summary = augment_scaffold(first, second, output, sample_id=args.sample_id)
        else:
            summary = materialize_scaffold(first, second, output, sample_id=args.sample_id, min_depth=args.min_depth, min_gq=args.min_gq)
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
