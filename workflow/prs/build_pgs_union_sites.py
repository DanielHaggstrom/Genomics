#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pathlib
import urllib.request

from workflow.core.pipeline_lib import ensure_dir, hash_strings, load_json, run_command, write_json


class IndexedFasta:
    def __init__(self, fasta_path: pathlib.Path, fai_path: pathlib.Path) -> None:
        self.handle = open(fasta_path, "rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with open(fai_path, "r", encoding="utf-8") as handle:
            for raw in handle:
                contig, length, offset, line_bases, line_width = raw.rstrip("\n").split("\t")[:5]
                self.index[contig] = (int(length), int(offset), int(line_bases), int(line_width))

    def get_base(self, contig: str, pos: int) -> str:
        length, offset, line_bases, line_width = self.index[contig]
        if pos < 1 or pos > length:
            raise ValueError(f"Position out of range for {contig}:{pos}")
        zero = pos - 1
        byte_offset = offset + (zero // line_bases) * line_width + (zero % line_bases)
        self.handle.seek(byte_offset)
        return self.handle.read(1).decode("ascii").upper()


def fetch_scorefile_url(pgs_id: str) -> str:
    with urllib.request.urlopen(f"https://www.pgscatalog.org/rest/score/{pgs_id}", timeout=60) as response:
        payload = json.load(response)
    return payload["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]


def download_scorefile(url: str, destination: pathlib.Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, destination)


def normalize_contig(raw: str, prefix: str) -> str:
    raw = raw.strip()
    if prefix == "chr" and not raw.startswith("chr"):
        return f"chr{raw}"
    if not prefix and raw.startswith("chr"):
        return raw[3:]
    return raw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--reference-stage", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    reference_stage = load_json(args.reference_stage)
    pgs_ids = context["manifest"]["analysis"]["prs"]["pgs_ids"]
    score_cache_root = pathlib.Path(context["cache_dir"]) / "references" / "pgs_scorefiles" / "GRCh38"
    ids_hash = hash_strings(pgs_ids)[:16]
    stage_root = ensure_dir(pathlib.Path(context["prs_cache_root"]) / ids_hash / "targets")
    stage_manifest_path = stage_root / "union_stage.json"

    if stage_manifest_path.exists():
        payload = load_json(stage_manifest_path)
        if payload.get("pgs_ids") == pgs_ids and pathlib.Path(payload["artifacts"]["targets_vcf"]).exists():
            write_json(args.output, payload)
            return

    fasta = IndexedFasta(pathlib.Path(reference_stage["fasta"]), pathlib.Path(reference_stage["fai"]))
    prefix = context.get("contig_prefix", "")
    site_map: dict[tuple[str, int], dict[str, object]] = {}
    per_score_counts: dict[str, int] = {}
    skipped = 0

    for pgs_id in pgs_ids:
        url = fetch_scorefile_url(pgs_id)
        score_path = score_cache_root / pathlib.Path(url).name
        if not score_path.exists():
            download_scorefile(url, score_path)

        per_score_counts[pgs_id] = 0
        with gzip.open(score_path, "rt", encoding="utf-8", errors="replace") as handle:
            header = None
            for raw in handle:
                if raw.startswith("#"):
                    continue
                header = raw.rstrip("\n").split("\t")
                break

            if not header:
                raise RuntimeError(f"Could not read header from {score_path}")

            reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
            for row in reader:
                chrom_raw = row.get("hm_chr") or row.get("chr_name")
                pos_raw = row.get("hm_pos") or row.get("chr_position")
                effect = (row.get("effect_allele") or "").upper()
                other = (row.get("other_allele") or "").upper()
                if not chrom_raw or not pos_raw:
                    skipped += 1
                    continue
                if len(effect) != 1 or len(other) != 1 or effect not in {"A", "C", "G", "T"} or other not in {"A", "C", "G", "T"}:
                    skipped += 1
                    continue

                chrom = normalize_contig(chrom_raw, prefix)
                pos = int(pos_raw)
                try:
                    ref = fasta.get_base(chrom, pos)
                except Exception:
                    skipped += 1
                    continue

                alleles = {effect, other}
                if ref not in alleles:
                    skipped += 1
                    continue
                alts = alleles - {ref}
                if not alts:
                    skipped += 1
                    continue

                key = (chrom, pos)
                entry = site_map.setdefault(
                    key,
                    {
                        "ref": ref,
                        "alts": set(),
                        "pgs_ids": set(),
                    },
                )
                entry["alts"].update(alts)
                entry["pgs_ids"].add(pgs_id)
                per_score_counts[pgs_id] += 1

    plain_vcf = stage_root / "pgs_union_sites.vcf"
    with open(plain_vcf, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, meta in fasta.index.items():
            out.write(f"##contig=<ID={contig},length={meta[0]}>\n")
        out.write("##INFO=<ID=PGS_IDS,Number=.,Type=String,Description=\"PGS IDs contributing to this union site\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (chrom, pos) in sorted(site_map.keys(), key=lambda item: (item[0], item[1])):
            entry = site_map[(chrom, pos)]
            alt = ",".join(sorted(entry["alts"]))
            info = "PGS_IDS=" + ",".join(sorted(entry["pgs_ids"]))
            out.write(f"{chrom}\t{pos}\t.\t{entry['ref']}\t{alt}\t.\tPASS\t{info}\n")

    targets_vcf = stage_root / "pgs_union_sites.vcf.gz"
    run_command(["bcftools", "view", "-Oz", "-o", str(targets_vcf), str(plain_vcf)], check=True)
    run_command(["bcftools", "index", "-t", str(targets_vcf)], check=True)

    payload = {
        "stage": "prs_union_sites",
        "status": "ok",
        "pgs_ids": pgs_ids,
        "score_hash": ids_hash,
        "summary": {
            "union_sites": len(site_map),
            "skipped_rows": skipped,
            "per_score_counts": per_score_counts,
        },
        "artifacts": {
            "targets_vcf": str(targets_vcf),
            "targets_vcf_index": str(targets_vcf) + ".tbi",
            "plain_vcf": str(plain_vcf),
        },
    }
    write_json(stage_manifest_path, payload)
    write_json(args.output, payload)


if __name__ == "__main__":
    main()
