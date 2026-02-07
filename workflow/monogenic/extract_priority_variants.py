#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import defaultdict


LOF_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
}

CLIN_SIG_TERMS = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
}


def parse_info(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-txt", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_txt), exist_ok=True)

    csq_fields = []
    total_records = 0
    priority_rows = []
    reason_counts = defaultdict(int)

    with gzip.open(args.input_vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                marker = "Format: "
                idx = line.find(marker)
                if idx != -1:
                    fmt = line[idx + len(marker):].rstrip('">')
                    csq_fields = fmt.split("|")
                continue
            if line.startswith("#CHROM"):
                continue
            if line.startswith("#"):
                continue

            total_records += 1
            cols = line.split("\t")
            if len(cols) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample = cols[:10]
            info_map = parse_info(info)
            csq_blob = info_map.get("CSQ", "")
            if not csq_blob or not csq_fields:
                continue

            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            sample_map = dict(zip(fmt_keys, sample_vals))
            gt = sample_map.get("GT", ".")
            dp = sample_map.get("DP", ".")
            gq = sample_map.get("GQ", ".")
            ad = sample_map.get("AD", ".")

            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                ann_map = dict(zip(csq_fields, parts))

                consequence = ann_map.get("Consequence", "")
                clin_sig = ann_map.get("CLIN_SIG", "")
                symbol = ann_map.get("SYMBOL", "")
                gene = ann_map.get("Gene", "")
                impact = ann_map.get("IMPACT", "")
                canonical = ann_map.get("CANONICAL", "")
                existing_var = ann_map.get("Existing_variation", "")
                hgvsc = ann_map.get("HGVSc", "")
                hgvsp = ann_map.get("HGVSp", "")

                cons_terms = {c.strip() for c in consequence.split("&") if c.strip()}
                clin_terms = {c.strip().lower() for c in clin_sig.replace(" ", "_").split("&") if c.strip()}

                lof_hit = any(c in LOF_TERMS for c in cons_terms)
                clin_hit = any(c in CLIN_SIG_TERMS for c in clin_terms)
                if not (lof_hit or clin_hit):
                    continue

                reasons = []
                if clin_hit:
                    reasons.append("clinvar_pathogenic_or_likely_pathogenic")
                    reason_counts["clinvar_pathogenic_or_likely_pathogenic"] += 1
                if lof_hit:
                    reasons.append("predicted_loss_of_function")
                    reason_counts["predicted_loss_of_function"] += 1

                priority_rows.append(
                    {
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": vid,
                        "REF": ref,
                        "ALT": alt,
                        "FILTER": flt,
                        "GT": gt,
                        "DP": dp,
                        "GQ": gq,
                        "AD": ad,
                        "GENE_SYMBOL": symbol,
                        "GENE_ID": gene,
                        "CONSEQUENCE": consequence,
                        "IMPACT": impact,
                        "CLIN_SIG": clin_sig,
                        "EXISTING_VARIATION": existing_var,
                        "HGVSc": hgvsc,
                        "HGVSp": hgvsp,
                        "CANONICAL": canonical,
                        "REASONS": ",".join(reasons),
                    }
                )

    # Deduplicate exact rows while preserving order.
    seen = set()
    dedup = []
    for row in priority_rows:
        key = tuple(row.values())
        if key in seen:
            continue
        seen.add(key)
        dedup.append(row)

    header = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "FILTER",
        "GT",
        "DP",
        "GQ",
        "AD",
        "GENE_SYMBOL",
        "GENE_ID",
        "CONSEQUENCE",
        "IMPACT",
        "CLIN_SIG",
        "EXISTING_VARIATION",
        "HGVSc",
        "HGVSp",
        "CANONICAL",
        "REASONS",
    ]

    with open(args.output_tsv, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(header) + "\n")
        for row in dedup:
            out.write("\t".join(row[h] for h in header) + "\n")

    with open(args.summary_txt, "w", encoding="utf-8") as out:
        out.write(f"input_vcf\t{args.input_vcf}\n")
        out.write(f"total_variant_records_seen\t{total_records}\n")
        out.write(f"priority_rows_before_dedup\t{len(priority_rows)}\n")
        out.write(f"priority_rows_after_dedup\t{len(dedup)}\n")
        for k in sorted(reason_counts.keys()):
            out.write(f"count_{k}\t{reason_counts[k]}\n")


if __name__ == "__main__":
    main()
