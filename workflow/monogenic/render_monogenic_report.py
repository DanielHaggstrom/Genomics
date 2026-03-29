#!/usr/bin/env python3
import argparse
import csv
from collections import Counter


NON_MONOGENIC_TAGS = {
    "affects",
    "association",
    "drug_response",
    "established_risk_allele",
    "likely_risk_allele",
    "protective",
    "risk_factor",
}


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def read_tsv(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_summary_txt(path: str) -> dict[str, str]:
    stats = {}
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            stats[key] = value
    return stats


def variant_gene_key(row: dict) -> tuple[str, ...]:
    return (
        row.get("CHROM", ""),
        row.get("POS", ""),
        row.get("REF", ""),
        row.get("ALT", ""),
        row.get("GT", ""),
        row.get("GENE_SYMBOL", ""),
        row.get("GENE_ID", ""),
    )


def sort_review_row(row: dict) -> tuple:
    impact_rank = 0 if row.get("IMPACT") == "HIGH" else 1
    review_rank = {"high": 0, "single_submitter": 1, "other": 2, "none": 3}.get(row.get("REVIEW_TIER", ""), 4)
    conflict_rank = 0 if row.get("CLINVAR_CONFLICT_FLAG") == "NO" else 1
    af_rank = parse_float(row.get("MAX_AF", ""))
    if af_rank is None:
        af_rank = 999.0
    return (
        conflict_rank,
        impact_rank,
        review_rank,
        af_rank,
        row.get("GENE_SYMBOL", ""),
        row.get("POS", ""),
    )


def sort_broad_row(row: dict) -> tuple:
    id_rank = 0 if row.get("ID", ".") != "." else 1
    return (
        row.get("GENE_SYMBOL", ""),
        id_rank,
        row.get("ID", ""),
        row.get("POS", ""),
    )


def format_variant_label(row: dict) -> str:
    gene = row.get("GENE_SYMBOL") or row.get("GENE_ID") or "UNKNOWN"
    variant_id = row.get("ID", ".")
    consequence = row.get("CONSEQUENCE", "")
    hgvsp = row.get("HGVSp", "")
    pieces = [f"`{gene}`"]
    if variant_id and variant_id != ".":
        pieces.append(f"`{variant_id}`")
    pieces.append(f"`{row.get('GT', '.')}`")
    if hgvsp:
        pieces.append(f"`{hgvsp}`")
    elif consequence:
        pieces.append(f"`{consequence}`")
    return " ".join(pieces)


def short_disease_list(value: str, limit: int = 3) -> str:
    items = [x.replace("_", " ") for x in value.split("&") if x and x not in {"not_provided", "not_specified"}]
    if not items:
        return "not specified"
    trimmed = items[:limit]
    if len(items) > limit:
        trimmed.append("...")
    return ", ".join(trimmed)


def clin_sig_terms(value: str) -> set[str]:
    return {x.strip().lower() for x in value.replace(" ", "_").split("&") if x.strip()}


def is_notable_marker(row: dict) -> bool:
    terms = clin_sig_terms(row.get("CLIN_SIG", ""))
    return any(tag in terms for tag in NON_MONOGENIC_TAGS)


def write_table(out, headers: list[str], rows: list[list[str]]) -> None:
    if not rows:
        return
    out.write("| " + " | ".join(headers) + " |\n")
    out.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
    for row in rows:
        out.write("| " + " | ".join(row) + " |\n")
    out.write("\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output-md", required=True)
    p.add_argument("--priority-summary", required=True)
    p.add_argument("--clinvar-unique-tsv", required=True)
    p.add_argument("--lof-unique-tsv", required=True)
    p.add_argument("--review-unique-tsv", required=True)
    p.add_argument("--mode", required=True)
    p.add_argument("--max-af", required=True)
    args = p.parse_args()

    priority_stats = read_summary_txt(args.priority_summary)
    clinvar_rows = sorted(read_tsv(args.clinvar_unique_tsv), key=sort_broad_row)
    lof_rows = read_tsv(args.lof_unique_tsv)
    review_rows = sorted(read_tsv(args.review_unique_tsv), key=sort_review_row)

    review_keys = {variant_gene_key(row) for row in review_rows}
    reviewed_rows = [
        row for row in review_rows
        if row.get("CLINVAR_CONFLICT_FLAG") == "NO" and row.get("REVIEW_TIER") == "high"
    ]
    review_queue_rows = [row for row in review_rows if row not in reviewed_rows]
    notable_marker_rows = [
        row for row in clinvar_rows
        if variant_gene_key(row) not in review_keys and is_notable_marker(row)
    ]
    broad_only_count = len([row for row in clinvar_rows if variant_gene_key(row) not in review_keys])
    lof_gene_counts = Counter(
        (row.get("GENE_SYMBOL") or row.get("GENE_ID") or "UNKNOWN") for row in lof_rows
    )
    top_lof_genes = [
        (gene, count) for gene, count in lof_gene_counts.most_common(10)
        if gene != "UNKNOWN"
    ]
    overlap_rows = [row for row in clinvar_rows if variant_gene_key(row) in review_keys]

    with open(args.output_md, "w", encoding="utf-8") as out:
        out.write(f"# Monogenic Report: {args.sample_id}\n\n")
        out.write("This report is the primary human-readable summary for the monogenic stage. ")
        out.write("It favors unique variant-gene hits over transcript-level row counts and keeps noisy detail in TSV artifacts.\n\n")

        out.write("## Overview\n")
        out.write(f"- Mode: `{args.mode}`\n")
        out.write(f"- Rare-review MAX_AF threshold: `{args.max_af}`\n")
        if "total_variant_records_seen" in priority_stats:
            out.write(f"- VCF records scanned: `{priority_stats['total_variant_records_seen']}`\n")
        if "priority_rows_after_dedup" in priority_stats:
            out.write(f"- Priority transcript rows after dedup: `{priority_stats['priority_rows_after_dedup']}`\n")
        out.write(f"- Broad ClinVar/protein-altering unique hits: `{len(clinvar_rows)}`\n")
        out.write(f"- Predicted LoF/protein-altering unique hits: `{len(lof_rows)}`\n")
        out.write(f"- Rare ClinVar review unique hits: `{len(review_rows)}`\n")
        out.write(f"- Reviewed, non-conflicting rare ClinVar hits: `{len(reviewed_rows)}`\n")
        out.write(f"- Overlap between broad ClinVar hits and rare review hits: `{len(overlap_rows)}`\n\n")

        out.write("## How To Read This\n")
        out.write("- `Reviewed rare ClinVar findings` are the strongest disease-style leads in this run.\n")
        out.write("- `Notable non-monogenic markers` are usually common risk, drug-response, or association alleles. They may matter, but they are not monogenic diagnoses.\n")
        out.write("- `Rare ClinVar review queue` is a manual-review list, not a diagnosis list. Most rows here still have conflicting ClinVar interpretations.\n")
        out.write("- `Predicted LoF` is the noisiest bucket. Use it for exploration, not for standalone conclusions.\n\n")

        out.write("## Why Overlap Is Small\n")
        out.write(f"- The broad ClinVar bucket contains `{len(clinvar_rows)}` unique hits and accepts any protein-altering row whose VEP `CLIN_SIG` mentions pathogenic or likely pathogenic, even when benign, risk-factor, or drug-response labels are also present.\n")
        out.write(f"- The rare review bucket contains `{len(review_rows)}` unique hits and adds rarity plus ClinVar review metadata, but in this sample `{len(review_queue_rows)}` of those remain conflict-heavy review items.\n")
        if overlap_rows:
            overlap_genes = ", ".join(sorted({row.get('GENE_SYMBOL') or row.get('GENE_ID') or 'UNKNOWN' for row in overlap_rows}))
            out.write(f"- Only `{len(overlap_rows)}` unique hit overlaps both buckets in this run: {overlap_genes}.\n\n")
        else:
            out.write("- No unique hits overlap both buckets in this run.\n\n")

        out.write("## Reviewed Rare ClinVar Findings\n")
        if reviewed_rows:
            rows = []
            for row in reviewed_rows:
                rows.append([
                    format_variant_label(row),
                    row.get("MAX_AF", "") or "n/a",
                    short_disease_list(row.get("CLINVAR_CLNDN", "")),
                    row.get("CLINVAR_CLNREVSTAT", "") or "n/a",
                ])
            write_table(out, ["Variant", "MAX_AF", "ClinVar conditions", "Review"], rows)
        else:
            out.write("No reviewed, non-conflicting rare ClinVar findings passed the current filter.\n\n")

        out.write("## Notable Non-Monogenic Markers\n")
        if notable_marker_rows:
            rows = []
            for row in notable_marker_rows:
                rows.append([
                    format_variant_label(row),
                    row.get("CLIN_SIG", "") or "n/a",
                    row.get("REASONS", "") or "n/a",
                ])
            write_table(out, ["Variant", "ClinVar labels", "Why it landed in the broad bucket"], rows)
        else:
            out.write("No risk/drug-response/association-style markers stood out outside the rare review queue.\n\n")

        out.write("## Rare ClinVar Review Queue\n")
        if review_queue_rows:
            rows = []
            for row in review_queue_rows:
                rows.append([
                    format_variant_label(row),
                    row.get("IMPACT", "") or "n/a",
                    row.get("MAX_AF", "") or "n/a",
                    row.get("CLINVAR_CLNSIG", "") or "n/a",
                    short_disease_list(row.get("CLINVAR_CLNDN", "")),
                ])
            write_table(out, ["Variant", "Impact", "MAX_AF", "ClinVar status", "ClinVar conditions"], rows)
        else:
            out.write("No additional rare ClinVar review items remained after pulling out the reviewed findings.\n\n")

        out.write("## Predicted LoF Bucket\n")
        out.write(f"- Unique predicted LoF/protein-altering hits: `{len(lof_rows)}`\n")
        out.write(f"- Genes shared between the broad ClinVar bucket and the rare review bucket: `{len({row.get('GENE_SYMBOL') or row.get('GENE_ID') for row in overlap_rows})}`\n")
        out.write(f"- Broad ClinVar hits that stayed outside the rare review bucket: `{broad_only_count}`\n")
        if top_lof_genes:
            formatted = ", ".join(f"`{gene}` ({count})" for gene, count in top_lof_genes)
            out.write(f"- Top genes by unique LoF hit count: {formatted}\n")
        out.write("\n")

        out.write("## Supporting Files\n")
        out.write("- `priority_variants.tsv`: broad transcript-level collection of ClinVar P/LP or predicted LoF annotations.\n")
        out.write("- `clinvar_plp_protein_altering.tsv`: broad ClinVar protein-altering transcript rows.\n")
        out.write("- `clinvar_plp_unique_variants.tsv`: deduplicated broad ClinVar variant-gene hits.\n")
        out.write("- `lof_protein_altering.tsv`: broad predicted LoF transcript rows.\n")
        out.write("- `lof_unique_variants.tsv`: deduplicated predicted LoF variant-gene hits.\n")
        out.write("- `rare_clinvar_review_queue.tsv`: rare ClinVar review transcript rows from the stricter triage step.\n")
        out.write("- `rare_clinvar_review_unique_variants.tsv`: deduplicated rare ClinVar review variant-gene hits.\n")


if __name__ == "__main__":
    main()
