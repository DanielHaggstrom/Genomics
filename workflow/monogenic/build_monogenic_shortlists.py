#!/usr/bin/env python3
import argparse
import csv
import os


PROTEIN_ALTERING = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
}

LOF_TERMS = {
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "transcript_ablation",
}


def has_any_term(consequence: str, terms: set[str]) -> bool:
    cset = set()
    for block in consequence.split("&"):
        block = block.strip()
        if block:
            cset.add(block)
    return any(t in cset for t in terms)


def write_tsv(path: str, rows: list[dict], header: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def dedup_variant_gene_rows(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, ...], list[dict]] = {}
    for row in rows:
        key = (
            row.get("CHROM", ""),
            row.get("POS", ""),
            row.get("REF", ""),
            row.get("ALT", ""),
            row.get("GENE_SYMBOL", ""),
            row.get("GENE_ID", ""),
        )
        grouped.setdefault(key, []).append(row)

    picked = []
    for key in sorted(grouped.keys()):
        bucket = grouped[key]
        canonical = [x for x in bucket if x.get("CANONICAL") == "YES"]
        picked.append(canonical[0] if canonical else bucket[0])
    return picked


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-tsv", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.input_tsv, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        header = r.fieldnames or []

    pass_rows = [x for x in rows if x.get("FILTER") == "PASS"]
    clinvar_plp = [
        x for x in pass_rows
        if "clinvar_pathogenic_or_likely_pathogenic" in (x.get("REASONS") or "")
        and has_any_term(x.get("CONSEQUENCE", ""), PROTEIN_ALTERING)
    ]
    lof = [
        x for x in pass_rows
        if has_any_term(x.get("CONSEQUENCE", ""), LOF_TERMS)
        and has_any_term(x.get("CONSEQUENCE", ""), PROTEIN_ALTERING)
    ]
    clinvar_unique = dedup_variant_gene_rows(clinvar_plp)
    lof_unique = dedup_variant_gene_rows(lof)

    write_tsv(os.path.join(args.outdir, "clinvar_plp_protein_altering.tsv"), clinvar_plp, header)
    write_tsv(os.path.join(args.outdir, "clinvar_plp_unique_variants.tsv"), clinvar_unique, header)
    write_tsv(os.path.join(args.outdir, "lof_protein_altering.tsv"), lof, header)
    write_tsv(os.path.join(args.outdir, "lof_unique_variants.tsv"), lof_unique, header)


if __name__ == "__main__":
    main()
