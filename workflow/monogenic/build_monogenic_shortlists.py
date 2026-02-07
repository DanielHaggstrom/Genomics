#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter


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

    write_tsv(os.path.join(args.outdir, "clinvar_plp_protein_altering.tsv"), clinvar_plp, header)
    write_tsv(os.path.join(args.outdir, "lof_protein_altering.tsv"), lof, header)

    gene_counts = Counter((x.get("GENE_SYMBOL") or x.get("GENE_ID") or "UNKNOWN") for x in clinvar_plp)
    top_genes = gene_counts.most_common(25)

    summary_path = os.path.join(args.outdir, "shortlist_summary.md")
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("# Monogenic Shortlist Summary\n\n")
        out.write(f"- Input rows: {len(rows)}\n")
        out.write(f"- PASS rows: {len(pass_rows)}\n")
        out.write(f"- ClinVar P/LP + protein-altering: {len(clinvar_plp)}\n")
        out.write(f"- LoF + protein-altering: {len(lof)}\n\n")
        out.write("## Top Genes (ClinVar P/LP + protein-altering)\n")
        for gene, n in top_genes:
            out.write(f"- {gene}: {n}\n")


if __name__ == "__main__":
    main()
