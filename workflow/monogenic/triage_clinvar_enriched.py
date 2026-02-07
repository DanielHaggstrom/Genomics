#!/usr/bin/env python3
import argparse
import gzip
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

HIGH_REVIEW_MARKERS = (
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "criteria_provided&_multiple_submitters&_no_conflicts",
)


def info_map(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out


def has_any_consequence(consequence: str, terms: set[str]) -> bool:
    vals = {x.strip() for x in consequence.split("&") if x.strip()}
    return any(x in vals for x in terms)


def is_plp(sig: str) -> bool:
    s = sig.lower()
    return ("pathogenic" in s) or ("likely_pathogenic" in s)


def has_conflict(sig: str, rev: str) -> bool:
    s = (sig + " " + rev).lower()
    return "conflicting_classifications" in s


def high_review(rev: str) -> bool:
    r = rev.lower()
    r_norm = r.replace("&", ",")
    return any((m in r) or (m in r_norm) for m in HIGH_REVIEW_MARKERS)


def parse_float(s: str) -> float | None:
    try:
        return float(s)
    except Exception:
        return None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-vcf", required=True)
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--summary-md", required=True)
    p.add_argument("--max-af", type=float, default=0.01)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_md), exist_ok=True)

    csq_fields = []
    rows = []

    n_total = 0
    c_pass = c_protein = c_plp = c_review = c_noconf = c_rare = 0

    with gzip.open(args.input_vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format: ", 1)[1].rstrip("\n").rstrip('">')
                csq_fields = fmt.split("|")
                continue
            if line.startswith("#"):
                continue
            n_total += 1

            chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample = line.rstrip("\n").split("\t")[:10]
            if flt != "PASS":
                continue
            c_pass += 1

            im = info_map(info)
            csq_blob = im.get("CSQ", "")
            if not csq_blob:
                continue

            fmt_keys = fmt.split(":")
            smp_vals = sample.split(":")
            smp = dict(zip(fmt_keys, smp_vals))

            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                m = dict(zip(csq_fields, parts))

                consequence = m.get("Consequence", "")
                if not has_any_consequence(consequence, PROTEIN_ALTERING):
                    continue
                c_protein += 1

                clnsig = m.get("ClinVar_CLNSIG", "")
                clnrev = m.get("ClinVar_CLNREVSTAT", "")
                if not is_plp(clnsig):
                    continue
                c_plp += 1

                if not high_review(clnrev):
                    continue
                c_review += 1

                if has_conflict(clnsig, clnrev):
                    continue
                c_noconf += 1

                max_af_str = m.get("MAX_AF", "")
                max_af = parse_float(max_af_str)
                if max_af is None or max_af > args.max_af:
                    continue
                c_rare += 1

                rows.append(
                    {
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": vid,
                        "REF": ref,
                        "ALT": alt,
                        "GT": smp.get("GT", "."),
                        "DP": smp.get("DP", "."),
                        "GQ": smp.get("GQ", "."),
                        "GENE_SYMBOL": m.get("SYMBOL", ""),
                        "GENE_ID": m.get("Gene", ""),
                        "CONSEQUENCE": consequence,
                        "IMPACT": m.get("IMPACT", ""),
                        "CANONICAL": m.get("CANONICAL", ""),
                        "HGVSc": m.get("HGVSc", ""),
                        "HGVSp": m.get("HGVSp", ""),
                        "EXISTING_VARIATION": m.get("Existing_variation", ""),
                        "CLINVAR_CLNSIG": clnsig,
                        "CLINVAR_CLNREVSTAT": clnrev,
                        "CLINVAR_CLNDN": m.get("ClinVar_CLNDN", ""),
                        "MAX_AF": max_af_str,
                        "MAX_AF_POPS": m.get("MAX_AF_POPS", ""),
                    }
                )

    # dedup exact rows
    seen = set()
    dedup = []
    for r in rows:
        k = tuple(r.values())
        if k in seen:
            continue
        seen.add(k)
        dedup.append(r)

    header = list(dedup[0].keys()) if dedup else [
        "CHROM", "POS", "ID", "REF", "ALT", "GT", "DP", "GQ", "GENE_SYMBOL", "GENE_ID",
        "CONSEQUENCE", "IMPACT", "CANONICAL", "HGVSc", "HGVSp", "EXISTING_VARIATION",
        "CLINVAR_CLNSIG", "CLINVAR_CLNREVSTAT", "CLINVAR_CLNDN", "MAX_AF", "MAX_AF_POPS"
    ]
    with open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in dedup:
            out.write("\t".join(r[h] for h in header) + "\n")

    top_genes = Counter((r["GENE_SYMBOL"] or r["GENE_ID"] or "UNKNOWN") for r in dedup).most_common(20)
    with open(args.summary_md, "w", encoding="utf-8") as out:
        out.write("# High-Confidence ClinVar Triage\n\n")
        out.write(f"- Total variants in VCF: {n_total}\n")
        out.write(f"- PASS variants: {c_pass}\n")
        out.write(f"- Protein-altering annotations: {c_protein}\n")
        out.write(f"- ClinVar P/LP annotations: {c_plp}\n")
        out.write(f"- High-review (2+ star equivalent) annotations: {c_review}\n")
        out.write(f"- Non-conflicting annotations: {c_noconf}\n")
        out.write(f"- Rare annotations (MAX_AF <= {args.max_af}): {c_rare}\n")
        out.write(f"- Final deduplicated shortlist rows: {len(dedup)}\n\n")
        out.write("## Top genes\n")
        for gene, n in top_genes:
            out.write(f"- {gene}: {n}\n")


if __name__ == "__main__":
    main()
