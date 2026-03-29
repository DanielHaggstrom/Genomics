#!/usr/bin/env python3
import argparse
import gzip
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

HIGH_REVIEW_MARKERS = (
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "criteria_provided&_multiple_submitters&_no_conflicts",
)

SINGLE_REVIEW_MARKER = "criteria_provided&_single_submitter"

LOF_TERMS = {
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "transcript_ablation",
}


def info_map(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out


def carried_alt_alleles(gt: str, alt: str) -> set[str]:
    alts = alt.split(",")
    carried = set()
    for token in gt.replace("|", "/").split("/"):
        if token in {"", ".", "0"}:
            continue
        try:
            idx = int(token) - 1
        except ValueError:
            continue
        if 0 <= idx < len(alts):
            carried.add(alts[idx])
    return carried


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


def single_submitter_review(rev: str) -> bool:
    return SINGLE_REVIEW_MARKER in rev.lower()


def review_tier(rev: str) -> str:
    if not rev:
        return "none"
    if high_review(rev):
        return "high"
    if single_submitter_review(rev):
        return "single_submitter"
    return "other"


def has_lof(consequence: str) -> bool:
    vals = {x.strip() for x in consequence.split("&") if x.strip()}
    return any(x in vals for x in LOF_TERMS)


def parse_float(s: str) -> float | None:
    try:
        return float(s)
    except Exception:
        return None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-vcf", required=True)
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--mode", choices=["strict", "moderate", "broad"], default="moderate")
    p.add_argument("--max-af", type=float, default=None)
    args = p.parse_args()

    if args.max_af is None:
        if args.mode == "strict":
            args.max_af = 0.01
        elif args.mode == "moderate":
            args.max_af = 0.05
        else:
            args.max_af = 0.10

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    csq_fields = []
    rows = []

    n_total = 0
    c_pass = c_protein = c_plp = c_review = c_noconf = c_rare = 0
    c_selected = 0

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
            carried_alts = carried_alt_alleles(smp.get("GT", "."), alt) if "," in alt else set()

            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                m = dict(zip(csq_fields, parts))

                if carried_alts:
                    allele = m.get("Allele", "")
                    if allele and allele not in carried_alts:
                        continue

                consequence = m.get("Consequence", "")
                if not has_any_consequence(consequence, PROTEIN_ALTERING):
                    continue
                c_protein += 1

                clnsig = m.get("ClinVar_CLNSIG", "")
                clnrev = m.get("ClinVar_CLNREVSTAT", "")
                if not is_plp(clnsig):
                    continue
                c_plp += 1

                rev_tier = review_tier(clnrev)
                if rev_tier == "high":
                    c_review += 1

                conflict = has_conflict(clnsig, clnrev)
                if not conflict:
                    c_noconf += 1

                consequence = m.get("Consequence", "")
                impact = m.get("IMPACT", "")
                canonical = m.get("CANONICAL", "")
                max_af_str = m.get("MAX_AF", "")
                max_af = parse_float(max_af_str)
                if max_af is not None and max_af <= args.max_af:
                    c_rare += 1

                # Selection logic by mode (all modes keep: PASS, protein-altering, ClinVar P/LP, non-conflicting)
                select = False
                reason = ""
                if args.mode == "strict":
                    if (not conflict) and rev_tier == "high" and max_af is not None and max_af <= args.max_af:
                        select = True
                        reason = "high_review_and_rare"
                elif args.mode == "moderate":
                    # Less restrictive but still selective:
                    # keep rare P/LP protein-altering variants, retain conflict flag for downstream review.
                    if max_af is not None and max_af <= args.max_af and (rev_tier != "none" or impact == "HIGH" or has_lof(consequence)):
                        select = True
                        if conflict:
                            reason = "rare_plp_conflicting"
                        elif rev_tier in {"high", "single_submitter"}:
                            reason = "rare_plp_reviewed"
                        else:
                            reason = "rare_plp_other_review"
                    elif max_af is None and (impact == "HIGH" or has_lof(consequence)) and rev_tier != "none":
                        select = True
                        if conflict:
                            reason = "high_impact_conflicting_missing_af"
                        else:
                            reason = "high_impact_reviewed_missing_af"
                else:
                    # Broad but still selective: allow any non-conflicting review tier,
                    # require either low-frequency evidence or severe consequence.
                    if (not conflict) and max_af is not None and max_af <= args.max_af:
                        select = True
                        reason = "nonconflicting_plp_low_frequency"
                    elif (not conflict) and max_af is None and (impact == "HIGH" or has_lof(consequence) or canonical == "YES"):
                        select = True
                        reason = "nonconflicting_plp_severe_missing_af"
                    elif conflict and max_af is not None and max_af <= 0.01 and (impact == "HIGH" or has_lof(consequence)):
                        select = True
                        reason = "conflicting_but_rare_high_impact"

                if not select:
                    continue
                c_selected += 1

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
                        "IMPACT": impact,
                        "CANONICAL": canonical,
                        "HGVSc": m.get("HGVSc", ""),
                        "HGVSp": m.get("HGVSp", ""),
                        "EXISTING_VARIATION": m.get("Existing_variation", ""),
                        "CLINVAR_CLNSIG": clnsig,
                        "CLINVAR_CLNREVSTAT": clnrev,
                        "CLINVAR_CLNDN": m.get("ClinVar_CLNDN", ""),
                        "MAX_AF": max_af_str,
                        "MAX_AF_POPS": m.get("MAX_AF_POPS", ""),
                        "CLINVAR_CONFLICT_FLAG": "YES" if conflict else "NO",
                        "REVIEW_TIER": rev_tier,
                        "SELECTION_MODE": args.mode,
                        "SELECTION_REASON": reason,
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
        "CLINVAR_CLNSIG", "CLINVAR_CLNREVSTAT", "CLINVAR_CLNDN", "MAX_AF", "MAX_AF_POPS",
        "CLINVAR_CONFLICT_FLAG",
        "REVIEW_TIER", "SELECTION_MODE", "SELECTION_REASON"
    ]
    with open(args.output_tsv, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in dedup:
            out.write("\t".join(r[h] for h in header) + "\n")


if __name__ == "__main__":
    main()
