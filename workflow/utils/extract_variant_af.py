#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-vcf", required=True)
    p.add_argument("--targets", required=True, help="Comma list: chr:pos:rsid")
    args = p.parse_args()

    targets = set()
    for t in args.targets.split(","):
        chrom, pos, rsid = t.split(":")
        targets.add((chrom, pos, rsid))

    csq_fields = []
    found = set()
    with gzip.open(args.input_vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format: ", 1)[1].rstrip("\n").rstrip('">')
                csq_fields = fmt.split("|")
                continue
            if line.startswith("#"):
                continue

            chrom, pos, rsid, ref, alt, qual, flt, info, *_ = line.rstrip("\n").split("\t")
            key = (chrom, pos, rsid)
            if key not in targets:
                continue
            found.add(key)

            csq_blob = ""
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq_blob = item[4:]
                    break

            afs = defaultdict(set)
            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                m = dict(zip(csq_fields, parts))
                for fld in ["AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "MAX_AF_POPS"]:
                    v = m.get(fld, "")
                    if v:
                        afs[fld].add(v)

            print(f"{chrom}:{pos}:{rsid} {ref}>{alt}")
            for fld in ["AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "MAX_AF_POPS"]:
                vals = sorted(afs[fld])
                print(f"  {fld}: {', '.join(vals) if vals else '.'}")

    missing = targets - found
    if missing:
        print("MISSING_TARGETS:")
        for m in sorted(missing):
            print(f"  {m[0]}:{m[1]}:{m[2]}")


if __name__ == "__main__":
    main()
