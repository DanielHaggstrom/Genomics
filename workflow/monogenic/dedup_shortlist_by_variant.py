#!/usr/bin/env python3
import argparse
import csv


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-tsv", required=True)
    p.add_argument("--output-tsv", required=True)
    args = p.parse_args()

    with open(args.input_tsv, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        rows = list(r)

    grouped = {}
    for row in rows:
        key = (row.get("CHROM", ""), row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))
        grouped.setdefault(key, []).append(row)

    picked = []
    for key in sorted(grouped.keys()):
        bucket = grouped[key]
        canonical = [x for x in bucket if x.get("CANONICAL") == "YES"]
        picked.append(canonical[0] if canonical else bucket[0])

    with open(args.output_tsv, "w", encoding="utf-8", newline="") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(picked)


if __name__ == "__main__":
    main()
