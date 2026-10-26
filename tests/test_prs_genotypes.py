from __future__ import annotations

import gzip
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from workflow.stages.prs_genotypes import augment_scaffold, materialize_scaffold


def vcf(rows: list[str], *, sample: str | None = "SAMPLE") -> StringIO:
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if sample is not None:
        columns += "\tFORMAT\t" + sample
    return StringIO("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" + columns + "\n" + "\n".join(rows) + "\n")


def records(text: str) -> dict[tuple[int, str], tuple[list[str], dict[str, str]]]:
    result = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            row = line.split("\t")
            result[(int(row[1]), row[4])] = (row, dict(zip(row[8].split(":"), row[9].split(":"))))
    return result


class PrsGenotypeTests(unittest.TestCase):
    def test_known_reference_and_alt_dosages_and_missingness(self) -> None:
        calls = vcf([
            "chr1\t10\t.\tA\t.\t.\t.\t.\tGT:DP:RGQ\t0/0:25:55",
            "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:20:60",
            "chr1\t30\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t1|1:30:80",
            "chr1\t40\t.\tT\tC\t.\t.\t.\tGT:DP:GQ\t./.:0:0",
            "chr1\t50\t.\tA\t.\t.\t.\t.\tGT:DP:RGQ\t0/0:2:80",
            "chr1\t60\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:25",
            "chr1\t70\t.\tG\tA\t.\t.\t.\tGT:DP:GQ\t0/1:25:5",
        ])
        targets = vcf([
            f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
            for pos, ref, alt in [(10, "A", "G"), (20, "C", "T"), (30, "G", "A"), (40, "T", "C"), (50, "A", "G"), (60, "C", "T"), (70, "G", "A"), (80, "A", "C")]
        ], sample=None)
        output = StringIO()
        summary = materialize_scaffold(calls, targets, output, sample_id="RENAMED")
        observed = records(output.getvalue())
        self.assertEqual(summary["callable_sites"], 3)
        self.assertEqual(summary["missing_sites"], 5)
        self.assertEqual([observed[key][1]["GT"] for key in [(10, "G"), (20, "T"), (30, "A")]], ["0/0", "0/1", "1|1"])
        self.assertEqual(observed[(10, "G")][1]["RGQ"], "55")
        self.assertIn("\tRENAMED\n", output.getvalue())
        self.assertEqual([observed[key][1]["FT"] for key in [(40, "C"), (50, "G"), (60, "T"), (70, "A"), (80, "C")]], ["no_call", "low_depth", "missing_quality", "low_genotype_quality", "no_record"])
        # Independent published-score arithmetic for this tiny fixture:
        # A(ref at 10)*-0.5 + T(alt at 20)*2 + A(alt at 30)*3 = -1+2+6=7.
        score = 0.0
        for key, effect, weight in [((10, "G"), "A", -0.5), ((20, "T"), "T", 2), ((30, "A"), "A", 3)]:
            row, values = observed[key]
            bases = [row[3], row[4]]
            genotype = values["GT"].replace("|", "/").split("/")
            score += sum(bases[int(index)] == effect for index in genotype) * weight
        self.assertEqual(score, 7.0)
        self.assertTrue(all(observed[key][1]["GT"] == "./." for key in [(40, "C"), (50, "G"), (60, "T"), (70, "A"), (80, "C")]))

    def test_indels_multiallelic_calls_and_ploidy_are_conservative(self) -> None:
        calls = vcf([
            "chr1\t10\t.\tA\t.\t.\t.\t.\tGT:DP:RGQ\t0/0:25:60",
            "chr1\t20\t.\tC\tCT\t50\tPASS\t.\tGT:DP:GQ\t0/1:25:60",
            "chr1\t30\t.\tA\tG,C\t50\tPASS\t.\tGT:DP:GQ\t1/2:25:60",
            "chr1\t40\t.\tT\tC\t50\tPASS\t.\tGT:DP:GQ\t1:25:60",
            "chr1\t50\t.\tG\tA\t50\tLowQual\t.\tGT:DP:GQ\t0/1:25:60",
            "chr1\t60\t.\tA\tC\t50\tPASS\t.\tGT:DP:GQ\t0/1:25:60",
        ])
        targets = vcf([
            f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
            for pos, ref, alt in [(10, "A", "AT"), (20, "C", "CT"), (30, "A", "G,C"), (40, "T", "C"), (50, "G", "A"), (60, "G", "C")]
        ], sample=None)
        output = StringIO()
        summary = materialize_scaffold(calls, targets, output)
        observed = records(output.getvalue())
        self.assertEqual(summary["total_sites"], 7)
        self.assertEqual(summary["callable_sites"], 1)
        self.assertEqual(observed[(20, "CT")][1]["GT"], "0/1")
        self.assertEqual(observed[(10, "AT")][1]["FT"], "unresolved_indel")
        self.assertEqual(observed[(30, "G")][1]["FT"], "other_allele")
        self.assertEqual(observed[(30, "C")][1]["FT"], "other_allele")
        self.assertEqual(observed[(40, "C")][1]["FT"], "unsupported_ploidy")
        self.assertEqual(observed[(50, "A")][1]["FT"], "filtered_call")
        self.assertEqual(observed[(60, "C")][1]["FT"], "reference_mismatch")

    def test_rejects_unfinalized_and_unsorted_inputs(self) -> None:
        target_text = "chr1\t10\t.\tA\tG\t.\t.\t."
        for rows in [
            ["chr1\t10\t.\tA\t<NON_REF>\t.\t.\t.\tGT:DP:GQ\t0/0:20:50"],
            ["chr1\t10\t.\tA\t.\t.\t.\tEND=20\tGT:DP:GQ\t0/0:20:50"],
            ["chr1\t20\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/0:20:50", "chr1\t10\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/0:20:50"],
        ]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                materialize_scaffold(vcf(rows), vcf([target_text], sample=None), StringIO())

    def test_imputation_augments_direct_calls_and_retains_direct_only_sites(self) -> None:
        direct = vcf([
            "chr1\t10\t.\tA\tG\t.\t.\t.\tGT:DP:GQ:FT\t0/0:25:60:PASS",
            "chr1\t20\t.\tC\tT\t.\t.\t.\tGT:DP:GQ:FT\t0/1:25:60:PASS",
            "chr1\t30\t.\tG\tA\t.\t.\t.\tGT:DP:GQ:FT\t./.:2:5:low_depth",
            "chr1\t50\t.\tT\tC\t.\t.\t.\tGT:DP:GQ:FT\t./.:0:.:no_record",
        ])
        imputed = vcf([
            "chr1\t10\t.\tA\tG\t.\t.\t.\tGT:DS\t1/1:1.8",
            "chr1\t30\t.\tG\tA\t.\t.\t.\tGT:DS\t0|1:0.8",
            "chr1\t40\t.\tA\tT\t.\t.\t.\tGT:DS\t0/0:0.02",
            "chr1\t50\t.\tT\tC\t.\t.\t.\tGT:DS\t0/1:nan",
        ])
        output = StringIO()
        summary = augment_scaffold(direct, imputed, output)
        observed = records(output.getvalue())
        self.assertEqual(summary["sources"], {"direct": 2, "imputed": 2, "missing": 1})
        self.assertEqual(observed[(10, "G")][1]["GT"], "0/0")
        self.assertEqual(observed[(10, "G")][1]["DS"], "0")
        self.assertEqual(observed[(20, "T")][1]["DS"], "1")
        self.assertEqual(observed[(30, "A")][1]["DS"], "0.8")
        self.assertEqual(observed[(40, "T")][1]["DS"], "0.02")
        self.assertEqual(observed[(50, "C")][1]["GT"], "./.")
        self.assertEqual(observed[(50, "C")][1]["DS"], ".")

    def test_merge_rejects_different_samples(self) -> None:
        with self.assertRaisesRegex(ValueError, "same sample"):
            augment_scaffold(vcf([], sample="ONE"), vcf([], sample="TWO"), StringIO())

    def test_imputation_cannot_invent_a_different_allele_at_a_called_position(self) -> None:
        direct = vcf(["chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:FT\t0/1:PASS"])
        imputed = vcf(["chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:DS\t0/1:0.9"])
        output = StringIO()
        augment_scaffold(direct, imputed, output)
        observed = records(output.getvalue())
        self.assertEqual(observed[(10, "G")][1]["GT"], "0/1")
        self.assertEqual(observed[(10, "C")][1]["GT"], "./.")

    def test_cli_reads_gzip_fixture_and_writes_qc_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = root / "calls.vcf.gz"
            targets = root / "targets.vcf"
            output = root / "output.vcf"
            summary = root / "summary.json"
            with gzip.open(calls, "wt") as handle:
                handle.write(vcf(["chr1\t10\t.\tA\t.\t.\t.\t.\tGT:DP:RGQ\t0/0:25:60"]).getvalue())
            targets.write_text(vcf(["chr1\t10\t.\tA\tG\t.\t.\t."], sample=None).getvalue())
            subprocess.run([sys.executable, "-m", "workflow.stages.prs_genotypes", "--calls", str(calls), "--targets", str(targets), "--output", str(output), "--summary", str(summary)], cwd=Path(__file__).resolve().parents[1], check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(summary.read_text())["callable_sites"], 1)
            self.assertEqual(records(output.read_text())[(10, "G")][1]["GT"], "0/0")
            before = calls.read_bytes()
            alias = root / "original-alias.vcf.gz"
            alias.hardlink_to(calls)
            rejected = subprocess.run([sys.executable, "-m", "workflow.stages.prs_genotypes", "--calls", str(calls), "--targets", str(targets), "--output", str(alias), "--summary", str(summary)], cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("Refusing to overwrite", rejected.stderr)
            self.assertEqual(calls.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
