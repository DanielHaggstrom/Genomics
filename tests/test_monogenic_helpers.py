from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from workflow.stages.monogenic import (
    extract_priority_variants,
    iter_annotated_variants,
    triage_clinvar,
)


CSQ_FIELDS = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "HGVSc",
    "HGVSp",
    "Existing_variation",
    "CANONICAL",
    "CLIN_SIG",
    "MAX_AF",
    "MAX_AF_POPS",
    "ClinVar_CLNSIG",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNDN",
]

TEST_VCF = """##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {fields}">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100\t.\tA\tG,T\t50\tPASS\tCSQ=G|missense_variant|MODERATE|GENE1|ENSG1|c.1A>G|p.Lys1Arg|rs1|YES|pathogenic|0.001|gnomAD|Pathogenic|practice_guideline|Disease1,T|missense_variant|MODERATE|GENE2|ENSG2|c.1A>T|p.Lys1Asn|rs2||benign|0.001|gnomAD|Benign|practice_guideline|Disease2\tGT:DP:GQ:AD\t1/1:30:99:10,20,0
chr2\t200\trs3\tC\tT\t60\tPASS\tCSQ=T|stop_gained|HIGH|GENE3|ENSG3|c.2C>T|p.Arg2Ter|rs3|YES|likely_pathogenic||.|Likely_pathogenic|criteria_provided&_single_submitter|Disease3\tGT:DP:GQ:AD\t0/1:25:80:12,13
""".format(fields="|".join(CSQ_FIELDS))


class MonogenicHelperTests(unittest.TestCase):
    def write_test_vcf(self, path: Path) -> None:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(TEST_VCF)

    def test_iter_annotated_variants_filters_multiallelic_annotations_to_carried_alt(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = Path(tmpdir) / "annotated.vcf.gz"
            self.write_test_vcf(vcf_path)

            records = list(iter_annotated_variants(vcf_path))

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["ALT"], "G,T")
        self.assertEqual(len(records[0]["annotations"]), 1)
        self.assertEqual(records[0]["annotations"][0]["Allele"], "G")
        self.assertEqual(records[0]["annotations"][0]["SYMBOL"], "GENE1")

    def test_extract_priority_variants_uses_shared_annotation_iterator(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf_path = tmp / "annotated.vcf.gz"
            output_tsv = tmp / "priority.tsv"
            summary_txt = tmp / "priority.summary.txt"
            self.write_test_vcf(vcf_path)

            summary = extract_priority_variants(vcf_path, output_tsv, summary_txt)

            with open(output_tsv, "r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary["total_variant_records_seen"], "2")
        self.assertEqual(summary["priority_rows_after_dedup"], "2")
        self.assertEqual(summary["count_clinvar_pathogenic_or_likely_pathogenic"], "2")
        self.assertEqual(summary["count_predicted_loss_of_function"], "1")
        self.assertEqual([row["GENE_SYMBOL"] for row in rows], ["GENE1", "GENE3"])
        self.assertEqual(rows[0]["CLIN_SIG"], "pathogenic")
        self.assertEqual(rows[1]["REASONS"], "clinvar_pathogenic_or_likely_pathogenic,predicted_loss_of_function")

    def test_triage_clinvar_uses_shared_annotation_iterator(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vcf_path = tmp / "annotated.vcf.gz"
            output_tsv = tmp / "review.tsv"
            self.write_test_vcf(vcf_path)

            triage_clinvar(vcf_path, output_tsv, mode="moderate", max_af=0.01)

            with open(output_tsv, "r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(rows), 2)
        self.assertEqual([row["GENE_SYMBOL"] for row in rows], ["GENE1", "GENE3"])
        self.assertEqual(rows[0]["REVIEW_TIER"], "high")
        self.assertEqual(rows[0]["SELECTION_REASON"], "rare_plp_reviewed")
        self.assertEqual(rows[1]["REVIEW_TIER"], "single_submitter")
        self.assertEqual(rows[1]["SELECTION_REASON"], "high_impact_missing_af")


if __name__ == "__main__":
    unittest.main()
