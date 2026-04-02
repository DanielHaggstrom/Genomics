from __future__ import annotations

import unittest
from pathlib import Path
import tempfile

from workflow.stages.prs import (
    canonical_score_id,
    complement_base_tokens,
    parse_base_tokens,
    parse_overlap_diagnostics,
    plink_sex_code,
    resolve_cached_scorefile,
    resolve_score_row_alleles,
    resolve_score_row_alleles_for_site,
    safe_sampleset_name,
)


class FakeFasta:
    def __init__(self, sequence_by_length: dict[tuple[str, int, int], str]) -> None:
        self.sequence_by_length = sequence_by_length

    def get_sequence(self, contig: str, pos: int, sequence_length: int) -> str:
        return self.sequence_by_length[(contig, pos, sequence_length)]


class PrsHelperTests(unittest.TestCase):
    def test_canonical_score_id_strips_harmonized_suffix(self) -> None:
        self.assertEqual(canonical_score_id("PGS999001_hmPOS_GRCh38"), "PGS999001")
        self.assertEqual(canonical_score_id("PGS999001_hmPOS_GRCh37"), "PGS999001")

    def test_resolve_cached_scorefile_prefers_exact_harmonized_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            exact = root / "PGS999001_hmPOS_GRCh38.txt.gz"
            other = root / "PGS999001_other.txt.gz"
            exact.write_text("x", encoding="utf-8")
            other.write_text("y", encoding="utf-8")
            self.assertEqual(resolve_cached_scorefile(root, "PGS999001"), exact)

    def test_safe_sampleset_name_normalizes_reserved_and_symbols(self) -> None:
        self.assertEqual(safe_sampleset_name("DEMO_SAMPLE_001"), "DEMO-SAMPLE-001")
        self.assertEqual(safe_sampleset_name("bad_name with spaces"), "bad-name-with-spaces")
        self.assertEqual(safe_sampleset_name("reference"), "sample")

    def test_parse_overlap_diagnostics(self) -> None:
        log = """
        Score PGS999001_hmPOS_GRCh38 fails minimum matching threshold (37.95% variants match)
        Score PGS999002_hmPOS_GRCh38 fails minimum matching threshold (56.25% variants match)
        Score PGS999001_hmPOS_GRCh38 fails minimum matching threshold (37.95% variants match)
        """
        parsed = parse_overlap_diagnostics(log)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0]["score"], "PGS999001_hmPOS_GRCh38")
        self.assertEqual(parsed[0]["match_percent"], "37.95")

    def test_plink_sex_code(self) -> None:
        self.assertEqual(plink_sex_code("male"), "M")
        self.assertEqual(plink_sex_code("female"), "F")
        self.assertIsNone(plink_sex_code(None))

    def test_parse_base_tokens_accepts_single_base_tokens_only(self) -> None:
        self.assertEqual(parse_base_tokens("G/T"), {"G", "T"})
        self.assertEqual(parse_base_tokens(" C / G / T "), {"C", "G", "T"})
        self.assertEqual(parse_base_tokens("NA"), set())
        self.assertEqual(parse_base_tokens("<DEL>"), set())

    def test_complement_base_tokens(self) -> None:
        self.assertEqual(complement_base_tokens({"A", "C"}), {"T", "G"})

    def test_resolve_score_row_alleles_prefers_explicit_other_allele(self) -> None:
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": "A",
                "other_allele": "G",
                "hm_inferOtherAllele": "C/T",
            },
            "G",
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["source"], "explicit_other_allele")
        self.assertEqual(resolved["alts"], ["A"])
        self.assertFalse(resolved["multiallelic"])
        self.assertFalse(resolved["strand_flipped"])

    def test_resolve_score_row_alleles_uses_hm_infer_other_allele(self) -> None:
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": "A",
                "hm_inferOtherAllele": "G/T",
            },
            "G",
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["source"], "hm_inferOtherAllele")
        self.assertEqual(resolved["alts"], ["A", "T"])
        self.assertTrue(resolved["multiallelic"])
        self.assertFalse(resolved["strand_flipped"])

    def test_resolve_score_row_alleles_rescues_complemented_orientation(self) -> None:
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": "A",
                "other_allele": "G",
            },
            "C",
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["effect_allele"], "T")
        self.assertEqual(resolved["other_alleles"], ["C"])
        self.assertEqual(resolved["alts"], ["T"])
        self.assertTrue(resolved["strand_flipped"])

    def test_resolve_score_row_alleles_reports_missing_other_allele(self) -> None:
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": "A",
            },
            "A",
        )
        self.assertIsNone(resolved)
        self.assertEqual(reason, "missing_other_allele")

    def test_resolve_score_row_alleles_reports_reference_mismatch(self) -> None:
        resolved, reason = resolve_score_row_alleles(
            {
                "effect_allele": "A",
                "hm_inferOtherAllele": "T",
            },
            "C",
        )
        self.assertIsNone(resolved)
        self.assertEqual(reason, "reference_not_in_score_alleles")

    def test_resolve_score_row_alleles_for_site_supports_indel_style_rows(self) -> None:
        fasta = FakeFasta(
            {
                ("chr10", 1560700, 1): "C",
                ("chr10", 1560700, 2): "CT",
            }
        )
        resolved, reason = resolve_score_row_alleles_for_site(
            {
                "effect_allele": "CT",
                "other_allele": "C",
            },
            fasta=fasta,
            contig="chr10",
            pos=1560700,
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["ref_allele"], "CT")
        self.assertEqual(resolved["alts"], ["C"])

    def test_resolve_score_row_alleles_for_site_supports_insertion_style_rows(self) -> None:
        fasta = FakeFasta(
            {
                ("chr10", 1560700, 1): "C",
                ("chr10", 1560700, 2): "CA",
            }
        )
        resolved, reason = resolve_score_row_alleles_for_site(
            {
                "effect_allele": "CT",
                "other_allele": "C",
            },
            fasta=fasta,
            contig="chr10",
            pos=1560700,
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["ref_allele"], "C")
        self.assertEqual(resolved["alts"], ["CT"])

    def test_resolve_score_row_alleles_for_site_supports_effect_only_rows(self) -> None:
        fasta = FakeFasta(
            {
                ("chr5", 12345, 1): "G",
            }
        )
        resolved, reason = resolve_score_row_alleles_for_site(
            {
                "effect_allele": "A",
            },
            fasta=fasta,
            contig="chr5",
            pos=12345,
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["source"], "effect_allele_only")
        self.assertEqual(resolved["ref_allele"], "G")
        self.assertEqual(resolved["alts"], ["A"])

    def test_resolve_score_row_alleles_for_site_supports_effect_only_reference_calls(self) -> None:
        fasta = FakeFasta(
            {
                ("chr5", 12345, 1): "A",
            }
        )
        resolved, reason = resolve_score_row_alleles_for_site(
            {
                "effect_allele": "A",
            },
            fasta=fasta,
            contig="chr5",
            pos=12345,
        )
        self.assertIsNone(reason)
        self.assertIsNotNone(resolved)
        assert resolved is not None
        self.assertEqual(resolved["source"], "effect_allele_only")
        self.assertEqual(resolved["ref_allele"], "A")
        self.assertEqual(resolved["alts"], [])


if __name__ == "__main__":
    unittest.main()
