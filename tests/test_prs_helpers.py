from __future__ import annotations

import unittest

from workflow.stages.prs import (
    parse_base_tokens,
    parse_overlap_diagnostics,
    plink_sex_code,
    resolve_score_row_alleles,
    safe_sampleset_name,
)


class PrsHelperTests(unittest.TestCase):
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
                "hm_inferOtherAllele": "G/T",
            },
            "C",
        )
        self.assertIsNone(resolved)
        self.assertEqual(reason, "reference_not_in_score_alleles")


if __name__ == "__main__":
    unittest.main()
