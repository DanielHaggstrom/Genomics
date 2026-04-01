from __future__ import annotations

import unittest

from workflow.stages.prs import parse_overlap_diagnostics, plink_sex_code, safe_sampleset_name


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


if __name__ == "__main__":
    unittest.main()
