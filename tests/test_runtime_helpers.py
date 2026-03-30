from __future__ import annotations

import unittest

from workflow.core.runtime import plan_contig_shards


class RuntimeHelperTests(unittest.TestCase):
    def test_plan_contig_shards_groups_non_primary_contigs_into_alt(self) -> None:
        shards = plan_contig_shards(["chr1", "chr2", "chrUn_KI270442v1", "chrEBV"], "chr")
        self.assertEqual(shards[0]["name"], "chr1")
        self.assertEqual(shards[1]["name"], "chr2")
        self.assertEqual(shards[2]["name"], "alt")
        self.assertEqual(shards[2]["intervals"], ["chrUn_KI270442v1", "chrEBV"])


if __name__ == "__main__":
    unittest.main()
