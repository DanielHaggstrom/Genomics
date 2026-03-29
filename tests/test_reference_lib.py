from __future__ import annotations

import unittest

from workflow.core.reference_lib import compare_sequence_catalogs, extract_reference_hints, parse_sq_entries


class ReferenceLibTests(unittest.TestCase):
    def test_compare_sequence_catalogs_requires_all_cram_contigs_to_match(self) -> None:
        cram_entries = parse_sq_entries(
            "@SQ\tSN:chr1\tLN:10\tM5:aaa\n"
            "@SQ\tSN:chr2\tLN:20\tM5:bbb\n"
        )
        candidate_entries = parse_sq_entries(
            "@SQ\tSN:chr1\tLN:10\tM5:aaa\n"
            "@SQ\tSN:chr2\tLN:20\tM5:ccc\n"
            "@SQ\tSN:chr3\tLN:30\tM5:ddd\n"
        )

        comparison = compare_sequence_catalogs(cram_entries, candidate_entries)
        self.assertFalse(comparison["exact_for_cram"])
        self.assertFalse(comparison["exact_dictionary_match"])
        self.assertEqual(comparison["md5_mismatches"][0]["contig"], "chr2")
        self.assertEqual(comparison["extra_contigs"], ["chr3"])

    def test_extract_reference_hints_finds_dragen_paths(self) -> None:
        pg_entries = [
            {
                "command_line": (
                    "dragen --ht-reference=/refs/hg38.fa --ref-dir /scratch/hg38-alt_masked/DRAGEN/11 "
                    "--variant-annotation-assembly GRCh38"
                )
            }
        ]

        hints = extract_reference_hints(pg_entries)
        self.assertEqual(hints["aligner"], "DRAGEN")
        self.assertEqual(hints["ht_reference"], "/refs/hg38.fa")
        self.assertEqual(hints["ref_dir"], "/scratch/hg38-alt_masked/DRAGEN/11")
        self.assertEqual(hints["variant_annotation_assembly"], "GRCh38")


if __name__ == "__main__":
    unittest.main()
