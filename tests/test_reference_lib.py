from __future__ import annotations

import unittest
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from workflow.core.reference import (
    compare_sequence_catalogs, ensure_bundle_files, extract_reference_hints,
    cmd_prepare_context, parse_sq_entries, prepare_import_artifact, validate_work_roots,
)
from workflow.core.manifest import apply_defaults
from workflow.core.runtime import write_json


class ReferenceLibTests(unittest.TestCase):
    def test_resolved_result_child_cannot_write_into_originals(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            originals = root / "originals"
            originals.mkdir()
            results = root / "results/SYNTHETIC"
            results.mkdir(parents=True)
            (results / "latest").symlink_to(originals, target_is_directory=True)
            manifest = apply_defaults({"sample_id": "SYNTHETIC",
                "inputs": {"cram": str(originals / "sample.cram"), "crai": str(originals / "sample.cram.crai")},
                "runtime": {"linux_cache_root": str(root / "cache")}, "outputs": {"results_dir": str(root / "results")}})
            validated = root / "validated.json"
            write_json(validated, {"manifest_path": str(root / "sample.yaml"), "manifest": manifest})
            with self.assertRaisesRegex(ValueError, "protected input"):
                cmd_prepare_context(argparse.Namespace(validated_manifest=str(validated), output=str(root / "context.json")))
            self.assertEqual(list(originals.iterdir()), [])

    def test_work_roots_reject_protected_tree_and_symlink_overlap(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            originals = tmp / "originals"
            originals.mkdir()
            alias = tmp / "alias"
            alias.symlink_to(originals, target_is_directory=True)
            inputs = [originals / "sample.cram", originals / "sample.cram.crai"]
            for unsafe in (originals, originals / "cache", tmp, alias / "cache"):
                with self.subTest(unsafe=unsafe), self.assertRaisesRegex(ValueError, "protected input"):
                    validate_work_roots(runtime_root=unsafe, results_root=tmp / "results", input_paths=inputs)
            validate_work_roots(runtime_root=tmp / "cache", results_root=tmp / "results", input_paths=inputs)
            with self.assertRaisesRegex(ValueError, "protected input"):
                validate_work_roots(runtime_root=tmp / "cache", results_root=tmp / "archive" / "results",
                                    input_paths=inputs, protected_roots=[tmp / "archive"])

    def test_missing_reference_indexes_are_created_only_on_cache_copy(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            originals = tmp / "originals"
            originals.mkdir()
            fasta = originals / "ref.fa"
            fasta.write_bytes(b">chr1\nACGT\n")
            before = fasta.read_bytes()

            def fake_index(command: list[str], **_: object) -> None:
                target = Path(command[3]) if command[1] == "dict" else Path(command[2] + ".fai")
                self.assertTrue(target.is_relative_to(tmp / "cache"))
                target.write_bytes(b"synthetic index")

            with patch("workflow.core.reference.run_command", side_effect=fake_index) as commands:
                active, fai, seq_dict, missing = ensure_bundle_files({"fasta": str(fasta)}, cache_dir=tmp / "cache")
                count = commands.call_count
                ensure_bundle_files({"fasta": str(fasta)}, cache_dir=tmp / "cache")
                self.assertEqual(commands.call_count, count)
                # A partial index left by an interrupted rebuild must not be reused.
                fai.write_bytes(b"partial index")
                ensure_bundle_files({"fasta": str(fasta)}, cache_dir=tmp / "cache")
                self.assertEqual(commands.call_count, count + 1)
            self.assertFalse(missing)
            self.assertNotEqual(active, fasta)
            self.assertTrue(fai.is_file())
            self.assertTrue(seq_dict.is_file())
            self.assertEqual(fasta.read_bytes(), before)
            self.assertEqual(list(originals.iterdir()), [fasta])

    def test_unindexed_import_is_indexed_on_copy_and_source_index_is_preserved(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            originals = tmp / "originals"
            originals.mkdir()
            vcf = originals / "sample.vcf.gz"
            vcf.write_bytes(b"synthetic input")
            source_index = Path(str(vcf) + ".tbi")
            source_index.touch()

            def fake_index(command: list[str], **_: object) -> None:
                active = Path(command[-1])
                self.assertTrue(active.is_relative_to(tmp / "cache"))
                Path(str(active) + ".tbi").write_bytes(b"synthetic index")

            with patch("workflow.core.reference.run_command", side_effect=fake_index):
                imported = prepare_import_artifact("called_vcf", str(vcf), manifest_dir=tmp,
                                                   staging_root=tmp / "cache" / "imports", linux_cache_root=tmp / "cache",
                                                   current_input_fingerprint_id="synthetic")
            self.assertTrue(imported["staged"])
            self.assertNotEqual(imported["active"], str(vcf))
            self.assertEqual(vcf.read_bytes(), b"synthetic input")
            self.assertEqual(source_index.stat().st_size, 0)

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
