from __future__ import annotations

import gzip
import subprocess
import unittest
from pathlib import Path
import tempfile

from workflow.stages.prs import (
    base_contig_name,
    beagle_failed_due_to_no_common_markers,
    beagle_log_has_no_common_markers,
    build_pgsc_calc_resource_config,
    build_pgsc_calc_samplesheet,
    canonical_score_id,
    complement_base_tokens,
    inject_missing_contig_headers,
    iter_reference_db_sites,
    target_has_imputable_markers,
    load_imputation_panel,
    parse_pvar_header,
    parse_base_tokens,
    parse_overlap_diagnostics,
    plink_sex_code,
    resolve_panel_resources_for_contig,
    resolve_cached_scorefile,
    resolve_score_row_alleles,
    resolve_score_row_alleles_for_site,
    safe_sampleset_name,
    strip_vcf_suffix,
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
            with gzip.open(exact, "wt", encoding="utf-8") as handle:
                handle.write("x\n")
            with gzip.open(other, "wt", encoding="utf-8") as handle:
                handle.write("y\n")
            self.assertEqual(resolve_cached_scorefile(root, "PGS999001"), exact)

    def test_safe_sampleset_name_normalizes_reserved_and_symbols(self) -> None:
        self.assertEqual(safe_sampleset_name("DEMO_SAMPLE_001"), "DEMO-SAMPLE-001")
        self.assertEqual(safe_sampleset_name("bad_name with spaces"), "bad-name-with-spaces")
        self.assertEqual(safe_sampleset_name("reference"), "sample")

    def test_strip_vcf_suffix_handles_bgz_and_vcf(self) -> None:
        self.assertEqual(strip_vcf_suffix(Path("/tmp/input.vcf.gz")), Path("/tmp/input"))
        self.assertEqual(strip_vcf_suffix(Path("/tmp/input.vcf")), Path("/tmp/input"))

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

    def test_build_pgsc_calc_samplesheet_for_vcf_with_dosage(self) -> None:
        text = build_pgsc_calc_samplesheet(
            sampleset="TEST",
            scoring_input_info={"path_prefix": "/cache/input", "format": "vcf"},
            dosage_field="DS",
        )
        self.assertEqual(
            text,
            "sampleset,path_prefix,chrom,format,vcf_genotype_field\n"
            "TEST,/cache/input,,vcf,DS\n",
        )

    def test_build_pgsc_calc_samplesheet_for_pfile_omits_dosage_field(self) -> None:
        text = build_pgsc_calc_samplesheet(
            sampleset="TEST",
            scoring_input_info={"path_prefix": "/cache/input", "format": "pfile"},
            dosage_field="DS",
        )
        self.assertEqual(
            text,
            "sampleset,path_prefix,chrom,format\n"
            "TEST,/cache/input,,pfile\n",
        )

    def test_build_pgsc_calc_resource_config_targets_match_and_score_processes(self) -> None:
        text = build_pgsc_calc_resource_config(max_memory="14.GB", match_cpus=1, score_cpus=1)
        self.assertIn("MATCH:MATCH_VARIANTS", text)
        self.assertIn("APPLY_SCORE:PLINK2_SCORE", text)
        self.assertIn("APPLY_SCORE:SCORE_AGGREGATE", text)
        self.assertIn("cpus = 1", text)
        self.assertIn("memory = '14.GB'", text)

    def test_beagle_failed_due_to_no_common_markers_detects_known_error(self) -> None:
        error = subprocess.CalledProcessError(
            1,
            ["java", "-jar", "beagle.jar"],
            stderr="java.lang.IllegalArgumentException: The reference and target VCF files contain no markers in common in the window:\nX:2781457-35596385\n",
        )
        self.assertTrue(beagle_failed_due_to_no_common_markers(error))

    def test_beagle_failed_due_to_no_common_markers_ignores_other_errors(self) -> None:
        error = subprocess.CalledProcessError(
            1,
            ["java", "-jar", "beagle.jar"],
            stderr="java.lang.OutOfMemoryError: Java heap space\n",
        )
        self.assertFalse(beagle_failed_due_to_no_common_markers(error))

    def test_beagle_log_has_no_common_markers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log = Path(tmpdir) / "chrX.beagle.log"
            log.write_text(
                "java.lang.IllegalArgumentException: The reference and target VCF files contain no markers in common in the window:\n",
                encoding="utf-8",
            )
            self.assertTrue(beagle_log_has_no_common_markers(log))

    def test_parse_base_tokens_accepts_single_base_tokens_only(self) -> None:
        self.assertEqual(parse_base_tokens("G/T"), {"G", "T"})
        self.assertEqual(parse_base_tokens(" C / G / T "), {"C", "G", "T"})
        self.assertEqual(parse_base_tokens("NA"), set())
        self.assertEqual(parse_base_tokens("<DEL>"), set())

    def test_complement_base_tokens(self) -> None:
        self.assertEqual(complement_base_tokens({"A", "C"}), {"T", "G"})

    def test_base_contig_name_strips_chr_prefix(self) -> None:
        self.assertEqual(base_contig_name("chr7", "chr"), "7")
        self.assertEqual(base_contig_name("7", ""), "7")

    def test_inject_missing_contig_headers_inserts_before_column_header(self) -> None:
        original = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        updated = inject_missing_contig_headers(original, ["5"])
        self.assertIn("##contig=<ID=5>\n#CHROM", updated)

    def test_inject_missing_contig_headers_is_idempotent_for_declared_contigs(self) -> None:
        original = "##fileformat=VCFv4.2\n##contig=<ID=5>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        updated = inject_missing_contig_headers(original, ["5"])
        self.assertEqual(updated, original)

    def test_parse_pvar_header_and_iter_reference_db_sites(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            pvar = Path(tmpdir) / "chr1.pvar"
            pvar.write_text(
                "\n".join(
                    [
                        "##INFO=test",
                        "#CHROM\tPOS\tID\tREF\tALT",
                        "1\t101\trs1\tA\tG",
                        "1\t102\trs2\tC\tT,G",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_pvar_header("#CHROM\tPOS\tID\tREF\tALT\n")["POS"], 1)
            records = list(iter_reference_db_sites(pvar, target_prefix="chr"))
            self.assertEqual(records[0]["chrom"], "chr1")
            self.assertEqual(records[0]["alts"], ["G"])
            self.assertEqual(records[1]["alts"], ["T", "G"])

    def test_target_has_imputable_markers_ignores_reference_only_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            vcf = root / "reference_only.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##contig=<ID=chrX>",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "chrX\t45216839\t.\tC\t.\t.\tPASS\t.\tGT\t0/0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertFalse(target_has_imputable_markers(target_vcf=vcf, shard_dir=root, shard_name="chrX"))

    def test_target_has_imputable_markers_accepts_variant_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            vcf = root / "variant.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##contig=<ID=chr5>",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "chr5\t12345\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertTrue(target_has_imputable_markers(target_vcf=vcf, shard_dir=root, shard_name="chr5"))

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

    def test_resolve_panel_resources_for_contig_uses_base_chromosome_placeholders(self) -> None:
        panel = {
            "contig_prefix": "",
            "base_chromosomes": ["1", "2", "X"],
            "chrom_vcf_pattern": "/ref/panel/{base_chrom}.vcf.gz",
            "map_pattern": "/maps/plink.{base_chrom}.map",
        }
        resolved = resolve_panel_resources_for_contig(panel, contig="chr2", target_prefix="chr")
        assert resolved is not None
        self.assertEqual(str(resolved["vcf"]), "/ref/panel/2.vcf.gz")
        self.assertEqual(str(resolved["map"]), "/maps/plink.2.map")
        self.assertEqual(resolved["panel_chrom"], "2")

    def test_load_imputation_panel_normalizes_patterns_and_chromosomes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            panel_dir = root / "imputation_panels"
            panel_dir.mkdir()
            panel_path = panel_dir / "demo.yaml"
            panel_path.write_text(
                "\n".join(
                    [
                        "panel_id: demo_panel",
                        "build: GRCh38",
                        "contig_prefix: ''",
                        "chromosomes:",
                        "  - 1",
                        "  - X",
                        "chrom_vcf_pattern: ./panel/chr{base_chrom}.vcf.gz",
                        "map_pattern: ./maps/plink.{base_chrom}.map",
                    ]
                ),
                encoding="utf-8",
            )
            manifest_path = root / "sample.yaml"
            manifest_path.write_text("sample_id: DEMO\n", encoding="utf-8")
            context = {
                "manifest_path": str(manifest_path),
                "manifest": {
                    "reference": {"build": "GRCh38"},
                    "analysis": {"prs": {"imputation": {"reference_panel": str(panel_path)}}},
                },
                "contig_prefix": "chr",
            }

            loaded = load_imputation_panel(context)
            self.assertEqual(loaded["panel_id"], "demo_panel")
            self.assertEqual(loaded["build"], "GRCh38")
            self.assertEqual(loaded["contig_prefix"], "")
            self.assertEqual(loaded["base_chromosomes"], ["1", "X"])
            self.assertTrue(str(panel_dir) in loaded["manifest_path"])
            self.assertTrue(loaded["chrom_vcf_pattern"].endswith("/panel/chr{base_chrom}.vcf.gz"))


if __name__ == "__main__":
    unittest.main()
