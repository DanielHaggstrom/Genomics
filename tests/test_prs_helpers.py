from __future__ import annotations

import gzip
import json
import subprocess
import sys
import unittest
from pathlib import Path
import tempfile
from unittest.mock import patch

from workflow.core.manifest import apply_defaults
from workflow.stages import prs

from workflow.stages.prs import (
    base_contig_name,
    beagle_failed_due_to_no_common_markers,
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
)


class FakeFasta:
    def __init__(self, sequence_by_length: dict[tuple[str, int, int], str]) -> None:
        self.sequence_by_length = sequence_by_length

    def get_sequence(self, contig: str, pos: int, sequence_length: int) -> str:
        return self.sequence_by_length[(contig, pos, sequence_length)]


class PrsHelperTests(unittest.TestCase):
    def test_scaffold_stage_materializes_reference_and_masks_uncertain_calls(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            originals = root / "originals"
            originals.mkdir()
            cram, reference = originals / "sample.cram", originals / "ref.fa"
            cram.write_bytes(b"synthetic alignment")
            reference.write_bytes(b">chr1\nACGT\n")
            before = {path: path.read_bytes() for path in (cram, reference)}
            header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            targets = root / "targets.vcf.gz"
            with gzip.open(targets, "wt") as handle:
                handle.write(header + "\nchr1\t10\t.\tA\tG\t.\t.\t.\nchr1\t20\t.\tC\tT\t.\t.\t.\n")
            calls = header + "\tFORMAT\tNEBULA_SAMPLE\nchr1\t10\t.\tA\t.\t.\t.\t.\tGT:DP:RGQ\t0/0:25:60\nchr1\t20\t.\tC\tT\t.\t.\t.\tGT:DP:GQ\t0/1:25:5\n"

            def fake_command(args, **kwargs):
                if args[0] == "docker":
                    mounts = {args[i + 1].split(":")[1]: Path(args[i + 1].split(":")[0]) for i, token in enumerate(args) if token == "-v"}
                    output = Path(args[args.index("-O") + 1])
                    destination = mounts["/" + output.parts[1]] / output.name
                    if "GenotypeGVCFs" in args:
                        self.assertIn("--include-non-variant-sites", args)
                        with gzip.open(destination, "wt") as handle:
                            handle.write(calls)
                    else:
                        destination.write_bytes(b"synthetic intermediate")
                    Path(str(destination) + ".tbi").write_bytes(b"synthetic index")
                elif args[:2] == ["tabix", "-f"]:
                    Path(args[-1] + ".tbi").write_bytes(b"synthetic index")
                elif args[:2] == ["bcftools", "sort"]:
                    with gzip.open(args[args.index("-o") + 1], "wt") as handle:
                        handle.write(Path(args[-1]).read_text())
                else:
                    Path(args[args.index("-o") + 1]).write_bytes(Path(args[-1]).read_bytes())
                return subprocess.CompletedProcess(args, 0, "", "")

            kwargs = dict(docker_bin="docker", stage_root=root / "cache", cram_path=cram,
                          reference_fasta=reference, shard={"name": "chr1", "intervals": ["chr1"]},
                          targets_vcf=targets, sample_id="TEST", pairhmm_threads=1, java_memory="1g")
            with patch.object(prs, "run_command", side_effect=fake_command) as runner:
                result = prs.build_dense_prs_shard(**kwargs)
                count = runner.call_count
                self.assertTrue(prs.build_dense_prs_shard(**kwargs)["reused"])
                self.assertEqual(runner.call_count, count)
            with gzip.open(result["dense_vcf"], "rt") as handle:
                rows = [line.split("\t") for line in handle if not line.startswith("#")]
            self.assertEqual(rows[0][3:5], ["A", "G"])
            self.assertTrue(rows[0][9].startswith("0/0:"))
            self.assertTrue(rows[1][9].startswith("./.:"))
            self.assertEqual(result["genotype_qc"]["callable_sites"], 1)
            self.assertEqual(result["genotype_qc"]["missing_sites"], 1)
            self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_conversion_cache_tracks_actual_input_kind_dosage_and_sex(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "original.vcf.gz"
            source.write_bytes(b"original-sequence-data")
            context = {"sample_id": "TEST", "host_checks": {"docker": "docker"},
                       "manifest": apply_defaults({"analysis": {"prs": {"pgs_ids": ["PGS999001"]}}})}
            commands = []

            def fake_command(args, **kwargs):
                commands.append(args)
                if args[:3] == ["bcftools", "query", "-l"]:
                    return subprocess.CompletedProcess(args, 0, "NEBULA_SAMPLE\n", "")
                if args[0] == "docker":
                    mount = next(args[i + 1] for i, item in enumerate(args) if item == "-v")
                    work = Path(mount.rsplit(":", 2)[0])
                    prefix = Path(args[args.index("--out") + 1]).name
                    for suffix in (".pgen", ".pvar.zst", ".psam"):
                        (work / (prefix + suffix)).write_bytes(b"generated")
                elif "-o" in args:
                    Path(args[args.index("-o") + 1]).write_bytes(b"generated-vcf")
                return subprocess.CompletedProcess(args, 0, "", "")

            with patch.object(prs, "ensure_docker_ready"), patch.object(prs, "run_command", side_effect=fake_command):
                first = prs.build_pgsc_input_pfile(context=context, stage_root=root / "cache", input_vcf=source, input_kind="dense_vcf")
                count = len(commands)
                repeated = prs.build_pgsc_input_pfile(context=context, stage_root=root / "cache", input_vcf=source, input_kind="dense_vcf")
                self.assertTrue(repeated["reused"])
                self.assertEqual(len(commands), count)
                imputed = prs.build_pgsc_input_pfile(context=context, stage_root=root / "cache", input_vcf=source, input_kind="imputed_vcf", dosage_field="DS")
                self.assertNotEqual(first["path_prefix"], imputed["path_prefix"])
                self.assertTrue(any('FMT/DS!="."' in command for command in commands))
                context["manifest"]["analysis"]["prs"]["sample_sex"] = "female"
                changed_sex = prs.build_pgsc_input_pfile(context=context, stage_root=root / "cache", input_vcf=source, input_kind="imputed_vcf", dosage_field="DS")
                self.assertNotEqual(imputed["path_prefix"], changed_sex["path_prefix"])
                source.write_bytes(b"different-genotype-data")
                changed_source = prs.build_pgsc_input_pfile(context=context, stage_root=root / "cache", input_vcf=source, input_kind="imputed_vcf", dosage_field="DS")
                self.assertNotEqual(changed_sex["path_prefix"], changed_source["path_prefix"])
                # Rebuilds touch generated cache only, never their input.
                self.assertEqual(source.read_bytes(), b"different-genotype-data")

    def test_imputation_threshold_and_panel_data_get_separate_cache_trees(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dense, reference, panel_vcf = (root / name for name in ("dense.vcf.gz", "ref.fa", "panel.vcf.gz"))
            for path in (dense, reference, panel_vcf):
                path.write_bytes(b"fixture")
            context = {"input_fingerprint": {"id": "sample"},
                       "manifest": apply_defaults({"analysis": {"prs": {"pgs_ids": ["PGS999001"]}}})}
            panel = {"panel_id": "panel", "manifest_path": "/fixture/panel.yaml", "manifest_fingerprint_id": "manifest"}
            shards = [{"panel_resource": {"vcf": panel_vcf}}]
            with patch.object(prs, "load_imputation_panel", return_value=panel), patch.object(prs, "build_imputation_shards", return_value=shards), patch.object(prs, "load_existing_imputation_stage", side_effect=lambda path, key: {"cache_key": key}) as lookup:
                kwargs = dict(context=context, stage_root=root / "cache", dense_info={"artifacts": {"dense_vcf": str(dense)}}, score_hash="scores", reference_fasta=reference)
                prs.build_or_reuse_imputed_vcf(**kwargs)
                first = lookup.call_args.args[0]
                context["manifest"]["analysis"]["prs"]["imputation"]["info_threshold"] = 0.9
                prs.build_or_reuse_imputed_vcf(**kwargs)
                second = lookup.call_args.args[0]
                self.assertNotEqual(first, second)
                panel_vcf.write_bytes(b"changed-panel-data")
                prs.build_or_reuse_imputed_vcf(**kwargs)
                self.assertNotEqual(second, lookup.call_args.args[0])

    def test_beagle_error_executes_real_passthrough_handler(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shard = {"name": "chr1", "intervals": ["chr1"], "imputation_mode": "impute", "panel_resource": {"panel_chrom": "chr1"}}
            error = subprocess.CalledProcessError(1, ["java"], stderr="no markers in common")
            with patch.object(prs, "extract_dense_contig_vcf", return_value=root / "target.vcf.gz"), patch.object(prs, "rename_contig_for_panel", return_value=root / "target.vcf.gz"), patch.object(prs, "target_has_imputable_markers", return_value=True), patch.object(prs, "stage_panel_resources", return_value={"panel_vcf": root / "panel.vcf.gz", "map": None}), patch.object(prs, "run_command", side_effect=error), patch.object(prs, "build_dense_passthrough_shard", return_value=(root / "passthrough.vcf.gz", False)):
                result = prs.build_imputed_prs_shard(context={"sample_id": "TEST", "runtime": {"linux_cache_root": str(root)}}, stage_root=root, dense_vcf=root / "dense.vcf.gz", reference_fasta=root / "ref.fa", shard=shard, panel={}, beagle_jar=root / "beagle.jar", java_memory="1g", beagle_threads=1, info_threshold=0.8)
            self.assertEqual(result["mode"], "passthrough")
            self.assertEqual(result["passthrough_reason"], "no_panel_overlap_window")

    def test_imported_scaffold_reaches_requested_imputation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dense, imputed, log = (root / name for name in ("dense.vcf.gz", "imputed.vcf.gz", "pgsc.log"))
            dense.write_bytes(b"direct")
            imputed.write_bytes(b"imputed")
            log.write_text("fixture failure after input selection\n")
            manifest = apply_defaults({"analysis": {"monogenic": {"enabled": False}, "prs": {"pgs_ids": ["PGS999001"], "imputation": {"enabled": True}}}})
            context = {"sample_id": "TEST", "input_fingerprint": {"id": "sample"}, "manifest": manifest,
                       "cache": {"prs_root": str(root / "cache")}, "reference_bundle": {"status": "ready", "fasta": str(root / "ref.fa")},
                       "reuse": {"prs_dense_vcf": {"active": str(dense), "original": str(dense), "fingerprint": {"id": "import"}, "sidecars": {"index": str(dense) + ".tbi"}, "provenance": "external", "trusted_provenance": False}}}
            path = root / "context.json"
            path.write_text(json.dumps(context))
            payload = {"cache_key": "imputed", "artifacts": {"imputed_vcf": str(imputed), "imputed_vcf_index": str(imputed) + ".tbi"}}
            pfile_info = {"artifacts": {}, "renamed_chr_prefix": None}
            with patch.object(sys, "argv", ["prs", "--context", str(path), "--output", str(root / "stage.json")]), patch.object(prs, "build_or_reuse_imputed_vcf", return_value=(payload, {"phase": "imputed"})) as impute, patch.object(prs, "run_pgsc_calc", return_value=(1, log, pfile_info)) as score:
                prs.main()
            self.assertEqual(impute.call_count, 1)
            self.assertEqual(score.call_args.kwargs["scoring_vcf"], imputed)
            self.assertEqual(score.call_args.kwargs["dosage_field"], "DS")

    def test_scoring_resume_keeps_existing_work_and_pins_release(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "work").mkdir()
            sentinel = root / "work" / "completed-task"
            sentinel.write_bytes(b"keep resumable work")
            manifest = apply_defaults({"analysis": {"prs": {"pgs_ids": ["PGS999001"]}}})
            context = {"sample_id": "TEST", "manifest": manifest, "runtime": {"nextflow_home": str(root / "home")}}
            info = {"path_prefix": str(root / "input"), "format": "pfile"}
            with patch.object(prs, "build_pgsc_input_pfile", return_value=info), \
                 patch.object(prs.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")), \
                 patch.object(prs, "run_command", return_value=subprocess.CompletedProcess([], 0, "", "")) as command:
                prs.run_pgsc_calc(context=context, stage_root=root, scoring_vcf=root / "fixture.vcf.gz", input_kind="dense_vcf")
            args = command.call_args.args[0]
            self.assertEqual(args[args.index("-r") + 1], prs.PGSC_CALC_REVISION)
            self.assertIn("-resume", args)
            self.assertEqual(command.call_args.kwargs["cwd"], root)
            self.assertEqual(sentinel.read_bytes(), b"keep resumable work")
            config = (root / "pgsc_calc.resources.config").read_text()
            self.assertIn("executor {", config)
            self.assertIn("cpus = 8", config)
            self.assertIn("--label genomics.scoring=", config)

    def test_scoring_interrupt_removes_only_its_labeled_containers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            context = {"sample_id": "TEST", "host_checks": {"docker": "docker"},
                       "manifest": apply_defaults({"analysis": {"prs": {"pgs_ids": ["PGS999001"]}}}),
                       "runtime": {"nextflow_home": str(root / "home")}}
            info = {"path_prefix": str(root / "input"), "format": "pfile"}
            with patch.object(prs, "build_pgsc_input_pfile", return_value=info), \
                 patch.object(prs, "run_command", side_effect=KeyboardInterrupt), \
                 patch.object(prs.subprocess, "run", side_effect=[subprocess.CompletedProcess([], 0, "owned-container\n", ""),
                                                                 subprocess.CompletedProcess([], 0, "", "")]) as cleanup:
                with self.assertRaises(KeyboardInterrupt):
                    prs.run_pgsc_calc(context=context, stage_root=root, scoring_vcf=root / "fixture.vcf.gz", input_kind="dense_vcf")
            query, removal = [call.args[0] for call in cleanup.call_args_list]
            self.assertEqual(query[:4], ["docker", "ps", "-aq", "--filter"])
            self.assertTrue(query[4].startswith("label=genomics.scoring="))
            self.assertEqual(removal, ["docker", "rm", "-f", "owned-container"])

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
            self.assertFalse(target_has_imputable_markers(vcf))

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
            self.assertTrue(target_has_imputable_markers(vcf))

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
