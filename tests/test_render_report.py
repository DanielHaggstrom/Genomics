from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class RenderReportTests(unittest.TestCase):
    def test_render_report_writes_markdown(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results_dir = tmp / "results" / "TEST" / "latest"
            report_dir = results_dir / "report"
            manifests_dir = results_dir / "manifests"
            report_dir.mkdir(parents=True, exist_ok=True)
            manifests_dir.mkdir(parents=True, exist_ok=True)

            scores_path = tmp / "aggregated_scores.txt.gz"
            with gzip.open(scores_path, "wt", encoding="utf-8") as handle:
                handle.write("sampleset\tIID\tPGS\tSUM\n")
                handle.write("test\tTEST\tPGS999001_hmPOS_GRCh38\t1.234\n")

            context = {
                "sample_id": "TEST",
                "execution_mode": "both",
                "input_fingerprint": {"id": "fingerprint"},
                "inputs": {
                    "original": {
                        "cram": "C:\\Data\\sample.cram",
                        "crai": "C:\\Data\\sample.cram.crai",
                    },
                    "active": {
                        "cram": "/home/user/.cache/genomics/staging/sample.cram",
                        "crai": "/home/user/.cache/genomics/staging/sample.cram.crai",
                    },
                    "staged": {"cram": True, "crai": True},
                },
                "runtime": {
                    "linux_cache_root": "/home/user/.cache/genomics",
                },
                "results": {
                    "run_dir": str(results_dir),
                    "report_dir": str(report_dir),
                    "manifests_dir": str(manifests_dir),
                },
                "manifest": {
                    "analysis": {
                        "monogenic": {"enabled": True},
                        "prs": {"enabled": True},
                    },
                },
                "reference_probe": {
                    "status": "ok",
                    "summary": {
                        "contig_count": 25,
                        "first_contigs": ["chr1", "chr2"],
                        "contig_prefix": "chr",
                        "all_contigs_have_md5": True,
                    },
                    "provenance_hints": {
                        "aligner": "DRAGEN",
                        "ref_dir": "/refs/dragen/11",
                    },
                    "remote_decode_probe": {"status": "ok"},
                },
                "reference_bundle": {
                    "status": "ready",
                    "fasta": "/home/user/.cache/genomics/reference/hg38.fa",
                    "comparison": {"exact_dictionary_match": True, "extra_contigs": []},
                },
                "upstream_sample_name": "TEST",
            }
            monogenic = {
                "status": "ok",
                "source_mode": "called_from_cram",
                "summary": {"priority_rows_after_dedup": "10", "total_variant_records_seen": "20"},
                "artifacts": {
                    "annotated_vcf": str(tmp / "cache" / "annotated.vcf.gz"),
                    "clinvar_unique_tsv": str(tmp / "cache" / "clinvar.tsv"),
                    "lof_unique_tsv": str(tmp / "cache" / "lof.tsv"),
                    "review_unique_tsv": str(tmp / "cache" / "review.tsv"),
                    "priority_summary": str(tmp / "cache" / "priority_variants.summary.txt"),
                    "priority_variants_tsv": str(tmp / "cache" / "priority_variants.tsv"),
                },
                "import_details": {
                    "provenance": "repo_stage_manifest",
                    "source": "/external/annotated.vcf.gz",
                    "active": str(tmp / "cache" / "annotated.vcf.gz"),
                    "trusted_provenance": True,
                    "stage_manifest": "/external/annotated.stage_manifest.json",
                },
            }
            prs = {
                "status": "ok",
                "dense_genotypes": {
                    "status": "ok",
                    "summary": {
                        "union_sites": 12345,
                        "accepted_strand_flip_rows": 222,
                        "accepted_rows_by_source": {
                            "explicit_other_allele": 12000,
                            "hm_inferOtherAllele": 345,
                        },
                        "skipped_reasons": {
                            "missing_other_allele": 12,
                            "reference_not_in_score_alleles": 34,
                        },
                    },
                    "artifacts": {
                        "dense_vcf": str(tmp / "cache" / "dense.vcf.gz"),
                        "targets_vcf": str(tmp / "cache" / "targets.vcf.gz"),
                        "pre_match_audit_tsv": str(tmp / "cache" / "pre_match_audit.tsv"),
                        "pre_match_audit_json": str(tmp / "cache" / "pre_match_audit.json"),
                    },
                    "import_details": {
                        "provenance": "external_import",
                        "source": "/external/dense.vcf.gz",
                        "active": str(tmp / "cache" / "dense.vcf.gz"),
                        "trusted_provenance": False,
                    },
                },
                "imputed_genotypes": {
                    "status": "ok",
                    "source_mode": "import_imputed_prs_vcf",
                    "import_details": {
                        "provenance": "external_import",
                        "source": "/external/imputed.vcf.gz",
                        "active": str(tmp / "cache" / "imputed.vcf.gz"),
                        "trusted_provenance": False,
                    },
                    "artifacts": {
                        "imputed_vcf": str(tmp / "cache" / "imputed.vcf.gz"),
                    },
                },
                "artifacts": {
                    "aggregated_scores": str(scores_path),
                    "match_diagnostics_json": str(tmp / "cache" / "match_diagnostics.json"),
                    "overlap_failures_tsv": str(tmp / "cache" / "overlap_failures.tsv"),
                    "match_summary": str(tmp / "cache" / "match_summary.csv"),
                    "pgsc_report_html": str(tmp / "cache" / "pgsc_report.html"),
                },
            }

            cache_dir = tmp / "cache"
            cache_dir.mkdir(exist_ok=True)
            for path in [
                Path(monogenic["artifacts"]["annotated_vcf"]),
                Path(monogenic["artifacts"]["clinvar_unique_tsv"]),
                Path(monogenic["artifacts"]["lof_unique_tsv"]),
                Path(monogenic["artifacts"]["review_unique_tsv"]),
                Path(monogenic["artifacts"]["priority_summary"]),
                Path(monogenic["artifacts"]["priority_variants_tsv"]),
                Path(prs["dense_genotypes"]["artifacts"]["pre_match_audit_tsv"]),
                Path(prs["dense_genotypes"]["artifacts"]["pre_match_audit_json"]),
                Path(prs["imputed_genotypes"]["artifacts"]["imputed_vcf"]),
                Path(prs["artifacts"]["match_diagnostics_json"]),
                Path(prs["artifacts"]["overlap_failures_tsv"]),
                Path(prs["artifacts"]["match_summary"]),
                Path(prs["artifacts"]["pgsc_report_html"]),
            ]:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("test\n", encoding="utf-8")
            Path(monogenic["artifacts"]["annotated_vcf"] + ".tbi").write_text("index\n", encoding="utf-8")

            context_path = tmp / "context.json"
            monogenic_path = tmp / "monogenic.json"
            prs_path = tmp / "prs.json"
            summary_path = tmp / "summary.json"
            for path, payload in [
                (context_path, context),
                (monogenic_path, monogenic),
                (prs_path, prs),
            ]:
                path.write_text(json.dumps(payload), encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    "-m",
                    "workflow.reporting.report",
                    "render",
                    "--context",
                    str(context_path),
                    "--monogenic-stage",
                    str(monogenic_path),
                    "--prs-stage",
                    str(prs_path),
                    "--output",
                    str(summary_path),
                ],
                cwd=REPO_ROOT,
                env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
                check=True,
            )

            report_path = report_dir / "report.md"
            self.assertTrue(report_path.exists())
            report_text = report_path.read_text(encoding="utf-8")
            self.assertIn("Genome Report: TEST", report_text)
            self.assertIn("Execution mode: `both`", report_text)
            self.assertIn("CRAM input: `C:\\Data\\sample.cram`", report_text)
            self.assertIn("Active CRAM path: `/home/user/.cache/genomics/staging/sample.cram`", report_text)
            self.assertIn("Reference probe: `ok`", report_text)
            self.assertIn(f"Dense PRS VCF: `{tmp / 'cache' / 'dense.vcf.gz'}`", report_text)
            self.assertIn("Accepted rows rescued by strand complement: `222`", report_text)
            self.assertIn(f"Pre-match audit TSV: `{tmp / 'cache' / 'pre_match_audit.tsv'}`", report_text)
            self.assertIn("Accepted score rows by source: `explicit_other_allele`=12000, `hm_inferOtherAllele`=345", report_text)
            self.assertIn("Import stage manifest: `/external/annotated.stage_manifest.json`", report_text)
            self.assertIn("Dense genotype import note: external import without repo-generated stage manifest", report_text)
            self.assertIn(f"Imputed PRS VCF: `{tmp / 'cache' / 'imputed.vcf.gz'}`", report_text)
            self.assertIn("Imputed genotype import note: external import without repo-generated stage manifest", report_text)
            self.assertIn(f"Match diagnostics JSON: `{tmp / 'cache' / 'match_diagnostics.json'}`", report_text)
            self.assertTrue((results_dir / "monogenic" / "annotated.vcf.gz").exists())
            self.assertTrue((results_dir / "monogenic" / "annotated.vcf.gz.tbi").exists())
            self.assertTrue((results_dir / "prs" / "pre_match_audit.tsv").exists())
            self.assertTrue((results_dir / "prs" / "match_diagnostics.json").exists())
            self.assertTrue((results_dir / "prs" / "imputed_prs.vcf.gz").exists())

    def test_render_report_prs_only_omits_monogenic_section(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results_dir = tmp / "results" / "TEST" / "latest"
            report_dir = results_dir / "report"
            manifests_dir = results_dir / "manifests"
            report_dir.mkdir(parents=True, exist_ok=True)
            manifests_dir.mkdir(parents=True, exist_ok=True)

            context = {
                "sample_id": "TEST",
                "execution_mode": "prs-only",
                "input_fingerprint": {"id": "fingerprint"},
                "inputs": {
                    "original": {"cram": "C:\\Data\\sample.cram", "crai": "C:\\Data\\sample.cram.crai"},
                    "active": {"cram": "/cache/sample.cram", "crai": "/cache/sample.cram.crai"},
                    "staged": {"cram": True, "crai": True},
                },
                "runtime": {"linux_cache_root": "/home/user/.cache/genomics"},
                "results": {
                    "run_dir": str(results_dir),
                    "report_dir": str(report_dir),
                    "manifests_dir": str(manifests_dir),
                },
                "manifest": {
                    "analysis": {
                        "monogenic": {"enabled": False},
                        "prs": {"enabled": True},
                    },
                },
                "reference_probe": {"status": "ok", "summary": {}, "provenance_hints": {}, "remote_decode_probe": {"status": "ok"}},
                "reference_bundle": {"status": "ready"},
            }
            monogenic = {"status": "skipped", "message": "Monogenic stage disabled."}
            prs = {
                "status": "failed",
                "dense_genotypes": {
                    "status": "ok",
                    "summary": {
                        "union_sites": 100,
                        "accepted_rows_by_source": {"explicit_other_allele": 90},
                        "skipped_reasons": {"missing_other_allele": 10},
                    },
                    "artifacts": {
                        "pre_match_audit_tsv": "/cache/pre_match_audit.tsv",
                        "pre_match_audit_json": "/cache/pre_match_audit.json",
                    },
                },
                "diagnostics": {
                    "overlap_failures": [
                        {"score": "PGS999001_hmPOS_GRCh38", "match_percent": "67.56"},
                        {"score": "PGS999002_hmPOS_GRCh38", "match_percent": "53.29"},
                    ],
                    "log_path": "/cache/pgsc_calc.log",
                },
                "artifacts": {
                    "match_diagnostics_json": "/cache/match_diagnostics.json",
                    "overlap_failures_tsv": "/cache/overlap_failures.tsv",
                    "match_summary": "/cache/match_summary.csv",
                },
            }

            context_path = tmp / "context.json"
            monogenic_path = tmp / "monogenic.json"
            prs_path = tmp / "prs.json"
            summary_path = tmp / "summary.json"
            for path, payload in [(context_path, context), (monogenic_path, monogenic), (prs_path, prs)]:
                path.write_text(json.dumps(payload), encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    "-m",
                    "workflow.reporting.report",
                    "render",
                    "--context",
                    str(context_path),
                    "--monogenic-stage",
                    str(monogenic_path),
                    "--prs-stage",
                    str(prs_path),
                    "--output",
                    str(summary_path),
                ],
                cwd=REPO_ROOT,
                env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
                check=True,
            )

            report_text = (report_dir / "report.md").read_text(encoding="utf-8")
            self.assertIn("PRS Report: TEST", report_text)
            self.assertIn("Execution mode: `prs-only`", report_text)
            self.assertNotIn("## Monogenic", report_text)
            self.assertIn("## PRS", report_text)
            self.assertIn("Match threshold summary TSV: `/cache/overlap_failures.tsv`", report_text)
            self.assertIn("| PGS999001_hmPOS_GRCh38 | 67.56 |", report_text)


if __name__ == "__main__":
    unittest.main()
