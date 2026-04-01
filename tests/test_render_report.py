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
                    "annotated_vcf": "/cache/annotated.vcf.gz",
                    "clinvar_unique_tsv": "/cache/clinvar.tsv",
                    "lof_unique_tsv": "/cache/lof.tsv",
                    "review_unique_tsv": "/cache/review.tsv",
                },
            }
            prs = {
                "status": "ok",
                "dense_genotypes": {
                    "status": "ok",
                    "summary": {"union_sites": 12345},
                    "artifacts": {
                        "dense_vcf": "/cache/dense.vcf.gz",
                        "targets_vcf": "/cache/targets.vcf.gz",
                    },
                },
                "artifacts": {
                    "aggregated_scores": str(scores_path),
                },
            }

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
            self.assertIn("Dense PRS VCF: `/cache/dense.vcf.gz`", report_text)

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
            prs = {"status": "failed", "dense_genotypes": {"status": "ok"}, "diagnostics": {}}

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


if __name__ == "__main__":
    unittest.main()
