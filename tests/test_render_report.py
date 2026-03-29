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
    def test_render_unified_report_writes_markdown(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results_dir = tmp / "results" / "TEST" / "latest"
            results_dir.mkdir(parents=True, exist_ok=True)

            scores_path = tmp / "aggregated_scores.txt.gz"
            with gzip.open(scores_path, "wt", encoding="utf-8") as handle:
                handle.write("sampleset\tIID\tPGS\tSUM\n")
                handle.write("test\tTEST\tPGS000014_hmPOS_GRCh38\t1.234\n")

            context = {
                "sample_id": "TEST",
                "results_run_dir": str(results_dir),
                "cache_dir": str(tmp / "cache"),
                "manifest": {
                    "inputs": {
                        "cram": "C:\\Data\\sample.cram",
                        "crai": "C:\\Data\\sample.cram.crai",
                    },
                    "analysis": {
                        "monogenic": {"enabled": True},
                        "prs": {"enabled": True},
                    },
                    "outputs": {
                        "base_dir": "output/results",
                        "cache_dir": "output/cache",
                        "run_id": "latest",
                    },
                },
                "input_fingerprint": {"id": "fingerprint"},
                "normalized_inputs": {
                    "cram": "/mnt/c/data/sample.cram",
                    "crai": "/mnt/c/data/sample.cram.crai",
                },
                "contig_prefix": "chr",
            }
            reference = {"status": "ready"}
            reference_probe = {
                "status": "ok",
                "summary": {
                    "contig_count": 25,
                    "first_contigs": ["chr1", "chr2"],
                    "all_contigs_have_md5": True,
                },
                "provenance_hints": {
                    "aligner": "DRAGEN",
                    "ht_reference": "/refs/hg38.fa",
                    "ref_dir": "/refs/dragen/11",
                },
                "remote_decode_probe": {
                    "attempted": False,
                    "status": "not_attempted",
                    "message": "disabled",
                },
            }
            monogenic = {
                "status": "ok",
                "summary": {"priority_rows_after_dedup": "10", "total_variant_records_seen": "20"},
                "artifacts": {},
            }
            prs = {
                "status": "ok",
                "artifacts": {
                    "aggregated_scores": str(scores_path),
                },
            }

            context_path = tmp / "context.json"
            reference_probe_path = tmp / "reference_probe.json"
            reference_path = tmp / "reference.json"
            monogenic_path = tmp / "monogenic.json"
            prs_path = tmp / "prs.json"
            summary_path = tmp / "summary.json"
            for path, payload in [
                (context_path, context),
                (reference_probe_path, reference_probe),
                (reference_path, reference),
                (monogenic_path, monogenic),
                (prs_path, prs),
            ]:
                path.write_text(json.dumps(payload), encoding="utf-8")

            subprocess.run(
                [
                    "python3",
                    "-m",
                    "workflow.core.render_unified_report",
                    "--context",
                    str(context_path),
                    "--reference-probe",
                    str(reference_probe_path),
                    "--reference-stage",
                    str(reference_path),
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

            report_path = results_dir / "report" / "report.md"
            self.assertTrue(report_path.exists())
            report_text = report_path.read_text(encoding="utf-8")
            self.assertIn("Genome Report: TEST", report_text)
            self.assertIn("CRAM: `C:\\Data\\sample.cram`", report_text)
            self.assertIn("Results directory: `output/results/TEST/latest`", report_text)
            self.assertIn("CRAM reference probe: `ok`", report_text)
            self.assertIn("Local reference bundle: validated for all CRAM contigs.", report_text)


if __name__ == "__main__":
    unittest.main()
