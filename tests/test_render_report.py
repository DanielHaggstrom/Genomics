from __future__ import annotations

import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from workflow.core.runtime import atomic_copy
from workflow.reporting.report import enforce_outcome, render_report


REPO_ROOT = Path(__file__).resolve().parents[1]


class RenderReportTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.cache = self.root / "cache"
        self.cache.mkdir()
        self.results = self.root / "results"
        self.context = {
            "sample_id": "TEST", "execution_mode": "prs-only",
            "manifest": {"analysis": {"monogenic": {"enabled": False}, "prs": {"enabled": True}}},
            "reference_bundle": {"status": "ready", "fasta": "/private/reference.fa"},
            "results": {"run_dir": str(self.results), "report_dir": str(self.results / "report")},
        }

    def fixture(self, name: str, content: str = "fixture\n") -> str:
        path = self.cache / name
        path.write_text(content)
        return str(path)

    def render(self, prs: dict, monogenic: dict | None = None) -> tuple[str, dict]:
        summary_path = self.root / "summary.json"
        render_report(self.context, monogenic or {"status": "skipped"}, prs, output_summary=summary_path)
        return (self.results / "report" / "report.md").read_text(), json.loads(summary_path.read_text())

    def test_complete_score_results_without_heavy_genotype_duplicates(self) -> None:
        self.context["execution_mode"] = "both"
        self.context["manifest"]["analysis"]["monogenic"]["enabled"] = True
        scores = self.cache / "scores.txt.gz"
        with gzip.open(scores, "wt") as handle:
            handle.write("sampleset\tPGS\tIID\tSUM\tDENOM\tAVG\ttrait\n")
            for index in range(1, 26):
                trait = "Example | trait" if index == 1 else ""
                handle.write(f"target\tPGS{999000 + index}_hmPOS_GRCh38\tTEST\t1.234\t196\t0.006296\t{trait}\n")
            handle.write("reference\tPGS999001_hmPOS_GRCh38\tREFERENCE_SAMPLE\t100\t2\t50\tReference\n")
        matching = self.fixture("match_summary.csv", "dataset,accession,match_status,count\n"
                                "target,PGS999001,matched,98\ntarget,PGS999001,excluded,1\n"
                                "target,PGS999001,unmatched,1\ntarget,PGS999001,not_best,20\n"
                                "reference,PGS999001,matched,75\nreference,PGS999001,unmatched,25\n")
        dense = self.fixture("dense.vcf.gz")
        imputed = self.fixture("imputed.vcf.gz")
        annotated = self.fixture("annotated.vcf.gz")
        monogenic = {"status": "ok", "summary": {"priority_rows_after_dedup": 2, "total_variant_records_seen": 5},
                     "artifacts": {"annotated_vcf": annotated, "priority_variants_tsv": self.fixture("priority.tsv")},
                     "import_details": {"source": "/originals/source.vcf.gz", "trusted_provenance": False}}
        prs = {
            "status": "ok", "dense_genotypes": {"artifacts": {"dense_vcf": dense}},
            "imputed_genotypes": {"artifacts": {"imputed_vcf": imputed}},
            "artifacts": {"aggregated_scores": str(scores), "match_summary": matching},
            "summary": {"num_requested_scores": 26, "genotype_qc": {"total_sites": 100, "callable_sites": 98, "missing_sites": 2},
                        "overlap_failures": [{"score": "PGS999999_hmPOS_GRCh38", "match_percent": "51.00"}]},
        }
        text, summary = self.render(prs, monogenic)
        self.assertIn("Genome Report: TEST", text)
        self.assertLess(text.index("## PRS Results"), text.index("## Execution Details"))
        self.assertIn("Calculated 25 score rows from 26 requested scores", text)
        self.assertIn("1 reference-panel rows", text)
        self.assertNotIn("REFERENCE_SAMPLE", text)
        self.assertIn("| PGS999025_hmPOS_GRCh38 |", text)
        self.assertEqual(text.count("| PGS999001_hmPOS_GRCh38 |"), 1)
        self.assertIn("| Example \\| trait | TEST | 1.234 | 196 | 0.006296 | 98 / 100 | 98.00 |", text)
        self.assertIn("| PGS999999_hmPOS_GRCh38 | 51.00 |", text)
        self.assertIn("No reference calibration was requested", text)
        self.assertIn("are not absolute disease probabilities", text)
        self.assertIn("Per-score genotype missingness was not reported", text)
        self.assertIn("98 callable, 2 missing out of 100 sites", text)
        self.assertIn("not per-score coverage", text)
        self.assertIn("## Monogenic", text)
        self.assertEqual(summary["overall_status"], "ok")
        self.assertTrue((self.results / "prs" / "aggregated_scores.txt.gz").is_file())
        self.assertTrue((self.results / "monogenic" / "priority_variants.tsv").is_file())
        self.assertEqual(list(self.results.rglob("*.vcf.gz")), [])
        self.assertEqual(list(self.results.rglob("*.tbi")), [])
        self.assertFalse((self.results / "manifests").exists())
        self.assertNotIn("/private/reference.fa", text)
        provenance = json.loads((self.results / "report" / "artifact_manifest.json").read_text())
        self.assertEqual(provenance["prs_stage"]["dense_genotypes"]["artifacts"]["dense_vcf"], dense)
        self.assertEqual(provenance["prs_stage"]["imputed_genotypes"]["artifacts"]["imputed_vcf"], imputed)
        self.assertEqual(provenance["monogenic_stage"]["artifacts"]["annotated_vcf"], annotated)
        self.assertEqual(provenance["monogenic_stage"]["import_details"], monogenic["import_details"])
        # Restarting report publication must not recopy unchanged companions.
        with patch("workflow.reporting.report.atomic_copy", wraps=atomic_copy) as copier:
            self.render(prs, monogenic)
        copier.assert_not_called()

    def test_cli_reports_failed_scores_and_enforces_failure(self) -> None:
        prs = {"status": "failed", "message": "No scores met the matching threshold.",
               "diagnostics": {"overlap_failures": [{"score": "PGS999001", "match_percent": "67.56"}]}}
        paths = [self.root / name for name in ("context.json", "monogenic.json", "prs.json")]
        for path, payload in zip(paths, [self.context, {"status": "skipped"}, prs]):
            path.write_text(json.dumps(payload))
        summary_path = self.root / "summary.json"
        result = subprocess.run([sys.executable, "-m", "workflow.reporting.report", "render", "--context", str(paths[0]),
                                 "--monogenic-stage", str(paths[1]), "--prs-stage", str(paths[2]), "--output", str(summary_path)],
                                cwd=REPO_ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        text = (self.results / "report" / "report.md").read_text()
        self.assertIn("No final PRS scores are available", text)
        self.assertIn("| PGS999001 | 67.56 |", text)
        self.assertNotIn("## Monogenic", text)
        summary = json.loads(summary_path.read_text())
        self.assertEqual(summary["failed_enabled_stages"], ["prs"])
        with self.assertRaises(SystemExit):
            enforce_outcome(summary)

    def test_incomplete_stage_cannot_be_reported_as_success(self) -> None:
        _, summary = self.render({"status": "interrupted"})
        self.assertEqual(summary["overall_status"], "failed")
        self.assertEqual(summary["failed_enabled_stages"], ["prs"])


if __name__ == "__main__":
    unittest.main()
