from __future__ import annotations

import argparse
import fcntl
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import yaml

from workflow import run as driver
from workflow.core.manifest import cmd_status, validate_manifest_file
from workflow.core.reference import project_cache_id
from workflow.core.runtime import checkpoint_complete, record_checkpoint, write_json


class RunTests(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, dict]:
        manifest = root / "sample.yaml"
        manifest.write_text(yaml.safe_dump({
            "sample_id": "SYNTHETIC", "inputs": {"cram": str(root / "originals/sample.cram"),
                                                  "crai": str(root / "originals/sample.cram.crai")},
            "reference": {"build": "GRCh38", "bundle_id": "synthetic"},
            "analysis": {"prs": {"pgs_ids": ["PGS999001"]}},
            "runtime": {"linux_cache_root": str(root / "cache")},
            "outputs": {"results_dir": str(root / "results")},
        }))
        context = {
            "sample_id": "SYNTHETIC", "manifest_path": str(manifest),
            "manifest": validate_manifest_file(manifest, driver.REPO_ROOT / "manifests/schema/run_manifest.schema.json"),
            "cache": {"prs_root": str(root / "cache/prs")},
            "results": {"run_dir": str(root / "results"), "report_dir": str(root / "results/report")},
        }
        return manifest, context

    def test_stop_then_resume_reuses_checkpoint_and_releases_lock(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, context = self.fixture(root)
            run_dir = root / "cache" / project_cache_id() / "run"
            artifact, receipt = root / "completed-shard", root / "completed.json"
            computations = []

            def stage(argv):
                if not checkpoint_complete(receipt, [artifact]):
                    computations.append("compute")
                    artifact.write_bytes(b"completed synthetic shard")
                    record_checkpoint(receipt, [artifact])
                    raise KeyboardInterrupt
                write_json(argv[-1], {"stage": "prs", "status": "ok"})

            def prepare(args):
                write_json(args.output, context)

            with patch.object(driver, "cmd_prepare_context", side_effect=prepare), \
                 patch.object(driver.prs, "main", side_effect=stage), patch.object(driver.monogenic, "main") as mono:
                self.assertEqual(driver.run(manifest), 130)
                self.assertEqual(json.loads((run_dir / "status.json").read_text())["status"], "stopped")
                self.assertEqual(driver.run(manifest), 0)
                self.assertEqual(computations, ["compute"])
                mono.assert_not_called()
            self.assertEqual(json.loads((run_dir / "status.json").read_text())["status"], "ok")
            self.assertTrue((root / "results/report/report.md").is_file())
            # Status reads saved pointers without needing even a real CRAM fixture.
            output = root / "status.json"
            cmd_status(argparse.Namespace(manifest=str(manifest), output=str(output)))
            self.assertEqual(json.loads(output.read_text())["status"], "ok")

    def test_lock_rejects_second_writer_before_preparation(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, _ = self.fixture(root)
            project = root / "cache" / project_cache_id()
            project.mkdir(parents=True)
            with (project / "run.lock").open("a") as lock, patch.object(driver, "cmd_prepare_context") as prepare:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                self.assertEqual(driver.run(manifest), 1)
                prepare.assert_not_called()
                self.assertFalse((project / "run/status.json").exists())

    def test_failed_stage_produces_failed_report_and_failure_status(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, context = self.fixture(root)
            with patch.object(driver, "cmd_prepare_context", side_effect=lambda args: write_json(args.output, context)), \
                 patch.object(driver.prs, "main", side_effect=lambda argv: write_json(argv[-1], {"status": "failed", "message": "synthetic failure"})):
                self.assertEqual(driver.run(manifest), 1)
            summary = json.loads((root / "cache" / project_cache_id() / "run/summary.json").read_text())
            self.assertEqual(summary["failed_enabled_stages"], ["prs"])
            self.assertIn("synthetic failure", (root / "results/report/report.md").read_text())


if __name__ == "__main__":
    unittest.main()
