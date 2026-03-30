from __future__ import annotations

import pathlib
import tempfile
import unittest

from workflow.core.validate_manifest import apply_defaults, expand_score_preset


class ValidateManifestTests(unittest.TestCase):
    def test_defaults_and_score_preset_expand(self) -> None:
        manifest = {
            "sample_id": "TEST_SAMPLE",
            "inputs": {
                "cram": "/mnt/c/data/sample.cram",
                "crai": "/mnt/c/data/sample.cram.crai",
            },
            "reference": {
                "build": "GRCh38",
                "bundle_id": "bundle",
            },
            "analysis": {
                "monogenic": {
                    "enabled": True,
                },
                "prs": {
                    "enabled": True,
                    "score_preset": "demo",
                },
            },
        }

        manifest = apply_defaults(manifest)
        manifest = expand_score_preset(manifest)

        self.assertEqual(manifest["runtime"]["profile"], "docker")
        self.assertEqual(manifest["outputs"]["run_id"], "latest")
        self.assertFalse(manifest["reference"]["remote_decode_probe"])
        self.assertEqual(manifest["analysis"]["monogenic"]["mode"], "moderate")
        self.assertEqual(manifest["analysis"]["prs"]["pgs_ids"], ["PGS999001"])

    def test_local_score_presets_override_examples(self) -> None:
        manifest = {
            "sample_id": "TEST_SAMPLE",
            "inputs": {
                "cram": "/mnt/c/data/sample.cram",
                "crai": "/mnt/c/data/sample.cram.crai",
            },
            "reference": {
                "build": "GRCh38",
                "bundle_id": "bundle",
            },
            "analysis": {
                "monogenic": {
                    "enabled": True,
                },
                "prs": {
                    "enabled": True,
                    "score_preset": "demo",
                },
            },
        }

        manifest = apply_defaults(manifest)
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest_path = pathlib.Path(tmpdir) / "sample.yaml"
            preset_path = pathlib.Path(tmpdir) / "score_presets.yaml"
            manifest_path.write_text("sample_id: TEST_SAMPLE\n", encoding="utf-8")
            preset_path.write_text("demo:\n  - PGS999002\n", encoding="utf-8")
            manifest = expand_score_preset(manifest, manifest_path)

        self.assertEqual(manifest["analysis"]["prs"]["pgs_ids"], ["PGS999002"])


if __name__ == "__main__":
    unittest.main()
