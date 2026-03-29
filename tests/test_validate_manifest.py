from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
