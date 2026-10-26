from __future__ import annotations

import pathlib
import tempfile
import unittest

from workflow.core.manifest import apply_defaults, derive_execution_mode, expand_score_preset, resolve_linux_cache_root, validate_manifest_file


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


class ValidateManifestTests(unittest.TestCase):
    def test_launcher_rejects_cache_in_originals_before_creating_directories(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            originals = root / "originals"
            originals.mkdir()
            manifest = root / "sample.yaml"
            manifest.write_text(f"inputs:\n  cram: {originals}/sample.cram\n  crai: {originals}/sample.cram.crai\nruntime:\n  linux_cache_root: {originals}/cache\n")
            with self.assertRaisesRegex(ValueError, "protected input"):
                resolve_linux_cache_root(manifest)
            self.assertEqual(list(originals.iterdir()), [])

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
        self.assertEqual(manifest["runtime"]["linux_cache_root"], "~/.cache/genomics")
        self.assertEqual(manifest["outputs"]["results_dir"], "output/results")
        self.assertEqual(manifest["outputs"]["run_id"], "latest")
        self.assertEqual(manifest["analysis"]["monogenic"]["mode"], "moderate")
        self.assertIsNone(manifest["analysis"]["monogenic"]["workers"])
        self.assertIsNone(manifest["analysis"]["monogenic"]["pairhmm_threads"])
        self.assertIsNone(manifest["analysis"]["prs"]["workers"])
        self.assertIsNone(manifest["analysis"]["prs"]["pairhmm_threads"])
        self.assertEqual(manifest["analysis"]["prs"]["reference_db"], "pgsc_HGDP+1kGP_v1")
        self.assertFalse(manifest["analysis"]["prs"]["imputation"]["enabled"])
        self.assertEqual(manifest["analysis"]["prs"]["imputation"]["backend"], "beagle")
        self.assertEqual(manifest["analysis"]["prs"]["imputation"]["dosage_field"], "DS")
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

    def test_derive_execution_mode(self) -> None:
        both = apply_defaults(
            {
                "sample_id": "TEST",
                "inputs": {"cram": "a", "crai": "b"},
                "reference": {"build": "GRCh38", "bundle_id": "bundle"},
                "analysis": {"monogenic": {"enabled": True}, "prs": {"enabled": True, "pgs_ids": ["PGS999001"]}},
            }
        )
        self.assertEqual(derive_execution_mode(both), "both")

        monogenic_only = apply_defaults(
            {
                "sample_id": "TEST",
                "inputs": {"cram": "a", "crai": "b"},
                "reference": {"build": "GRCh38", "bundle_id": "bundle"},
                "analysis": {"monogenic": {"enabled": True}, "prs": {"enabled": False}},
            }
        )
        self.assertEqual(derive_execution_mode(monogenic_only), "monogenic-only")

        prs_only = apply_defaults(
            {
                "sample_id": "TEST",
                "inputs": {"cram": "a", "crai": "b"},
                "reference": {"build": "GRCh38", "bundle_id": "bundle"},
                "analysis": {"monogenic": {"enabled": False}, "prs": {"enabled": True, "pgs_ids": ["PGS999001"]}},
            }
        )
        self.assertEqual(derive_execution_mode(prs_only), "prs-only")

    def test_validate_manifest_rejects_both_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            manifest_path = tmp / "sample.yaml"
            manifest_path.write_text(
                "\n".join(
                    [
                        "sample_id: TEST_SAMPLE",
                        "inputs:",
                        "  cram: /mnt/c/data/sample.cram",
                        "  crai: /mnt/c/data/sample.cram.crai",
                        "reference:",
                        "  build: GRCh38",
                        "  bundle_id: bundle",
                        "analysis:",
                        "  monogenic:",
                        "    enabled: false",
                        "  prs:",
                        "    enabled: false",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            schema_path = REPO_ROOT / "manifests" / "schema" / "run_manifest.schema.json"
            with self.assertRaises(SystemExit) as ctx:
                validate_manifest_file(manifest_path, schema_path)
        self.assertIn("At least one analysis stage must be enabled", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
