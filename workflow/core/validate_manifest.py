#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys

from workflow.core.pipeline_lib import REPO_ROOT, load_yaml, validate_payload, write_json


DEFAULTS = {
    "reuse": {},
    "runtime": {
        "profile": "docker",
        "max_memory": "14.GB",
        "max_cpus": 8,
        "min_free_gb": 40,
    },
    "outputs": {
        "base_dir": "output/results",
        "cache_dir": "output/cache",
        "run_id": "latest",
    },
}

MONOGENIC_DEFAULTS = {
    "enabled": True,
    "mode": "moderate",
    "max_af": 0.05,
    "vep_forks": 8,
}

PRS_DEFAULTS = {
    "enabled": True,
    "run_ancestry": False,
    "reference_db": None,
    "min_overlap": 0.75,
}


def load_presets() -> dict[str, list[str]]:
    presets_path = REPO_ROOT / "manifests" / "score_presets.yaml"
    data = load_yaml(presets_path) or {}
    return {str(key): list(value or []) for key, value in data.items()}


def apply_defaults(manifest: dict) -> dict:
    out = dict(manifest)

    for section, defaults in DEFAULTS.items():
        merged = dict(defaults)
        merged.update(out.get(section, {}) or {})
        out[section] = merged

    analysis = dict(out.get("analysis", {}) or {})
    monogenic = dict(MONOGENIC_DEFAULTS)
    monogenic.update(analysis.get("monogenic", {}) or {})
    prs = dict(PRS_DEFAULTS)
    prs.update(analysis.get("prs", {}) or {})
    analysis["monogenic"] = monogenic
    analysis["prs"] = prs
    out["analysis"] = analysis

    reference = dict(out.get("reference", {}) or {})
    reference.setdefault("autodownload", False)
    reference.setdefault("remote_decode_probe", False)
    out["reference"] = reference

    return out


def expand_score_preset(manifest: dict) -> dict:
    prs = manifest["analysis"]["prs"]
    if not prs.get("enabled"):
        return manifest

    if prs.get("pgs_ids"):
        prs["pgs_ids"] = [str(x).strip().upper() for x in prs["pgs_ids"]]
        return manifest

    preset_name = prs.get("score_preset")
    if not preset_name:
        return manifest

    presets = load_presets()
    if preset_name not in presets:
        raise SystemExit(f"Unknown PRS score_preset: {preset_name}")

    prs["pgs_ids"] = presets[preset_name]
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--schema", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    manifest_path = pathlib.Path(args.manifest).resolve()
    schema = load_yaml(args.schema)
    manifest = load_yaml(manifest_path)
    if not isinstance(manifest, dict):
        raise SystemExit("Manifest must be a YAML object at the top level")

    manifest = apply_defaults(manifest)
    manifest = expand_score_preset(manifest)

    errors = validate_payload(manifest, schema)
    if errors:
        raise SystemExit("Manifest validation failed:\n- " + "\n- ".join(errors))

    payload = {
        "manifest_version": 1,
        "manifest_path": str(manifest_path),
        "manifest": manifest,
    }
    write_json(args.output, payload)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        print(f"Manifest validation crashed: {exc}", file=sys.stderr)
        raise
