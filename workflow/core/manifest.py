#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Any

from workflow.core.runtime import REPO_ROOT, load_json, load_yaml, normalize_path, validate_payload, write_json


DEFAULTS = {
    "reuse": {},
    "runtime": {
        "profile": "docker",
        "max_memory": "14.GB",
        "max_cpus": 8,
        "min_free_gb": 40,
        "linux_cache_root": "~/.cache/genomics",
    },
    "outputs": {
        "results_dir": "output/results",
        "run_id": "latest",
    },
}

MONOGENIC_DEFAULTS = {
    "enabled": False,
    "mode": "moderate",
    "max_af": 0.05,
    "vep_forks": 8,
    "workers": None,
    "pairhmm_threads": None,
}

PRS_DEFAULTS = {
    "enabled": True,
    "run_ancestry": False,
    "reference_db": "pgsc_HGDP+1kGP_v1",
    "min_overlap": 0.75,
    "sample_sex": None,
    "workers": None,
    "pairhmm_threads": None,
    "min_depth": 10,
    "min_gq": 20,
    "imputation": {
        "enabled": False,
        "backend": "beagle",
        "reference_panel": None,
        "info_threshold": 0.8,
        "dosage_field": "DS",
    },
}


def derive_execution_mode(manifest: dict[str, Any]) -> str:
    monogenic_enabled = bool(((manifest.get("analysis") or {}).get("monogenic") or {}).get("enabled"))
    prs_enabled = bool(((manifest.get("analysis") or {}).get("prs") or {}).get("enabled"))
    if monogenic_enabled and prs_enabled:
        return "both"
    if monogenic_enabled:
        return "monogenic-only"
    if prs_enabled:
        return "prs-only"
    raise SystemExit("At least one analysis stage must be enabled: monogenic, prs, or both.")


def load_presets(manifest_path: pathlib.Path | None = None) -> dict[str, list[str]]:
    candidates = [REPO_ROOT / "manifests" / "examples" / "score_presets.yaml"]
    if manifest_path is not None:
        candidates.append(manifest_path.parent / "score_presets.yaml")

    merged: dict[str, list[str]] = {}
    for candidate in candidates:
        if not candidate.exists():
            continue
        payload = load_yaml(candidate) or {}
        merged.update({str(key): [str(item).strip().upper() for item in (value or [])] for key, value in payload.items()})
    return merged


def apply_defaults(manifest: dict[str, Any]) -> dict[str, Any]:
    out = dict(manifest)
    out["inputs"] = dict(out.get("inputs") or {})
    out["inputs"].setdefault("protected_roots", [])

    for section, defaults in DEFAULTS.items():
        merged = dict(defaults)
        merged.update(out.get(section, {}) or {})
        out[section] = merged

    analysis = dict(out.get("analysis", {}) or {})

    monogenic = dict(MONOGENIC_DEFAULTS)
    monogenic.update(analysis.get("monogenic", {}) or {})
    analysis["monogenic"] = monogenic

    prs = dict(PRS_DEFAULTS)
    prs.update(analysis.get("prs", {}) or {})
    if not prs.get("reference_db"):
        prs["reference_db"] = PRS_DEFAULTS["reference_db"]
    prs_imputation = dict(PRS_DEFAULTS["imputation"])
    prs_imputation.update((analysis.get("prs", {}) or {}).get("imputation", {}) or {})
    prs["imputation"] = prs_imputation
    analysis["prs"] = prs

    out["analysis"] = analysis

    reference = dict(out.get("reference", {}) or {})
    out["reference"] = reference

    return out


def expand_score_preset(manifest: dict[str, Any], manifest_path: pathlib.Path | None = None) -> dict[str, Any]:
    prs = manifest["analysis"]["prs"]
    if not prs.get("enabled"):
        return manifest

    if prs.get("pgs_ids"):
        prs["pgs_ids"] = [str(item).strip().upper() for item in prs["pgs_ids"]]
        return manifest

    preset_name = prs.get("score_preset")
    if not preset_name:
        return manifest

    presets = load_presets(manifest_path)
    if preset_name not in presets:
        raise SystemExit(f"Unknown PRS score_preset: {preset_name}")

    prs["pgs_ids"] = presets[preset_name]
    return manifest


def validate_manifest_file(manifest_path: pathlib.Path, schema_path: pathlib.Path) -> dict[str, Any]:
    manifest = load_yaml(manifest_path)
    if not isinstance(manifest, dict):
        raise SystemExit("Manifest must be a YAML object at the top level")

    manifest = apply_defaults(manifest)
    manifest = expand_score_preset(manifest, manifest_path)

    schema = load_yaml(schema_path)
    errors = validate_payload(manifest, schema)
    if errors:
        raise SystemExit("Manifest validation failed:\n- " + "\n- ".join(errors))

    derive_execution_mode(manifest)

    return manifest


def resolve_linux_cache_root(manifest_path: pathlib.Path | None = None) -> str:
    if manifest_path is None:
        return normalize_path(DEFAULTS["runtime"]["linux_cache_root"]) or str(pathlib.Path.home() / ".cache" / "genomics")

    raw_manifest = load_yaml(manifest_path)
    if not isinstance(raw_manifest, dict):
        raise SystemExit("Manifest must be a YAML object at the top level")

    manifest = apply_defaults(raw_manifest)
    cache_root = normalize_path(manifest["runtime"]["linux_cache_root"], base_dir=manifest_path.parent)
    from workflow.core.reference import validate_work_roots
    validate_work_roots(
        runtime_root=cache_root,
        results_root=normalize_path(manifest["outputs"]["results_dir"], base_dir=REPO_ROOT),
        input_paths=[normalize_path(manifest["inputs"][key], base_dir=manifest_path.parent) for key in ("cram", "crai")],
        protected_roots=[normalize_path(path, base_dir=manifest_path.parent) for path in manifest["inputs"].get("protected_roots", [])],
    )
    return cache_root


def cmd_validate(args: argparse.Namespace) -> int:
    manifest_path = pathlib.Path(args.manifest).resolve()
    schema_path = pathlib.Path(args.schema).resolve()
    manifest = validate_manifest_file(manifest_path, schema_path)
    payload = {
        "manifest_version": 2,
        "manifest_path": str(manifest_path),
        "schema_path": str(schema_path),
        "manifest": manifest,
        "execution_mode": derive_execution_mode(manifest),
    }
    write_json(args.output, payload)
    return 0


def cmd_print_project_cache_root(args: argparse.Namespace) -> int:
    from workflow.core.reference import project_cache_id

    manifest_path = pathlib.Path(args.manifest).resolve()
    cache_root = pathlib.Path(resolve_linux_cache_root(manifest_path))
    print(cache_root / project_cache_id())
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    from workflow.core.reference import project_cache_id
    import json

    manifest_path = pathlib.Path(args.manifest).resolve()
    run_dir = pathlib.Path(resolve_linux_cache_root(manifest_path)) / project_cache_id() / "run"
    state = run_dir / "status.json"
    payload = load_json(state) if state.exists() else {"status": "not_started"}
    payload["status_file"] = str(state)
    # Read existing stage pointers; status must not fingerprint or scan genomic inputs.
    context_path = run_dir / "context.json"
    if payload.get("manifest_path") == str(manifest_path) and context_path.exists():
        context = load_json(context_path)
        prs_status = pathlib.Path(context["cache"]["prs_root"]) / "stage_status.json"
        if payload.get("phase") == "prs" and prs_status.exists():
            payload["prs"] = load_json(prs_status)
    elif payload.get("manifest_path") not in (None, str(manifest_path)):
        payload = {"status": "not_started", "status_file": str(state),
                   "message": "The cache currently records a different manifest."}
    if args.output:
        write_json(args.output, payload)
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("--manifest", required=True)
    validate_parser.add_argument("--schema", required=True)
    validate_parser.add_argument("--output", required=True)
    validate_parser.set_defaults(func=cmd_validate)

    project_cache_root_parser = subparsers.add_parser("print-project-cache-root")
    project_cache_root_parser.add_argument("--manifest", required=True)
    project_cache_root_parser.set_defaults(func=cmd_print_project_cache_root)

    status_parser = subparsers.add_parser("status")
    status_parser.add_argument("--manifest", required=True)
    status_parser.add_argument("--output")
    status_parser.set_defaults(func=cmd_status)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    raise SystemExit(args.func(args))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        print(f"Manifest command crashed: {exc}", file=sys.stderr)
        raise
