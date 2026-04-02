#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Any

from workflow.core.runtime import REPO_ROOT, hash_strings, load_json, load_yaml, normalize_path, sample_file_fingerprint, validate_payload, write_json


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
    "enabled": True,
    "mode": "moderate",
    "max_af": 0.05,
    "vep_forks": 8,
    "workers": None,
    "pairhmm_threads": None,
}

PRS_DEFAULTS = {
    "enabled": True,
    "run_ancestry": False,
    "reference_db": None,
    "min_overlap": 0.75,
    "sample_sex": None,
    "workers": None,
    "pairhmm_threads": None,
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
    analysis["prs"] = prs

    out["analysis"] = analysis

    reference = dict(out.get("reference", {}) or {})
    reference.setdefault("autodownload", False)
    reference.setdefault("remote_decode_probe", False)
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
    return normalize_path(manifest["runtime"]["linux_cache_root"], base_dir=manifest_path.parent) or str(pathlib.Path.home() / ".cache" / "genomics")


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


def cmd_print_linux_cache_root(args: argparse.Namespace) -> int:
    manifest_path = pathlib.Path(args.manifest).resolve()
    print(resolve_linux_cache_root(manifest_path))
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
    schema_path = REPO_ROOT / "manifests" / "schema" / "run_manifest.schema.json"
    manifest = validate_manifest_file(manifest_path, schema_path)

    linux_cache_root = pathlib.Path(resolve_linux_cache_root(manifest_path))
    project_cache_root = linux_cache_root / project_cache_id()
    normalized_cram = pathlib.Path(normalize_path(manifest["inputs"]["cram"], base_dir=manifest_path.parent))
    normalized_crai = pathlib.Path(normalize_path(manifest["inputs"]["crai"], base_dir=manifest_path.parent))

    payload: dict[str, Any] = {
        "manifest_path": str(manifest_path),
        "execution_mode": derive_execution_mode(manifest),
        "project_cache_root": str(project_cache_root),
        "paths": {
            "nextflow_log": str(project_cache_root / "nextflow" / "nextflow.log"),
        },
        "inputs": {
            "cram": str(normalized_cram),
            "crai": str(normalized_crai),
            "cram_exists": normalized_cram.exists(),
            "crai_exists": normalized_crai.exists(),
        },
    }

    if normalized_cram.exists() and normalized_crai.exists():
        def load_status_json(path: pathlib.Path) -> dict[str, Any] | None:
            if not path.exists():
                return None
            try:
                return load_json(path)
            except Exception as exc:
                return {"status": "unreadable", "message": str(exc), "path": str(path)}

        cram_fingerprint = sample_file_fingerprint(normalized_cram)
        crai_fingerprint = sample_file_fingerprint(normalized_crai)
        input_fingerprint_id = hash_strings([cram_fingerprint["id"], crai_fingerprint["id"]])
        sample_cache_root = project_cache_root / "samples" / input_fingerprint_id
        monogenic_status_path = sample_cache_root / "monogenic" / "stage_manifest.json"
        prs_status_path = sample_cache_root / "prs" / "stage_status.json"
        payload["input_fingerprint_id"] = input_fingerprint_id
        payload["paths"].update(
            {
                "monogenic_stage_status": str(monogenic_status_path),
                "prs_stage_status": str(prs_status_path),
            }
        )
        payload["stage_statuses"] = {
            "monogenic": load_status_json(monogenic_status_path),
            "prs": load_status_json(prs_status_path),
        }

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

    cache_root_parser = subparsers.add_parser("print-linux-cache-root")
    cache_root_parser.add_argument("--manifest", required=True)
    cache_root_parser.set_defaults(func=cmd_print_linux_cache_root)

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
