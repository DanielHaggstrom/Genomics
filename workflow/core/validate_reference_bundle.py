#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import pathlib
import shutil
import urllib.request

from workflow.core.pipeline_lib import (
    REPO_ROOT,
    load_json,
    load_yaml,
    normalize_path,
    run_command,
    stage_payload,
    write_json,
)
from workflow.core.reference_lib import compare_sequence_catalogs, parse_sequence_dict


def resolve_bundle_manifest(context: dict) -> pathlib.Path | None:
    ref_cfg = context["manifest"]["reference"]
    if ref_cfg.get("bundle_manifest"):
        candidate = pathlib.Path(normalize_path(ref_cfg["bundle_manifest"], base_dir=REPO_ROOT))
        return candidate if candidate.exists() else None

    bundle_id = ref_cfg["bundle_id"]
    candidates = [
        REPO_ROOT / "manifests" / "local" / "reference_bundles" / f"{bundle_id}.yaml",
        REPO_ROOT / "manifests" / "examples" / "reference_bundles" / f"{bundle_id}.yaml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def maybe_download(bundle: dict, bundle_cache_dir: pathlib.Path) -> dict:
    download_cfg = bundle.get("download") or {}
    url = download_cfg.get("url")
    if not url:
        return bundle

    fasta = pathlib.Path(bundle["fasta"]) if bundle.get("fasta") else None
    if fasta and fasta.exists():
        return bundle

    bundle_cache_dir.mkdir(parents=True, exist_ok=True)
    archive_path = bundle_cache_dir / pathlib.Path(url).name
    urllib.request.urlretrieve(url, archive_path)

    if archive_path.suffix == ".gz":
        fasta_path = bundle_cache_dir / archive_path.stem
        with gzip.open(archive_path, "rb") as src, open(fasta_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        bundle["fasta"] = str(fasta_path)
    else:
        bundle["fasta"] = str(archive_path)

    fasta = pathlib.Path(bundle["fasta"])
    bundle.setdefault("fai", str(fasta.parent / f"{fasta.name}.fai"))
    bundle.setdefault("dict", str(fasta.with_suffix(".dict")))
    return bundle


def ensure_bundle_files(bundle: dict) -> tuple[pathlib.Path | None, pathlib.Path | None, pathlib.Path | None, list[str]]:
    fasta = pathlib.Path(bundle["fasta"]) if bundle.get("fasta") else None
    fai = pathlib.Path(bundle["fai"]) if bundle.get("fai") else (pathlib.Path(str(fasta) + ".fai") if fasta else None)
    seq_dict = pathlib.Path(bundle["dict"]) if bundle.get("dict") else (fasta.with_suffix(".dict") if fasta else None)

    missing: list[str] = []
    if fasta is None or not fasta.exists():
        missing.append("fasta")
    if fai is None or not fai.exists():
        if fasta and fasta.exists():
            run_command(["samtools", "faidx", str(fasta)], check=True)
        if fai is not None and not fai.exists():
            missing.append("fai")
    if seq_dict is None or not seq_dict.exists():
        if fasta and fasta.exists():
            run_command(["samtools", "dict", "-o", str(seq_dict), str(fasta)], check=True)
        if seq_dict is not None and not seq_dict.exists():
            missing.append("dict")
    return fasta, fai, seq_dict, missing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--reference-probe", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    probe = load_json(args.reference_probe)
    required = context["stage_requirements"]["any_reference_required"]
    bundle_id = context["manifest"]["reference"]["bundle_id"]
    bundle_path = resolve_bundle_manifest(context)
    bundle_cache_dir = pathlib.Path(context["cache_dir"]) / "references" / bundle_id

    if probe.get("status") != "ok":
        payload = stage_payload(
            "reference_bundle",
            "failed" if required else "skipped",
            required=required,
            message="CRAM reference probe did not succeed, so bundle validation could not continue.",
            probe_status=probe.get("status"),
        )
        write_json(args.output, payload)
        return

    if bundle_path is None:
        status = "failed" if required else "skipped"
        message = f"No bundle manifest found for reference bundle_id={bundle_id}"
        write_json(
            args.output,
            stage_payload(
                "reference_bundle",
                status,
                required=required,
                message=message,
                probe_summary=probe.get("summary"),
                remote_decode_probe=probe.get("remote_decode_probe"),
            ),
        )
        return

    bundle = load_yaml(bundle_path) or {}
    if bundle.get("fasta"):
        bundle["fasta"] = normalize_path(bundle["fasta"], base_dir=bundle_path.parent)
    if bundle.get("fai"):
        bundle["fai"] = normalize_path(bundle["fai"], base_dir=bundle_path.parent)
    if bundle.get("dict"):
        bundle["dict"] = normalize_path(bundle["dict"], base_dir=bundle_path.parent)

    if context["manifest"]["reference"].get("autodownload"):
        bundle = maybe_download(bundle, bundle_cache_dir)

    try:
        fasta, fai, seq_dict, missing = ensure_bundle_files(bundle)
    except Exception as exc:
        payload = stage_payload(
            "reference_bundle",
            "failed" if required else "skipped",
            required=required,
            bundle_id=bundle_id,
            bundle_manifest=str(bundle_path),
            message=f"Reference bundle preparation failed: {exc}",
        )
        write_json(args.output, payload)
        return

    if missing:
        status = "failed" if required else "skipped"
        payload = stage_payload(
            "reference_bundle",
            status,
            required=required,
            bundle_id=bundle_id,
            bundle_manifest=str(bundle_path),
            message=(
                f"Reference bundle {bundle_id} is missing required files ({', '.join(missing)}). "
                "CRAM-based stages need a validated local bundle."
            ),
            remote_decode_probe=probe.get("remote_decode_probe"),
        )
        write_json(args.output, payload)
        return

    cram_sq_entries = (probe.get("cram_header") or {}).get("sq_entries") or []
    candidate_sq_entries = parse_sequence_dict(seq_dict)
    comparison = compare_sequence_catalogs(cram_sq_entries, candidate_sq_entries)
    expected_prefix = bundle.get("expected_contig_prefix", "")
    observed_prefix = (probe.get("summary") or {}).get("contig_prefix", "")
    prefix_matches = expected_prefix in {"", observed_prefix}

    if comparison["exact_for_cram"]:
        payload = stage_payload(
            "reference_bundle",
            "ready",
            required=required,
            bundle_id=bundle_id,
            bundle_manifest=str(bundle_path),
            fasta=str(fasta),
            fai=str(fai),
            dict=str(seq_dict),
            expected_contig_prefix=expected_prefix,
            observed_contig_prefix=observed_prefix,
            prefix_matches=prefix_matches,
            comparison=comparison,
            provenance_hints=probe.get("provenance_hints"),
            remote_decode_probe=probe.get("remote_decode_probe"),
        )
        write_json(args.output, payload)
        return

    status = "failed" if required else "skipped"
    payload = stage_payload(
        "reference_bundle",
        status,
        required=required,
        bundle_id=bundle_id,
        bundle_manifest=str(bundle_path),
        fasta=str(fasta),
        fai=str(fai),
        dict=str(seq_dict),
        expected_contig_prefix=expected_prefix,
        observed_contig_prefix=observed_prefix,
        prefix_matches=prefix_matches,
        comparison=comparison,
        provenance_hints=probe.get("provenance_hints"),
        remote_decode_probe=probe.get("remote_decode_probe"),
        message="Local reference bundle did not validate against the CRAM header.",
    )
    write_json(args.output, payload)


if __name__ == "__main__":
    main()
