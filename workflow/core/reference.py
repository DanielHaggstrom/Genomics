#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import pathlib
import re
import shutil
import sys
import urllib.request
from typing import Any

from workflow.core.runtime import (
    REPO_ROOT,
    detect_executable,
    ensure_dir,
    hash_strings,
    load_json,
    load_yaml,
    needs_staging,
    normalize_path,
    plan_contig_shards,
    run_command,
    sample_file_fingerprint,
    stage_file,
    stage_payload,
    write_json,
)


PROJECT_CACHE_VERSION = 1


def parse_tagged_header_line(line: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for token in line.rstrip().split("\t")[1:]:
        if ":" not in token:
            continue
        key, value = token.split(":", 1)
        fields[key] = value.strip()
    return fields


def parse_sq_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "name": fields.get("SN"),
                "length": int(fields["LN"]) if fields.get("LN") else None,
                "md5": fields.get("M5"),
                "uri": fields.get("UR"),
            }
        )
    return entries


def parse_program_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@PG"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "id": fields.get("ID"),
                "name": fields.get("PN"),
                "version": fields.get("VN"),
                "command_line": fields.get("CL"),
                "previous": fields.get("PP"),
                "description": fields.get("DS"),
            }
        )
    return entries


def parse_read_group_entries(header_text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for line in header_text.splitlines():
        if not line.startswith("@RG"):
            continue
        fields = parse_tagged_header_line(line)
        entries.append(
            {
                "id": fields.get("ID"),
                "sample": fields.get("SM"),
                "platform_unit": fields.get("PU"),
                "library": fields.get("LB"),
            }
        )
    return entries


def parse_sequence_dict(path: str | pathlib.Path) -> list[dict[str, Any]]:
    return parse_sq_entries(pathlib.Path(path).read_text(encoding="utf-8"))


def index_sq_entries(entries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {entry["name"]: entry for entry in entries if entry.get("name")}


def summarize_contigs(sq_entries: list[dict[str, Any]]) -> dict[str, Any]:
    names = [entry["name"] for entry in sq_entries if entry.get("name")]
    contig_prefix = "chr" if names and names[0].startswith("chr") else ""
    all_have_md5 = all(entry.get("md5") for entry in sq_entries if entry.get("name"))
    return {
        "count": len(names),
        "first_contigs": names[:10],
        "contig_prefix": contig_prefix,
        "all_contigs_have_md5": all_have_md5,
    }


def extract_reference_hints(pg_entries: list[dict[str, Any]]) -> dict[str, Any]:
    hints: dict[str, Any] = {
        "aligner": None,
        "ht_reference": None,
        "ref_dir": None,
        "variant_annotation_assembly": None,
    }
    combined = "\n".join(entry.get("command_line") or "" for entry in pg_entries)
    lowered = combined.lower()
    if "dragen" in lowered:
        hints["aligner"] = "DRAGEN"

    patterns = {
        "ht_reference": [r"--ht-reference=(\S+)", r"--ht-reference\s+(\S+)"],
        "ref_dir": [r"--ref-dir=(\S+)", r"--ref-dir\s+(\S+)"],
        "variant_annotation_assembly": [
            r"--variant-annotation-assembly=(\S+)",
            r"--variant-annotation-assembly\s+(\S+)",
        ],
    }
    for key, regexes in patterns.items():
        for regex in regexes:
            match = re.search(regex, combined)
            if match:
                hints[key] = match.group(1)
                break
    return hints


def compare_sequence_catalogs(cram_sq_entries: list[dict[str, Any]], candidate_sq_entries: list[dict[str, Any]]) -> dict[str, Any]:
    cram_map = index_sq_entries(cram_sq_entries)
    candidate_map = index_sq_entries(candidate_sq_entries)

    missing_contigs: list[str] = []
    length_mismatches: list[dict[str, Any]] = []
    md5_mismatches: list[dict[str, Any]] = []
    md5_missing_in_candidate: list[str] = []

    for name, cram_entry in cram_map.items():
        candidate_entry = candidate_map.get(name)
        if candidate_entry is None:
            missing_contigs.append(name)
            continue
        if candidate_entry.get("length") != cram_entry.get("length"):
            length_mismatches.append(
                {
                    "contig": name,
                    "cram_length": cram_entry.get("length"),
                    "candidate_length": candidate_entry.get("length"),
                }
            )
        cram_md5 = cram_entry.get("md5")
        candidate_md5 = candidate_entry.get("md5")
        if cram_md5 and not candidate_md5:
            md5_missing_in_candidate.append(name)
        elif cram_md5 and candidate_md5 and cram_md5 != candidate_md5:
            md5_mismatches.append(
                {
                    "contig": name,
                    "cram_md5": cram_md5,
                    "candidate_md5": candidate_md5,
                }
            )

    extra_contigs = sorted(name for name in candidate_map.keys() if name not in cram_map)
    exact_for_cram = not missing_contigs and not length_mismatches and not md5_mismatches and not md5_missing_in_candidate
    exact_dictionary_match = exact_for_cram and not extra_contigs and len(candidate_map) == len(cram_map)
    return {
        "compared_contigs": len(cram_map),
        "candidate_contigs": len(candidate_map),
        "exact_for_cram": exact_for_cram,
        "exact_dictionary_match": exact_dictionary_match,
        "missing_contigs": missing_contigs,
        "length_mismatches": length_mismatches,
        "md5_mismatches": md5_mismatches,
        "md5_missing_in_candidate": md5_missing_in_candidate,
        "extra_contigs": extra_contigs,
    }


def project_cache_id() -> str:
    digest = hash_strings([str(REPO_ROOT.resolve()), f"v{PROJECT_CACHE_VERSION}"])
    return f"{REPO_ROOT.name.lower()}-{digest[:12]}"


def resolve_bundle_manifest(reference_cfg: dict[str, Any]) -> pathlib.Path | None:
    if reference_cfg.get("bundle_manifest"):
        candidate = pathlib.Path(normalize_path(reference_cfg["bundle_manifest"], base_dir=REPO_ROOT))
        return candidate if candidate.exists() else None

    bundle_id = reference_cfg["bundle_id"]
    candidates = [
        REPO_ROOT / "manifests" / "local" / "reference_bundles" / f"{bundle_id}.yaml",
        REPO_ROOT / "manifests" / "examples" / "reference_bundles" / f"{bundle_id}.yaml",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def maybe_download_bundle(bundle: dict[str, Any], bundle_cache_dir: pathlib.Path) -> dict[str, Any]:
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

    if archive_path.suffix == ".gz" and archive_path.name.endswith(".fa.gz"):
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


def ensure_bundle_files(bundle: dict[str, Any]) -> tuple[pathlib.Path | None, pathlib.Path | None, pathlib.Path | None, list[str]]:
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


def stage_maybe(source_path: pathlib.Path, destination_dir: pathlib.Path, linux_cache_root: pathlib.Path) -> tuple[str, bool]:
    if needs_staging(source_path, linux_cache_root):
        fingerprint = sample_file_fingerprint(source_path)
        return stage_file(source_path, destination_dir, fingerprint=fingerprint), True
    return str(source_path), False


def infer_reuse_requirements(manifest: dict[str, Any]) -> dict[str, Any]:
    reuse = manifest.get("reuse") or {}
    monogenic_enabled = manifest["analysis"]["monogenic"]["enabled"]
    prs_enabled = manifest["analysis"]["prs"]["enabled"]
    monogenic_reference_required = bool(monogenic_enabled and not (reuse.get("annotated_vcf") or reuse.get("called_vcf")))
    prs_reference_required = bool(prs_enabled and not (reuse.get("prs_dense_vcf") or reuse.get("imputed_prs_vcf")))
    return {
        "monogenic_reference_required": monogenic_reference_required,
        "prs_reference_required": prs_reference_required,
        "any_reference_required": monogenic_reference_required or prs_reference_required,
    }


def detect_upstream_sample_name(cram_path: str) -> str | None:
    completed = run_command(["samtools", "view", "-H", cram_path], check=True)
    for line in completed.stdout.splitlines():
        if not line.startswith("@RG"):
            continue
        match = re.search(r"\bSM:([^\t]+)", line)
        if match:
            return match.group(1)
    return None


def run_remote_decode_probe(cram_path: str, first_contig: str, cache_root: pathlib.Path, repo_root: pathlib.Path) -> dict[str, Any]:
    ensure_dir(cache_root)
    region = f"{first_contig}:1-1000"
    env = {
        "REF_CACHE": str(cache_root / "%2s" / "%2s" / "%s"),
        "REF_PATH": f"{cache_root}/%2s/%2s/%s:https://www.ebi.ac.uk/ena/cram/md5/%s",
    }
    try:
        completed = run_command(
            ["samtools", "view", cram_path, region],
            cwd=repo_root,
            env=env,
            check=False,
            timeout=120,
        )
    except Exception as exc:
        return {
            "attempted": True,
            "status": "failed",
            "region": region,
            "cache_dir": str(cache_root),
            "message": f"Remote MD5 decode probe crashed: {exc}",
        }

    if completed.returncode == 0:
        return {
            "attempted": True,
            "status": "ok",
            "region": region,
            "cache_dir": str(cache_root),
            "message": "samtools successfully decoded CRAM slices via REF_PATH/ENA MD5 resolution.",
        }
    return {
        "attempted": True,
        "status": "failed",
        "region": region,
        "cache_dir": str(cache_root),
        "message": "samtools could not decode CRAM slices via REF_PATH/ENA MD5 resolution.",
        "exit_code": completed.returncode,
        "stderr_tail": (completed.stderr or "").splitlines()[-10:],
    }


def locate_stage_manifest_for_import(target_path: pathlib.Path) -> pathlib.Path | None:
    for candidate_dir in [target_path.parent, *target_path.parents[:3]]:
        candidate = candidate_dir / "stage_manifest.json"
        if candidate.exists():
            return candidate
    return None


def prepare_import_artifact(
    label: str,
    raw_path: str | None,
    *,
    manifest_dir: pathlib.Path,
    staging_root: pathlib.Path,
    linux_cache_root: pathlib.Path,
    current_input_fingerprint_id: str,
) -> dict[str, Any] | None:
    if not raw_path:
        return None

    normalized = pathlib.Path(normalize_path(raw_path, base_dir=manifest_dir))
    if not normalized.exists():
        raise SystemExit(f"Configured reuse artifact not found: {normalized}")
    fingerprint = sample_file_fingerprint(normalized)

    active, staged = stage_maybe(normalized, staging_root / label, linux_cache_root)
    active_path = pathlib.Path(active)
    sidecars: dict[str, str] = {}
    if active_path.name.endswith(".vcf.gz"):
        index_path = active_path.parent / f"{active_path.name}.tbi"
        source_index = pathlib.Path(f"{normalized}.tbi")
        if staged and source_index.exists() and not index_path.exists():
            stage_maybe(source_index, staging_root / label, linux_cache_root)
        if index_path.exists():
            sidecars["index"] = str(index_path)
        elif source_index.exists() and not staged:
            sidecars["index"] = str(source_index)
        else:
            run_command(["tabix", "-f", "-p", "vcf", str(active_path)], check=True)
            sidecars["index"] = str(index_path)

    stage_manifest = locate_stage_manifest_for_import(normalized)
    provenance = "external_import"
    trusted = False
    if stage_manifest is not None:
        payload = load_json(stage_manifest)
        if payload.get("input_fingerprint_id") == current_input_fingerprint_id:
            provenance = "repo_import_verified"
            trusted = True
        else:
            provenance = "repo_import_unverified"

    return {
        "original": str(normalized),
        "active": str(active_path),
        "staged": staged,
        "fingerprint": fingerprint,
        "provenance": provenance,
        "trusted_provenance": trusted,
        "sidecars": sidecars,
        "stage_manifest": str(stage_manifest) if stage_manifest else None,
    }


def detect_required_executables() -> dict[str, str | None]:
    names = ["docker", "java", "nextflow", "python3", "bcftools", "samtools", "tabix"]
    return {name: detect_executable(name) for name in names}


def cmd_prepare_context(args: argparse.Namespace) -> int:
    validated = load_json(args.validated_manifest)
    manifest_path = pathlib.Path(validated["manifest_path"]).resolve()
    manifest = validated["manifest"]
    execution_mode = validated.get("execution_mode") or "both"

    linux_cache_root = pathlib.Path(normalize_path(manifest["runtime"]["linux_cache_root"], base_dir=manifest_path.parent))
    project_cache_root = ensure_dir(linux_cache_root / project_cache_id())
    staging_root = ensure_dir(project_cache_root / "staging")
    nextflow_root = ensure_dir(project_cache_root / "nextflow")
    tool_cache_root = ensure_dir(project_cache_root / "tool_cache")

    results_root = pathlib.Path(normalize_path(manifest["outputs"]["results_dir"], base_dir=REPO_ROOT))
    results_run_dir = ensure_dir(results_root / manifest["sample_id"] / manifest["outputs"]["run_id"])
    report_dir = ensure_dir(results_run_dir / "report")
    manifests_dir = ensure_dir(results_run_dir / "manifests")
    monogenic_results_dir = ensure_dir(results_run_dir / "monogenic")
    prs_results_dir = ensure_dir(results_run_dir / "prs")

    normalized_inputs = {
        "cram": pathlib.Path(normalize_path(manifest["inputs"]["cram"], base_dir=manifest_path.parent)),
        "crai": pathlib.Path(normalize_path(manifest["inputs"]["crai"], base_dir=manifest_path.parent)),
    }
    for label, path in normalized_inputs.items():
        if not path.exists():
            raise SystemExit(f"Missing {label.upper()}: {path}")

    cram_fingerprint = sample_file_fingerprint(normalized_inputs["cram"])
    crai_fingerprint = sample_file_fingerprint(normalized_inputs["crai"])
    input_fingerprint = {
        "id": hash_strings([cram_fingerprint["id"], crai_fingerprint["id"]]),
        "cram": cram_fingerprint,
        "crai": crai_fingerprint,
    }
    sample_cache_root = ensure_dir(project_cache_root / "samples" / input_fingerprint["id"])
    monogenic_cache_root = ensure_dir(sample_cache_root / "monogenic")
    prs_cache_root = ensure_dir(sample_cache_root / "prs")

    staged_cram, staged_cram_flag = stage_maybe(
        normalized_inputs["cram"],
        staging_root / "inputs" / input_fingerprint["id"],
        linux_cache_root,
    )
    staged_crai, staged_crai_flag = stage_maybe(
        normalized_inputs["crai"],
        staging_root / "inputs" / input_fingerprint["id"],
        linux_cache_root,
    )

    reuse_cfg = manifest.get("reuse") or {}
    prepared_reuse = {
        key: prepare_import_artifact(
            key,
            reuse_cfg.get(key),
            manifest_dir=manifest_path.parent,
            staging_root=staging_root / "imports" / input_fingerprint["id"],
            linux_cache_root=linux_cache_root,
            current_input_fingerprint_id=input_fingerprint["id"],
        )
        for key in ("called_vcf", "annotated_vcf", "prs_dense_vcf", "imputed_prs_vcf")
    }

    executable_checks = detect_required_executables()
    missing = [name for name, value in executable_checks.items() if value is None]
    if missing:
        raise SystemExit("Missing required executables: " + ", ".join(sorted(missing)))

    disk_usage = shutil.disk_usage(project_cache_root)
    free_gb = disk_usage.free / (1024 ** 3)
    min_free_gb = manifest["runtime"]["min_free_gb"]
    if free_gb < min_free_gb:
        raise SystemExit(
            f"Only {free_gb:.1f} GiB free on {project_cache_root}, but manifest requires at least {min_free_gb} GiB"
        )

    header_completed = run_command(["samtools", "view", "-H", staged_cram], check=True)
    header_text = header_completed.stdout
    sq_entries = parse_sq_entries(header_text)
    pg_entries = parse_program_entries(header_text)
    rg_entries = parse_read_group_entries(header_text)
    contig_summary = summarize_contigs(sq_entries)
    provenance_hints = extract_reference_hints(pg_entries)
    contig_names = [entry["name"] for entry in sq_entries if entry.get("name")]
    shards = plan_contig_shards(contig_names, contig_summary["contig_prefix"])

    remote_probe = {
        "attempted": False,
        "status": "not_attempted",
        "message": "Remote MD5 decode probe disabled in manifest.",
    }
    if manifest["reference"].get("remote_decode_probe") and contig_summary["first_contigs"]:
        remote_probe = run_remote_decode_probe(
            staged_cram,
            contig_summary["first_contigs"][0],
            project_cache_root / "references" / "remote_md5_cache",
            REPO_ROOT,
        )

    stage_requirements = infer_reuse_requirements(manifest)
    reference_cfg = manifest["reference"]
    bundle_id = reference_cfg["bundle_id"]
    bundle_manifest_path = resolve_bundle_manifest(reference_cfg)
    reference_stage: dict[str, Any]

    if bundle_manifest_path is None:
        status = "failed" if stage_requirements["any_reference_required"] else "skipped"
        reference_stage = stage_payload(
            "reference_bundle",
            status,
            required=stage_requirements["any_reference_required"],
            bundle_id=bundle_id,
            message=f"No bundle manifest found for reference bundle_id={bundle_id}",
            remote_decode_probe=remote_probe,
        )
    else:
        bundle = load_yaml(bundle_manifest_path) or {}
        for key in ("fasta", "fai", "dict"):
            if bundle.get(key):
                bundle[key] = normalize_path(bundle[key], base_dir=bundle_manifest_path.parent)
        if reference_cfg.get("autodownload"):
            bundle = maybe_download_bundle(bundle, project_cache_root / "references" / bundle_id)

        try:
            fasta, fai, seq_dict, missing_bundle_files = ensure_bundle_files(bundle)
        except Exception as exc:
            fasta = fai = seq_dict = None
            missing_bundle_files = ["prepare_failed"]
            prepare_error = str(exc)
        else:
            prepare_error = None

        if missing_bundle_files:
            status = "failed" if stage_requirements["any_reference_required"] else "skipped"
            reference_stage = stage_payload(
                "reference_bundle",
                status,
                required=stage_requirements["any_reference_required"],
                bundle_id=bundle_id,
                bundle_manifest=str(bundle_manifest_path),
                message=prepare_error or f"Reference bundle is missing required files: {', '.join(missing_bundle_files)}",
                remote_decode_probe=remote_probe,
            )
        else:
            reference_stage_root = ensure_dir(project_cache_root / "references" / bundle_id)
            staged_fasta, fasta_staged = stage_maybe(fasta, reference_stage_root, linux_cache_root)
            staged_fai, fai_staged = stage_maybe(fai, reference_stage_root, linux_cache_root)
            staged_dict, dict_staged = stage_maybe(seq_dict, reference_stage_root, linux_cache_root)
            comparison = compare_sequence_catalogs(sq_entries, parse_sequence_dict(staged_dict))
            expected_prefix = bundle.get("expected_contig_prefix", "")
            observed_prefix = contig_summary["contig_prefix"]
            prefix_matches = expected_prefix in {"", observed_prefix}
            status = "ready" if comparison["exact_for_cram"] else ("failed" if stage_requirements["any_reference_required"] else "skipped")
            reference_stage = stage_payload(
                "reference_bundle",
                status,
                required=stage_requirements["any_reference_required"],
                bundle_id=bundle_id,
                bundle_manifest=str(bundle_manifest_path),
                fasta=staged_fasta,
                fai=staged_fai,
                dict=staged_dict,
                staged={"fasta": fasta_staged, "fai": fai_staged, "dict": dict_staged},
                expected_contig_prefix=expected_prefix,
                observed_contig_prefix=observed_prefix,
                prefix_matches=prefix_matches,
                comparison=comparison,
                provenance_hints=provenance_hints,
                remote_decode_probe=remote_probe,
                message=None if comparison["exact_for_cram"] else "Local reference bundle did not validate against the CRAM header.",
            )

    context = {
        "context_version": 2,
        "repo_root": str(REPO_ROOT),
        "manifest_path": str(manifest_path),
        "manifest": manifest,
        "execution_mode": execution_mode,
        "sample_id": manifest["sample_id"],
        "input_fingerprint": input_fingerprint,
        "inputs": {
            "original": {key: str(value) for key, value in normalized_inputs.items()},
            "active": {"cram": staged_cram, "crai": staged_crai},
            "staged": {"cram": staged_cram_flag, "crai": staged_crai_flag},
        },
        "reuse": prepared_reuse,
        "runtime": {
            **manifest["runtime"],
            "linux_cache_root": str(linux_cache_root),
            "project_cache_root": str(project_cache_root),
            "nextflow_home": str(ensure_dir(nextflow_root / "home")),
            "nextflow_launcher": str(ensure_dir(nextflow_root / "launcher")),
            "nextflow_work": str(ensure_dir(nextflow_root / "work")),
            "tool_cache_root": str(tool_cache_root),
        },
        "results": {
            "results_root": str(results_root),
            "run_dir": str(results_run_dir),
            "report_dir": str(report_dir),
            "manifests_dir": str(manifests_dir),
            "monogenic_dir": str(monogenic_results_dir),
            "prs_dir": str(prs_results_dir),
        },
        "cache": {
            "sample_root": str(sample_cache_root),
            "staging_root": str(staging_root),
            "monogenic_root": str(monogenic_cache_root),
            "prs_root": str(prs_cache_root),
        },
        "host_checks": executable_checks,
        "disk": {
            "usage_target": str(project_cache_root),
            "free_gb": round(free_gb, 2),
            "min_free_gb": min_free_gb,
        },
        "stage_requirements": stage_requirements,
        "upstream_sample_name": detect_upstream_sample_name(staged_cram),
        "contig_prefix": contig_summary["contig_prefix"],
        "contigs": contig_names,
        "shards": {
            "monogenic": shards,
            "prs": shards,
        },
        "reference_probe": stage_payload(
            "reference_probe",
            "ok",
            summary={
                "contig_count": contig_summary["count"],
                "first_contigs": contig_summary["first_contigs"],
                "contig_prefix": contig_summary["contig_prefix"],
                "all_contigs_have_md5": contig_summary["all_contigs_have_md5"],
                "program_count": len(pg_entries),
                "read_group_count": len(rg_entries),
            },
            cram_header={"sq_entries": sq_entries, "pg_entries": pg_entries, "rg_entries": rg_entries},
            provenance_hints=provenance_hints,
            remote_decode_probe=remote_probe,
        ),
        "reference_bundle": reference_stage,
    }
    write_json(args.output, context)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare-context")
    prepare_parser.add_argument("--validated-manifest", required=True)
    prepare_parser.add_argument("--output", required=True)
    prepare_parser.set_defaults(func=cmd_prepare_context)

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
        print(f"Reference command crashed: {exc}", file=sys.stderr)
        raise
