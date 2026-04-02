#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import urllib.request
from collections import defaultdict
from typing import Any

from workflow.core.runtime import (
    artifacts_exist,
    build_docker_run_command,
    complete_artifact,
    ensure_dir,
    ensure_docker_ready,
    init_shard_progress,
    load_json,
    plan_worker_pool,
    record_shard_completion,
    remove_if_exists,
    run_command,
    run_tasks_in_parallel,
    stage_payload,
    write_json,
)


STAGE_VERSION = 3
GATK_IMAGE = "broadinstitute/gatk:4.6.1.0"
VEP_IMAGE = "ensemblorg/ensembl-vep:release_115.0"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_TBI_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"

PROTEIN_ALTERING = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
}

LOF_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
}

CLIN_SIG_TERMS = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
}

HIGH_REVIEW_MARKERS = (
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "criteria_provided&_multiple_submitters&_no_conflicts",
)

SINGLE_REVIEW_MARKER = "criteria_provided&_single_submitter"


def parse_info(info: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = True
    return out


def carried_alt_alleles(gt: str, alt: str) -> set[str]:
    alts = alt.split(",")
    carried: set[str] = set()
    for token in gt.replace("|", "/").split("/"):
        if token in {"", ".", "0"}:
            continue
        try:
            idx = int(token) - 1
        except ValueError:
            continue
        if 0 <= idx < len(alts):
            carried.add(alts[idx])
    return carried


def has_any_term(consequence: str, terms: set[str]) -> bool:
    return any(term in {block.strip() for block in consequence.split("&") if block.strip()} for term in terms)


def is_plp(sig: str) -> bool:
    lowered = sig.lower()
    return ("pathogenic" in lowered) or ("likely_pathogenic" in lowered)


def has_conflict(sig: str, rev: str) -> bool:
    return "conflicting_classifications" in f"{sig} {rev}".lower()


def high_review(rev: str) -> bool:
    lowered = rev.lower()
    normalized = lowered.replace("&", ",")
    return any(marker in lowered or marker in normalized for marker in HIGH_REVIEW_MARKERS)


def single_submitter_review(rev: str) -> bool:
    return SINGLE_REVIEW_MARKER in rev.lower()


def review_tier(rev: str) -> str:
    if not rev:
        return "none"
    if high_review(rev):
        return "high"
    if single_submitter_review(rev):
        return "single_submitter"
    return "other"


def parse_float(raw: str) -> float | None:
    try:
        return float(raw)
    except Exception:
        return None


def load_existing_stage(stage_manifest_path: pathlib.Path, cache_key: str, required_artifacts: list[str]) -> dict[str, Any] | None:
    if not stage_manifest_path.exists():
        return None
    payload = load_json(stage_manifest_path)
    if payload.get("cache_key") != cache_key or payload.get("status") != "ok":
        return None
    artifacts = payload.get("artifacts") or {}
    if not artifacts_exist([artifacts.get(key) for key in required_artifacts]):
        return None
    payload["reused"] = True
    return payload


def write_tsv(path: pathlib.Path, rows: list[dict[str, str]], header: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def extract_priority_variants(input_vcf: pathlib.Path, output_tsv: pathlib.Path, summary_txt: pathlib.Path) -> dict[str, str]:
    csq_fields: list[str] = []
    total_records = 0
    priority_rows: list[dict[str, str]] = []
    reason_counts: defaultdict[str, int] = defaultdict(int)

    with gzip.open(input_vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                marker = "Format: "
                idx = line.find(marker)
                if idx != -1:
                    csq_fields = line[idx + len(marker):].rstrip('">').split("|")
                continue
            if line.startswith("#"):
                continue

            total_records += 1
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample = cols[:10]
            info_blob = parse_info(info)
            csq_blob = info_blob.get("CSQ", "")
            if not csq_blob or not csq_fields:
                continue

            fmt_keys = fmt.split(":")
            sample_values = sample.split(":")
            sample_map = dict(zip(fmt_keys, sample_values))
            gt = sample_map.get("GT", ".")
            carried_alts = carried_alt_alleles(gt, alt) if "," in alt else set()

            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                ann_map = dict(zip(csq_fields, parts))

                if carried_alts:
                    allele = ann_map.get("Allele", "")
                    if allele and allele not in carried_alts:
                        continue

                consequence = ann_map.get("Consequence", "")
                clin_sig = ann_map.get("CLIN_SIG", "")
                cons_terms = {item.strip() for item in consequence.split("&") if item.strip()}
                clin_terms = {item.strip().lower() for item in clin_sig.replace(" ", "_").split("&") if item.strip()}
                lof_hit = any(term in LOF_TERMS for term in cons_terms)
                clin_hit = any(term in CLIN_SIG_TERMS for term in clin_terms)
                if not (lof_hit or clin_hit):
                    continue

                reasons = []
                if clin_hit:
                    reasons.append("clinvar_pathogenic_or_likely_pathogenic")
                    reason_counts["clinvar_pathogenic_or_likely_pathogenic"] += 1
                if lof_hit:
                    reasons.append("predicted_loss_of_function")
                    reason_counts["predicted_loss_of_function"] += 1

                priority_rows.append(
                    {
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": vid,
                        "REF": ref,
                        "ALT": alt,
                        "FILTER": flt,
                        "GT": gt,
                        "DP": sample_map.get("DP", "."),
                        "GQ": sample_map.get("GQ", "."),
                        "AD": sample_map.get("AD", "."),
                        "GENE_SYMBOL": ann_map.get("SYMBOL", ""),
                        "GENE_ID": ann_map.get("Gene", ""),
                        "CONSEQUENCE": consequence,
                        "IMPACT": ann_map.get("IMPACT", ""),
                        "CLIN_SIG": clin_sig,
                        "EXISTING_VARIATION": ann_map.get("Existing_variation", ""),
                        "HGVSc": ann_map.get("HGVSc", ""),
                        "HGVSp": ann_map.get("HGVSp", ""),
                        "CANONICAL": ann_map.get("CANONICAL", ""),
                        "REASONS": ",".join(reasons),
                    }
                )

    seen = set()
    deduped: list[dict[str, str]] = []
    for row in priority_rows:
        key = tuple(row.values())
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)

    header = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "FILTER",
        "GT",
        "DP",
        "GQ",
        "AD",
        "GENE_SYMBOL",
        "GENE_ID",
        "CONSEQUENCE",
        "IMPACT",
        "CLIN_SIG",
        "EXISTING_VARIATION",
        "HGVSc",
        "HGVSp",
        "CANONICAL",
        "REASONS",
    ]
    write_tsv(output_tsv, deduped, header)

    summary = {
        "input_vcf": str(input_vcf),
        "total_variant_records_seen": str(total_records),
        "priority_rows_before_dedup": str(len(priority_rows)),
        "priority_rows_after_dedup": str(len(deduped)),
    }
    for key in sorted(reason_counts):
        summary[f"count_{key}"] = str(reason_counts[key])
    summary_txt.write_text("".join(f"{key}\t{value}\n" for key, value in summary.items()), encoding="utf-8")
    return summary


def dedup_variant_gene_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        key = (
            row.get("CHROM", ""),
            row.get("POS", ""),
            row.get("REF", ""),
            row.get("ALT", ""),
            row.get("GENE_SYMBOL", ""),
            row.get("GENE_ID", ""),
        )
        grouped.setdefault(key, []).append(row)
    picked: list[dict[str, str]] = []
    for key in sorted(grouped):
        bucket = grouped[key]
        canonical = [item for item in bucket if item.get("CANONICAL") == "YES"]
        picked.append(canonical[0] if canonical else bucket[0])
    return picked


def build_shortlists(input_tsv: pathlib.Path, outdir: pathlib.Path) -> None:
    with open(input_tsv, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []

    pass_rows = [row for row in rows if row.get("FILTER") == "PASS"]
    clinvar_plp = [
        row
        for row in pass_rows
        if "clinvar_pathogenic_or_likely_pathogenic" in (row.get("REASONS") or "")
        and has_any_term(row.get("CONSEQUENCE", ""), PROTEIN_ALTERING)
    ]
    lof_rows = [
        row
        for row in pass_rows
        if has_any_term(row.get("CONSEQUENCE", ""), LOF_TERMS)
        and has_any_term(row.get("CONSEQUENCE", ""), PROTEIN_ALTERING)
    ]
    write_tsv(outdir / "clinvar_plp_protein_altering.tsv", clinvar_plp, header)
    write_tsv(outdir / "clinvar_plp_unique_variants.tsv", dedup_variant_gene_rows(clinvar_plp), header)
    write_tsv(outdir / "lof_protein_altering.tsv", lof_rows, header)
    write_tsv(outdir / "lof_unique_variants.tsv", dedup_variant_gene_rows(lof_rows), header)


def triage_clinvar(input_vcf: pathlib.Path, output_tsv: pathlib.Path, *, mode: str, max_af: float) -> None:
    csq_fields: list[str] = []
    rows: list[dict[str, str]] = []

    with gzip.open(input_vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format: ", 1)[1].rstrip('">').split("|")
                continue
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample = line.split("\t")[:10]
            if flt != "PASS":
                continue

            info_blob = parse_info(info)
            csq_blob = info_blob.get("CSQ", "")
            if not csq_blob:
                continue

            fmt_keys = fmt.split(":")
            sample_values = sample.split(":")
            sample_map = dict(zip(fmt_keys, sample_values))
            carried_alts = carried_alt_alleles(sample_map.get("GT", "."), alt) if "," in alt else set()

            for ann in csq_blob.split(","):
                parts = ann.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                ann_map = dict(zip(csq_fields, parts))

                if carried_alts:
                    allele = ann_map.get("Allele", "")
                    if allele and allele not in carried_alts:
                        continue

                consequence = ann_map.get("Consequence", "")
                if not has_any_term(consequence, PROTEIN_ALTERING):
                    continue

                clnsig = ann_map.get("ClinVar_CLNSIG", "")
                clnrev = ann_map.get("ClinVar_CLNREVSTAT", "")
                if not is_plp(clnsig):
                    continue

                rev_tier = review_tier(clnrev)
                conflict = has_conflict(clnsig, clnrev)
                impact = ann_map.get("IMPACT", "")
                canonical = ann_map.get("CANONICAL", "")
                max_af_str = ann_map.get("MAX_AF", "")
                max_af_value = parse_float(max_af_str)
                select = False
                reason = ""

                if mode == "strict":
                    if (not conflict) and rev_tier == "high" and max_af_value is not None and max_af_value <= max_af:
                        select = True
                        reason = "high_review_and_rare"
                elif mode == "moderate":
                    if max_af_value is not None and max_af_value <= max_af and (rev_tier != "none" or impact == "HIGH" or has_any_term(consequence, LOF_TERMS)):
                        select = True
                        reason = "rare_plp_conflicting" if conflict else "rare_plp_reviewed"
                    elif max_af_value is None and (impact == "HIGH" or has_any_term(consequence, LOF_TERMS)) and rev_tier != "none":
                        select = True
                        reason = "high_impact_missing_af"
                else:
                    if (not conflict) and max_af_value is not None and max_af_value <= max_af:
                        select = True
                        reason = "nonconflicting_plp_low_frequency"
                    elif (not conflict) and max_af_value is None and (impact == "HIGH" or has_any_term(consequence, LOF_TERMS) or canonical == "YES"):
                        select = True
                        reason = "nonconflicting_plp_severe_missing_af"
                    elif conflict and max_af_value is not None and max_af_value <= 0.01 and (impact == "HIGH" or has_any_term(consequence, LOF_TERMS)):
                        select = True
                        reason = "conflicting_but_rare_high_impact"

                if not select:
                    continue

                rows.append(
                    {
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": vid,
                        "REF": ref,
                        "ALT": alt,
                        "GT": sample_map.get("GT", "."),
                        "DP": sample_map.get("DP", "."),
                        "GQ": sample_map.get("GQ", "."),
                        "GENE_SYMBOL": ann_map.get("SYMBOL", ""),
                        "GENE_ID": ann_map.get("Gene", ""),
                        "CONSEQUENCE": consequence,
                        "IMPACT": impact,
                        "CANONICAL": canonical,
                        "HGVSc": ann_map.get("HGVSc", ""),
                        "HGVSp": ann_map.get("HGVSp", ""),
                        "EXISTING_VARIATION": ann_map.get("Existing_variation", ""),
                        "CLINVAR_CLNSIG": clnsig,
                        "CLINVAR_CLNREVSTAT": clnrev,
                        "CLINVAR_CLNDN": ann_map.get("ClinVar_CLNDN", ""),
                        "MAX_AF": max_af_str,
                        "MAX_AF_POPS": ann_map.get("MAX_AF_POPS", ""),
                        "CLINVAR_CONFLICT_FLAG": "YES" if conflict else "NO",
                        "REVIEW_TIER": rev_tier,
                        "SELECTION_MODE": mode,
                        "SELECTION_REASON": reason,
                    }
                )

    seen = set()
    deduped: list[dict[str, str]] = []
    for row in rows:
        key = tuple(row.values())
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)

    header = list(deduped[0].keys()) if deduped else [
        "CHROM", "POS", "ID", "REF", "ALT", "GT", "DP", "GQ", "GENE_SYMBOL", "GENE_ID",
        "CONSEQUENCE", "IMPACT", "CANONICAL", "HGVSc", "HGVSp", "EXISTING_VARIATION",
        "CLINVAR_CLNSIG", "CLINVAR_CLNREVSTAT", "CLINVAR_CLNDN", "MAX_AF", "MAX_AF_POPS",
        "CLINVAR_CONFLICT_FLAG", "REVIEW_TIER", "SELECTION_MODE", "SELECTION_REASON",
    ]
    write_tsv(output_tsv, deduped, header)


def dedup_tsv(input_tsv: pathlib.Path, output_tsv: pathlib.Path) -> None:
    with open(input_tsv, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []
    write_tsv(output_tsv, dedup_variant_gene_rows(rows), header)


def build_gatk_command(
    *,
    docker_bin: str,
    cram_path: pathlib.Path,
    reference_fasta: pathlib.Path,
    output_dir: pathlib.Path,
    image: str,
    inner_args: list[str],
) -> list[str]:
    return build_docker_run_command(
        image=image,
        mounts={
            cram_path.parent: "/input",
            reference_fasta.parent: "/ref",
            output_dir: "/out",
        },
        inner_args=inner_args,
        docker_bin=docker_bin,
    )


def call_monogenic_shard(
    *,
    docker_bin: str,
    stage_root: pathlib.Path,
    cram_path: pathlib.Path,
    reference_fasta: pathlib.Path,
    shard: dict[str, Any],
    sample_id: str,
    pairhmm_threads: int,
    java_memory: str,
) -> dict[str, Any]:
    shard_dir = ensure_dir(stage_root / "shards" / shard["name"])
    log_dir = ensure_dir(stage_root / "logs")
    gvcf = shard_dir / f"{sample_id}.{shard['name']}.g.vcf.gz"
    raw_vcf = shard_dir / f"{sample_id}.{shard['name']}.raw.vcf.gz"
    interval_args = [item for interval in shard["intervals"] for item in ("-L", interval)]
    gvcf_reused = complete_artifact(gvcf)
    raw_vcf_reused = complete_artifact(raw_vcf)

    if not gvcf_reused:
        remove_if_exists(gvcf, f"{gvcf}.tbi")
        run_command(
            build_gatk_command(
                docker_bin=docker_bin,
                cram_path=cram_path,
                reference_fasta=reference_fasta,
                output_dir=stage_root,
                image=GATK_IMAGE,
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "HaplotypeCaller",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-I", f"/input/{cram_path.name}",
                    "-O", f"/out/shards/{shard['name']}/{gvcf.name}",
                    "-ERC", "GVCF",
                    "--native-pair-hmm-threads", str(pairhmm_threads),
                    *interval_args,
                ],
            ),
            check=True,
            log_path=log_dir / f"{shard['name']}.haplotypecaller.log",
        )

    if not raw_vcf_reused:
        remove_if_exists(raw_vcf, f"{raw_vcf}.tbi")
        run_command(
            build_gatk_command(
                docker_bin=docker_bin,
                cram_path=cram_path,
                reference_fasta=reference_fasta,
                output_dir=stage_root,
                image=GATK_IMAGE,
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "GenotypeGVCFs",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/shards/{shard['name']}/{gvcf.name}",
                    "-O", f"/out/shards/{shard['name']}/{raw_vcf.name}",
                    *interval_args,
                ],
            ),
            check=True,
            log_path=log_dir / f"{shard['name']}.genotypegvcfs.log",
        )

    return {
        "shard_name": shard["name"],
        "raw_vcf": str(raw_vcf),
        "gvcf": str(gvcf),
        "reused": raw_vcf_reused,
        "reused_gvcf": gvcf_reused,
        "reused_raw_vcf": raw_vcf_reused,
    }


def concat_and_filter_monogenic(stage_root: pathlib.Path, reference_fasta: pathlib.Path, sample_id: str, shard_vcfs: list[str], docker_bin: str) -> tuple[str, str]:
    raw_vcf = stage_root / f"{sample_id}.raw.vcf.gz"
    filtered_vcf = stage_root / f"{sample_id}.hard-filtered.vcf.gz"
    tmp_snp = stage_root / f"{sample_id}.snps.filtered.vcf.gz"
    tmp_indel = stage_root / f"{sample_id}.indels.filtered.vcf.gz"
    java_memory = "4g"
    log_dir = ensure_dir(stage_root / "logs")

    if not complete_artifact(raw_vcf):
        remove_if_exists(raw_vcf, f"{raw_vcf}.tbi")
        run_command(
            ["bcftools", "concat", "-a", "-Oz", "-o", str(raw_vcf), *shard_vcfs],
            check=True,
            log_path=log_dir / "raw.concat.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(raw_vcf)],
            check=True,
            log_path=log_dir / "raw.concat.index.log",
        )

    if not complete_artifact(tmp_snp):
        remove_if_exists(stage_root / f"{sample_id}.snps.vcf.gz", stage_root / f"{sample_id}.snps.vcf.gz.tbi", tmp_snp, f"{tmp_snp}.tbi")
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={reference_fasta.parent: "/ref", stage_root: "/out"},
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "SelectVariants",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/{raw_vcf.name}",
                    "--select-type-to-include", "SNP",
                    "-O", f"/out/{sample_id}.snps.vcf.gz",
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / "select_snps.log",
        )
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={reference_fasta.parent: "/ref", stage_root: "/out"},
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "VariantFiltration",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/{sample_id}.snps.vcf.gz",
                    "-O", f"/out/{tmp_snp.name}",
                    "--filter-name", "SNP_HARD_FILTER",
                    "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0",
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / "filter_snps.log",
        )

    if not complete_artifact(tmp_indel):
        remove_if_exists(stage_root / f"{sample_id}.indels.vcf.gz", stage_root / f"{sample_id}.indels.vcf.gz.tbi", tmp_indel, f"{tmp_indel}.tbi")
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={reference_fasta.parent: "/ref", stage_root: "/out"},
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "SelectVariants",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/{raw_vcf.name}",
                    "--select-type-to-include", "INDEL",
                    "-O", f"/out/{sample_id}.indels.vcf.gz",
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / "select_indels.log",
        )
        run_command(
            build_docker_run_command(
                image=GATK_IMAGE,
                mounts={reference_fasta.parent: "/ref", stage_root: "/out"},
                inner_args=[
                    "gatk",
                    "--java-options",
                    f"-Xmx{java_memory}",
                    "VariantFiltration",
                    "-R", f"/ref/{reference_fasta.name}",
                    "-V", f"/out/{sample_id}.indels.vcf.gz",
                    "-O", f"/out/{tmp_indel.name}",
                    "--filter-name", "INDEL_HARD_FILTER",
                    "--filter-expression", "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0",
                ],
                docker_bin=docker_bin,
            ),
            check=True,
            log_path=log_dir / "filter_indels.log",
        )

    if not complete_artifact(filtered_vcf):
        remove_if_exists(filtered_vcf, f"{filtered_vcf}.tbi")
        run_command(
            ["bcftools", "concat", "-a", "-Oz", "-o", str(filtered_vcf), str(tmp_snp), str(tmp_indel)],
            check=True,
            log_path=log_dir / "filtered.concat.log",
        )
        run_command(
            ["tabix", "-f", "-p", "vcf", str(filtered_vcf)],
            check=True,
            log_path=log_dir / "filtered.concat.index.log",
        )

    return str(raw_vcf), str(filtered_vcf)


def ensure_vep_resources(tool_cache_root: pathlib.Path) -> pathlib.Path:
    vep_root = ensure_dir(tool_cache_root / "vep")
    cache_root = ensure_dir(vep_root / "cache")
    cache_info = cache_root / "homo_sapiens" / "115_GRCh38" / "info.txt"
    if not cache_info.exists():
        tarball = ensure_dir(cache_root / "tmp") / "homo_sapiens_vep_115_GRCh38.tar.gz"
        urllib.request.urlretrieve(
            "https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz",
            tarball,
        )
        run_command(["tar", "-xzf", str(tarball), "-C", str(cache_root)], check=True)
        if not cache_info.exists():
            raise RuntimeError(f"VEP cache extraction did not produce {cache_info}")
        tarball.unlink(missing_ok=True)

    custom_root = ensure_dir(vep_root / "custom")
    clinvar_vcf = custom_root / "clinvar.vcf.gz"
    clinvar_tbi = custom_root / "clinvar.vcf.gz.tbi"
    if not clinvar_vcf.exists():
        urllib.request.urlretrieve(CLINVAR_URL, clinvar_vcf)
    if not clinvar_tbi.exists():
        urllib.request.urlretrieve(CLINVAR_TBI_URL, clinvar_tbi)
    return vep_root


def annotate_with_vep(
    *,
    docker_bin: str,
    input_vcf: pathlib.Path,
    output_vcf: pathlib.Path,
    vep_data_root: pathlib.Path,
    forks: int,
) -> str:
    if not complete_artifact(input_vcf):
        run_command(["tabix", "-f", "-p", "vcf", str(input_vcf)], check=True)
    if complete_artifact(output_vcf):
        return str(output_vcf)

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    warning_path = output_vcf.parent / f"{output_vcf.name}_warnings.txt"
    log_path = output_vcf.parent / "vep.annotate.log"
    work_root = output_vcf.parent
    run_command(
        build_docker_run_command(
            image=VEP_IMAGE,
            mounts={
                work_root: "/work",
                vep_data_root: "/vep_data",
            },
            inner_args=[
                "vep",
                "--cache",
                "--offline",
                "--assembly", "GRCh38",
                "--species", "homo_sapiens",
                "--dir_cache", "/vep_data/cache",
                "--input_file", f"/work/{input_vcf.name}",
                "--output_file", f"/work/{output_vcf.name}",
                "--warning_file", f"/work/{warning_path.name}",
                "--format", "vcf",
                "--vcf",
                "--compress_output", "bgzip",
                "--force_overwrite",
                "--fork", str(forks),
                "--no_stats",
                "--symbol",
                "--canonical",
                "--check_existing",
                "--max_af",
                "--custom", "file=/vep_data/custom/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN",
            ],
            docker_bin=docker_bin,
        ),
        check=True,
        log_path=log_path,
    )
    run_command(["tabix", "-f", "-p", "vcf", str(output_vcf)], check=True, log_path=output_vcf.parent / "vep.index.log")
    return str(output_vcf)


def build_import_summary(reuse_entry: dict[str, Any] | None) -> dict[str, Any] | None:
    if not reuse_entry:
        return None
    return {
        "source": reuse_entry["original"],
        "active": reuse_entry["active"],
        "provenance": reuse_entry["provenance"],
        "trusted_provenance": reuse_entry["trusted_provenance"],
        "stage_manifest": reuse_entry.get("stage_manifest"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    context = load_json(args.context)
    cfg = context["manifest"]["analysis"]["monogenic"]
    stage_root = ensure_dir(pathlib.Path(context["cache"]["monogenic_root"]))
    stage_manifest_path = stage_root / "stage_manifest.json"
    cache_key = ":".join(
        [
            f"v{STAGE_VERSION}",
            context["input_fingerprint"]["id"],
            context["reference_bundle"].get("bundle_id", "no-reference"),
            cfg["mode"],
            str(cfg["max_af"]),
            str(cfg["vep_forks"]),
            str(bool(context["reuse"].get("called_vcf"))),
            str(bool(context["reuse"].get("annotated_vcf"))),
        ]
    )

    try:
        worker_plan: dict[str, Any] | None = None
        progress: dict[str, Any] | None = None

        def write_stage_status(status: str, *, phase: str | None = None, message: str | None = None, artifacts: dict[str, Any] | None = None) -> None:
            payload = stage_payload(
                "monogenic",
                status,
                cache_key=cache_key,
                input_fingerprint_id=context["input_fingerprint"]["id"],
                source_mode=source_mode or None,
                import_details=import_details,
                worker_plan=worker_plan,
                progress=progress,
                phase=phase,
                message=message,
                artifacts=artifacts or {"log_dir": str(stage_root / "logs")},
            )
            write_json(stage_manifest_path, {key: value for key, value in payload.items() if value is not None})

        if not cfg["enabled"]:
            payload = stage_payload("monogenic", "skipped", message="Monogenic stage disabled.")
            write_json(args.output, payload)
            return

        reused = load_existing_stage(
            stage_manifest_path,
            cache_key,
            [
                "annotated_vcf",
                "priority_variants_tsv",
                "clinvar_unique_tsv",
                "lof_unique_tsv",
                "review_unique_tsv",
            ],
        )
        if reused:
            write_json(args.output, reused)
            return

        docker_bin = context["host_checks"]["docker"] or "docker"
        ensure_docker_ready(docker_bin)

        called_vcf = context["reuse"].get("called_vcf")
        annotated_vcf_import = context["reuse"].get("annotated_vcf")
        reference_stage = context["reference_bundle"]
        source_mode = ""
        import_details = None

        if annotated_vcf_import:
            annotated_vcf = pathlib.Path(annotated_vcf_import["active"])
            source_mode = "import_annotated_vcf"
            import_details = build_import_summary(annotated_vcf_import)
            write_stage_status("running", phase="using_imported_annotated_vcf")
        else:
            if called_vcf:
                filtered_vcf = pathlib.Path(called_vcf["active"])
                source_mode = "import_called_vcf"
                import_details = build_import_summary(called_vcf)
                write_stage_status("running", phase="annotating_imported_called_vcf")
            else:
                if reference_stage.get("status") != "ready":
                    payload = stage_payload(
                        "monogenic",
                        "failed",
                        cache_key=cache_key,
                        input_fingerprint_id=context["input_fingerprint"]["id"],
                        message="CRAM-based monogenic calling requires a validated local reference bundle.",
                    )
                    write_json(stage_manifest_path, payload)
                    write_json(args.output, payload)
                    return

                sample_id = context["sample_id"]
                cram_path = pathlib.Path(context["inputs"]["active"]["cram"])
                reference_fasta = pathlib.Path(reference_stage["fasta"])
                shards = context["shards"]["monogenic"]
                worker_plan = plan_worker_pool(
                    stage_name="monogenic",
                    shard_count=len(shards),
                    max_cpus=int(context["manifest"]["runtime"]["max_cpus"]),
                    max_memory=context["manifest"]["runtime"]["max_memory"],
                    requested_workers=cfg.get("workers"),
                    requested_pairhmm_threads=cfg.get("pairhmm_threads"),
                    default_worker_cap=4,
                    default_pairhmm_cap=2,
                    min_memory_per_worker_gb=3.0,
                    free_gb=context["disk"].get("free_gb"),
                    min_free_gb=context["disk"].get("min_free_gb"),
                )
                progress = init_shard_progress(phase="calling_shards", total_shards=len(shards))
                write_stage_status("running", phase="calling_shards")

                tasks = [
                    (
                        lambda shard=shard: call_monogenic_shard(
                            docker_bin=docker_bin,
                            stage_root=stage_root,
                            cram_path=cram_path,
                            reference_fasta=reference_fasta,
                            shard=shard,
                            sample_id=sample_id,
                            pairhmm_threads=worker_plan["pairhmm_threads"],
                            java_memory=worker_plan["java_memory"],
                        )
                    )
                    for shard in shards
                ]
                shard_names = [shard["name"] for shard in shards]

                def on_shard_complete(index: int, result: Any) -> None:
                    assert progress is not None
                    shard_name = result.get("shard_name") or shard_names[index]
                    record_shard_completion(progress, shard_name=shard_name, reused=bool(result.get("reused")))
                    write_stage_status("running", phase="calling_shards")

                shard_results = run_tasks_in_parallel(tasks, max_workers=worker_plan["workers"], on_complete=on_shard_complete)
                shard_vcfs = [str(item["raw_vcf"]) for item in shard_results]
                progress = {
                    **(progress or {}),
                    "phase": "concat_and_filter",
                    "raw_shard_vcf_count": len(shard_vcfs),
                    "reused_raw_shard_vcfs": sum(1 for item in shard_results if item.get("reused")),
                }
                write_stage_status("running", phase="concat_and_filter")
                _, filtered_vcf_path = concat_and_filter_monogenic(stage_root, reference_fasta, sample_id, shard_vcfs, docker_bin)
                filtered_vcf = pathlib.Path(filtered_vcf_path)
                source_mode = "called_from_cram"
            if not called_vcf:
                import_details = None

            annotated_vcf = stage_root / f"{context['sample_id']}.vep.clinvar_enriched.vcf.gz"
            vep_root = ensure_vep_resources(pathlib.Path(context["runtime"]["tool_cache_root"]))
            progress = {
                **(progress or {}),
                "phase": "annotating_with_vep",
            }
            write_stage_status("running", phase="annotating_with_vep")
            annotate_with_vep(
                docker_bin=docker_bin,
                input_vcf=filtered_vcf,
                output_vcf=annotated_vcf,
                vep_data_root=vep_root,
                forks=int(cfg["vep_forks"]),
            )

        priority_tsv = stage_root / "priority_variants.tsv"
        priority_summary = stage_root / "priority_variants.summary.txt"
        progress = {
            **(progress or {}),
            "phase": "triaging_variants",
        }
        write_stage_status("running", phase="triaging_variants")
        summary = extract_priority_variants(annotated_vcf, priority_tsv, priority_summary)
        build_shortlists(priority_tsv, stage_root)
        review_queue = stage_root / "rare_clinvar_review_queue.tsv"
        review_unique = stage_root / "rare_clinvar_review_unique_variants.tsv"
        triage_clinvar(annotated_vcf, review_queue, mode=cfg["mode"], max_af=float(cfg["max_af"]))
        dedup_tsv(review_queue, review_unique)

        payload = stage_payload(
            "monogenic",
            "ok",
            cache_key=cache_key,
            input_fingerprint_id=context["input_fingerprint"]["id"],
            reused=False,
            source_mode=source_mode,
            import_details=import_details,
            worker_plan=worker_plan if not called_vcf and not annotated_vcf_import else None,
            progress=progress,
            summary=summary,
            artifacts={
                "annotated_vcf": str(annotated_vcf),
                "priority_variants_tsv": str(priority_tsv),
                "priority_summary": str(priority_summary),
                "clinvar_unique_tsv": str(stage_root / "clinvar_plp_unique_variants.tsv"),
                "lof_unique_tsv": str(stage_root / "lof_unique_variants.tsv"),
                "review_unique_tsv": str(review_unique),
            },
        )
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)
    except Exception as exc:
        payload = stage_payload(
            "monogenic",
            "failed",
            cache_key=cache_key,
            input_fingerprint_id=context["input_fingerprint"]["id"],
            message=str(exc),
            worker_plan=worker_plan,
            progress=progress,
            artifacts={"log_dir": str(stage_root / "logs")},
        )
        write_json(stage_manifest_path, payload)
        write_json(args.output, payload)


if __name__ == "__main__":
    main()
