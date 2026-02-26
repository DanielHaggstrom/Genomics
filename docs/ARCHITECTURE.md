# Project Architecture

## Goal

Run reproducible personal-genomics analysis on immutable WGS inputs:

1. Monogenic triage from ClinVar-enriched VEP annotations.
2. Polygenic scoring (PRS) from selected PGS IDs against the same sample.

## Data Model

1. Inputs (immutable, outside repo):
   - `*.vcf.gz`, `*.tbi`, `*.cram`, `*.crai` in `C:/Users/Daniel/Salud/WGS`.
2. Config (local + versioned templates):
   - `config/paths.env` (private)
   - `config/prs_samplesheet.csv` (private)
   - `config/prs_pgs_ids.txt` (private)
   - `config/prs_update_sex.tsv` (private; strict 3-column format)
3. Derived outputs:
   - `results/vep/*`
   - `results/monogenic/*`
   - `results/prs/*` (when PRS succeeds)
4. Audit logs:
   - `logs/input_checksums_pre.txt`
   - `logs/input_checksums_post.txt`
   - `logs/input_checksums_recheck_2026-02-08.txt`

## Workflow Layers

1. Annotation layer:
   - Ensembl VEP (offline cache, GRCh38)
   - Output: ClinVar-enriched annotated VCF
2. Monogenic triage layer:
   - `workflow/monogenic/triage_clinvar_enriched.py`
   - `workflow/monogenic/dedup_shortlist_by_variant.py`
3. PRS orchestration layer:
   - `workflow/prs/prepare_prs_samplesheet.sh`
   - `workflow/prs/run_pgsc_calc.sh`
   - Nextflow pipeline: `pgscatalog/pgsc_calc`

## Runtime Strategy

1. Container-first:
   - `docker` -> `apptainer`/`singularity`
2. `conda` fallback only.
3. Long-running PRS execution in `tmux` (not shell job-control background loops).

## Operational Boundaries

1. This repo never edits source genomic files.
2. Cleanup targets only transient runtime artifacts (`work/`, `.nextflow/`, failed run snapshots).
3. Reference build for this project is GRCh38 end-to-end.
