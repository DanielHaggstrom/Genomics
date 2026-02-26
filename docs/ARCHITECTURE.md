# Project Architecture

## Goal

Run reproducible personal-genomics analysis on immutable WGS inputs:

1. Monogenic triage from ClinVar-enriched VEP annotations.
2. Polygenic scoring (PRS) from selected PGS IDs.

## Data Model

1. Inputs (immutable, outside repo):
   - `*.vcf.gz`, `*.tbi`, `*.cram`, `*.crai` under your local WGS data directory.
2. Config:
   - Tracked templates: `config/*.example*`, `.env.example`
   - Local/private copies: `.env`, `config/pipeline.env`, `config/secrets.env`, `config/paths.env`, PRS local config files
3. Derived outputs:
   - `results/monogenic/<SAMPLE_ID>/`
   - `results/prs/pgsc_calc/`
4. Runtime artifacts:
   - `work/`, `.nextflow/`, `.nextflow.log*`
5. Operational logs:
   - `logs/`

## Workflow Layers

1. Orchestration layer:
   - `workflow/run_pipeline.sh`
2. Monogenic layer:
   - `workflow/monogenic/extract_priority_variants.py`
   - `workflow/monogenic/build_monogenic_shortlists.py`
   - `workflow/monogenic/triage_clinvar_enriched.py`
   - `workflow/monogenic/dedup_shortlist_by_variant.py`
3. PRS layer:
   - `workflow/prs/prepare_prs_samplesheet.sh`
   - `workflow/prs/run_pgsc_calc.sh`
   - Nextflow pipeline: `pgscatalog/pgsc_calc`

## Runtime Policy

1. Container-first PRS execution: `docker`, then `apptainer`/`singularity`.
2. `conda` fallback only when containers are unavailable.
3. Long runs in `tmux`.

## Publish-Safe Policy

1. Never commit local/private env files.
2. Commit only example templates with fake values.
3. Never commit genomic source data or generated run outputs.