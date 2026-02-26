# Genomics WGS Risk Analysis Pipeline

Local, privacy-preserving workflow for personal whole-genome sequencing (WGS) analysis.

Primary objective:
- monogenic variant triage with ClinVar evidence and transcript-aware filtering.

Planned next objective:
- polygenic risk score (PRS) computation using PGS Catalog tooling.

## Strategy Update

This project now uses a strict runtime policy for PRS:

1. Prefer containerized execution (`docker`, then `apptainer`/`singularity`).
2. Use `conda` only as fallback when containers are unavailable.
3. Run long jobs in `tmux` (not fragile shell job control / ad-hoc `nohup`).

See `docs/PRS_EXECUTION_RUNBOOK.md` for the operational playbook.
See `docs/ARCHITECTURE.md` for the high-level system view.
See `docs/TARGETED_GENOTYPING_PLAN.md` for the CRAM -> PRS-targeted VCF strategy.

## Data Safety Policy

1. Source data in the WGS folder is treated as immutable input.
2. Original `*.vcf.gz`, `*.tbi`, `*.cram`, `*.crai` must never be modified.
3. Derived outputs are written only under `results/` and logs under `logs/`.
4. Pre/post SHA256 manifests are retained to verify input integrity.

## Repository Layout

- `config/`
- `config/paths.env.example`: template for local paths.
- `config/paths.env`: local/private paths (gitignored).
- `workflow/monogenic/`
- Core monogenic pipeline scripts.
- `workflow/utils/`
- Auxiliary utility scripts.
- `results/vep/`
- ClinVar-enriched VEP output(s).
- `results/monogenic/`
- Current triage outputs.
- `logs/`
- Checksum manifests and active run logs.

## Current Retained Outputs

- `results/vep/SGIHZBQI9ZH.vep.clinvar_enriched.vcf.gz`
- `results/vep/SGIHZBQI9ZH.vep.clinvar_enriched.vcf.gz_summary.html`
- `results/monogenic/high_confidence_clinvar_summary.md`
- `results/monogenic/high_confidence_clinvar_shortlist.tsv`
- `results/monogenic/high_confidence_clinvar_unique_variants.tsv`

First-pass non-enriched VEP output and redundant intermediate monogenic tables were removed during cleanup.

## Monogenic Workflow

1. Extract broad priority annotations from VEP output:

```bash
python3 workflow/monogenic/extract_priority_variants.py \
  --input-vcf results/vep/<sample>.vep.clinvar_enriched.vcf.gz \
  --output-tsv results/monogenic/priority_variants.tsv \
  --summary-txt results/monogenic/priority_variants.summary.txt
```

2. Build broad monogenic shortlists:

```bash
python3 workflow/monogenic/build_monogenic_shortlists.py \
  --input-tsv results/monogenic/priority_variants.tsv \
  --outdir results/monogenic
```

3. Build high-confidence ClinVar shortlist:

```bash
python3 workflow/monogenic/triage_clinvar_enriched.py \
  --input-vcf results/vep/<sample>.vep.clinvar_enriched.vcf.gz \
  --output-tsv results/monogenic/high_confidence_clinvar_shortlist.tsv \
  --summary-md results/monogenic/high_confidence_clinvar_summary.md \
  --max-af 0.01
```

4. Collapse to unique variants (canonical-preferred):

```bash
python3 workflow/monogenic/dedup_shortlist_by_variant.py \
  --input-tsv results/monogenic/high_confidence_clinvar_shortlist.tsv \
  --output-tsv results/monogenic/high_confidence_clinvar_unique_variants.tsv
```

## Environment

Recommended runtime:
- Windows 10/11 + WSL2 (Ubuntu)
- `python3`, `vep`, `bcftools`, `plink2`

## PRS Stage (Next)

This repo now includes a starter PRS workflow for `pgsc_calc`:

- `config/prs_samplesheet.example.csv`
- `config/prs_pgs_ids.example.txt`
- `workflow/prs/prepare_prs_samplesheet.sh`
- `workflow/prs/run_pgsc_calc.sh`

### 1) Prepare config files

```bash
cp config/prs_samplesheet.example.csv config/prs_samplesheet.csv
cp config/prs_pgs_ids.example.txt config/prs_pgs_ids.txt
```

Or auto-build the samplesheet from `config/paths.env`:

```bash
bash workflow/prs/prepare_prs_samplesheet.sh
```

Note:
- `prepare_prs_samplesheet.sh` now sanitizes `sampleset` to safe IDs (`A-Z`, `a-z`, `0-9`, `_`, `-`).

### 2) Run pgsc_calc

```bash
bash workflow/prs/run_pgsc_calc.sh
```

Run with a specific ID file (for batching):

```bash
bash workflow/prs/run_pgsc_calc.sh config/prs_batches/pgs_ids_batch_01.txt
```

Notes:
- Runtime profile is auto-selected by the script in this order:
  - `docker` -> `apptainer` -> `singularity` -> `conda`
- You can force a profile:
  - `PGSC_PROFILE=docker bash workflow/prs/run_pgsc_calc.sh`
  - `PGSC_PROFILE=conda bash workflow/prs/run_pgsc_calc.sh`
- You can adjust memory:
  - `PGSC_MAX_MEMORY=12.GB bash workflow/prs/run_pgsc_calc.sh`
- You can lower the match threshold (default `0.75`) when many scores fail matching:
  - `PGSC_MIN_OVERLAP=0.50 bash workflow/prs/run_pgsc_calc.sh config/prs_batches/pgs_ids_batch_01.txt`
- You can disable Nextflow cache reuse for a clean rerun:
  - `PGSC_RESUME=0 bash workflow/prs/run_pgsc_calc.sh config/prs_batches/pgs_ids_batch_01.txt`
- If using `config/prs_update_sex.tsv`, file must be strict `FID IID SEX` (no comments/blank lines).

Optional batching helper:

```bash
bash workflow/prs/split_pgs_ids.sh
```

Defaults:
- input: `config/prs_pgs_ids.txt`
- batch size: `10`
- output dir: `config/prs_batches/`

Outputs are written to:
- `results/prs/pgsc_calc/`

Note:
- `pgsc_calc` support and expectations for WGS/VCF inputs can evolve by version.
- Use `docs/PRS_EXECUTION_RUNBOOK.md` for failure recovery (`-resume`, lock cleanup flow, chrX handling context).

## Disclaimer

This repository supports research/hobby analysis only and is not a medical diagnostic system.
