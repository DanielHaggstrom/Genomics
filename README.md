# Genomics WGS Risk Analysis Pipeline

Local, privacy-preserving workflow for personal whole-genome sequencing (WGS) analysis.

Primary objective:
- monogenic variant triage with ClinVar evidence and transcript-aware filtering.

Planned next objective:
- polygenic risk score (PRS) computation using PGS Catalog tooling.

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

- `results/vep/DEMO_SAMPLE_001.vep.clinvar_enriched.vcf.gz`
- `results/vep/DEMO_SAMPLE_001.vep.clinvar_enriched.vcf.gz_summary.html`
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

## Next Planned Step

Initialize PRS stage:
- target genotype conversion with `plink2`
- `pgsc_calc` setup
- first trait batch: CAD, T2D, BMI, Alzheimer’s, lipids

## Disclaimer

This repository supports research/hobby analysis only and is not a medical diagnostic system.
