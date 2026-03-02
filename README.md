# Genomics WGS Risk Analysis Pipeline

Structured, reproducible local pipeline for personal WGS analysis with two stages:
- Monogenic triage from ClinVar-enriched VEP output.
- Polygenic risk scoring (PRS) via `pgscatalog/pgsc_calc`.

This repository is publish-safe by default:
- Real secrets and local machine paths are ignored by git.
- Example config files with fake values are tracked.

## Project Layout

- `workflow/run_pipeline.sh`: single entrypoint orchestrating stages.
- `workflow/monogenic/`: monogenic analysis scripts.
- `workflow/prs/`: PRS prep and execution scripts.
- `workflow/utils/`: utility helpers.
- `config/*.example*`: tracked templates.
- `docs/`: architecture and runbooks.
- `results/`, `logs/`, `work/`, `.nextflow/`: generated artifacts (ignored).

## Quick Start

1. Bootstrap local config files:

```bash
make bootstrap
```

2. Edit local files (ignored by git):
- `config/pipeline.env`
- `config/secrets.env`
- `.env`
- `config/prs_samplesheet.csv` (or set auto-build)
- `config/prs_pgs_ids.txt`

3. Run validation checks:

```bash
make validate
make prepublish
```

4. Run pipeline:

```bash
# full pipeline
make pipeline

# monogenic only
make monogenic

# PRS only
make prs
```

## Pipeline Configuration

Primary runtime config: `config/pipeline.env`

Key variables:
- `SAMPLE_ID`: logical sample identifier for outputs.
- `ANNOTATED_VCF_GZ`: path to ClinVar-enriched VEP VCF (`.vcf.gz`).
- `RUN_MONOGENIC` / `RUN_PRS`: enable stages in `--stage all` mode.
- `MONOGENIC_MODE`, `MONOGENIC_MAX_AF`: triage strictness.
- `PRS_AUTOBUILD_SAMPLESHEET`: set `1` to build from `config/paths.env`.
- `PRS_SAMPLESHEET`, `PRS_PGS_IDS_FILE`: PRS inputs.
- `PGSC_PROFILE`, `PGSC_MAX_MEMORY`, `PGSC_MIN_OVERLAP`, `PGSC_RESUME`: pgsc_calc runtime controls.

Optional secret files:
- `.env` (from `.env.example`)
- `config/secrets.env` (from `config/secrets.env.example`)

Both are auto-loaded if present.

## Security and Publishing

Ignored by default:
- Real env/secrets files (`.env`, `config/secrets.env`, `config/pipeline.env`, etc.)
- Local path configs (`config/paths.env`, `config/prs_*.csv|txt|tsv` local variants)
- Runtime/cache/output artifacts

Before publishing:
1. Confirm `git status` does not include private/local files.
2. Keep only `*.example*` templates for configuration.
3. Do not commit generated `results/` or `logs/`.

## CI

GitHub Actions workflow at `.github/workflows/ci.yml` runs:
- Python syntax checks
- Shell lint and syntax checks
- Pipeline help/entrypoint check

## Disclaimer

This project is for research/hobby analysis and is not a medical diagnostic system.

## License

MIT. See [LICENSE](LICENSE).
