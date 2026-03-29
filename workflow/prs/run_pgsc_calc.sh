#!/usr/bin/env bash
set -euo pipefail

# Run pgsc_calc from explicit inputs.
#
# Required:
# - samplesheet CSV
# - PGS IDs file
# - nextflow in PATH
# - runtime for selected profile:
#   - docker (preferred) OR
#   - apptainer/singularity OR
#   - conda (fallback)
#
# Usage:
#   bash workflow/prs/run_pgsc_calc.sh --samplesheet sample.csv --pgs-ids-file ids.txt [options]

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

SAMPLESHEET="${PGSC_SAMPLESHEET:-}"
PGS_IDS_FILE="${PGSC_PGS_IDS_FILE:-}"
SEX_UPDATE_FILE="${PGSC_SEX_UPDATE_FILE:-}"
OUTDIR="${PGSC_OUTDIR:-output/results/prs/pgsc_calc}"
WORKDIR="${PGSC_WORKDIR:-output/work/pgsc_calc}"
PROFILE="${PGSC_PROFILE:-auto}"
MAX_MEMORY="${PGSC_MAX_MEMORY:-14.GB}"
MIN_OVERLAP="${PGSC_MIN_OVERLAP:-0.75}"
RESUME_FLAG="${PGSC_RESUME:-1}"
GENOTYPES_CACHE="${PGSC_GENOTYPES_CACHE:-}"

usage() {
  cat <<EOF
Usage: bash workflow/prs/run_pgsc_calc.sh --samplesheet sample.csv --pgs-ids-file ids.txt [options]

Options:
  --samplesheet PATH
  --pgs-ids-file PATH
  --sex-update-file PATH
  --outdir PATH
  --workdir PATH
  --profile docker|apptainer|singularity|conda|auto
  --max-memory 14.GB
  --min-overlap 0.75
  --resume 0|1
  --genotypes-cache PATH
EOF
}

POSITIONAL_PGS_IDS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet)
      SAMPLESHEET="${2:-}"
      shift 2
      ;;
    --pgs-ids-file)
      PGS_IDS_FILE="${2:-}"
      shift 2
      ;;
    --sex-update-file)
      SEX_UPDATE_FILE="${2:-}"
      shift 2
      ;;
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --workdir)
      WORKDIR="${2:-}"
      shift 2
      ;;
    --profile)
      PROFILE="${2:-}"
      shift 2
      ;;
    --max-memory)
      MAX_MEMORY="${2:-}"
      shift 2
      ;;
    --min-overlap)
      MIN_OVERLAP="${2:-}"
      shift 2
      ;;
    --resume)
      RESUME_FLAG="${2:-}"
      shift 2
      ;;
    --genotypes-cache)
      GENOTYPES_CACHE="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "$POSITIONAL_PGS_IDS" ]]; then
        POSITIONAL_PGS_IDS="$1"
        shift
      else
        echo "Unknown argument: $1"
        usage
        exit 1
      fi
      ;;
  esac
done

if [[ -n "$POSITIONAL_PGS_IDS" && -z "$PGS_IDS_FILE" ]]; then
  PGS_IDS_FILE="$POSITIONAL_PGS_IDS"
fi

if [[ -z "$SAMPLESHEET" ]]; then
  echo "--samplesheet is required"
  usage
  exit 1
fi

if [[ -z "$PGS_IDS_FILE" ]]; then
  echo "--pgs-ids-file is required"
  usage
  exit 1
fi

if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "Missing $SAMPLESHEET"
  exit 1
fi

if [[ ! -f "$PGS_IDS_FILE" ]]; then
  echo "Missing $PGS_IDS_FILE"
  exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
  echo "nextflow not found in PATH"
  exit 1
fi

resolve_profile() {
  local requested="$1"
  case "$requested" in
    docker|apptainer|singularity|conda)
      echo "$requested"
      return 0
      ;;
    auto)
      if command -v docker >/dev/null 2>&1; then
        echo "docker"
        return 0
      fi
      if command -v apptainer >/dev/null 2>&1; then
        echo "apptainer"
        return 0
      fi
      if command -v singularity >/dev/null 2>&1; then
        echo "singularity"
        return 0
      fi
      if command -v conda >/dev/null 2>&1; then
        echo "conda"
        return 0
      fi
      echo ""
      return 1
      ;;
    *)
      echo "Unsupported PGSC_PROFILE: $requested"
      echo "Use one of: auto, docker, apptainer, singularity, conda"
      return 1
      ;;
  esac
}

PROFILE="$(resolve_profile "$PROFILE")" || {
  echo "No supported runtime found."
  echo "Install one runtime and retry:"
  echo "  - docker (recommended)"
  echo "  - apptainer/singularity"
  echo "  - conda (fallback)"
  exit 1
}

case "$PROFILE" in
  docker)
    if ! command -v docker >/dev/null 2>&1; then
      echo "docker profile selected but docker not found in PATH"
      exit 1
    fi
    ;;
  apptainer)
    if ! command -v apptainer >/dev/null 2>&1; then
      echo "apptainer profile selected but apptainer not found in PATH"
      exit 1
    fi
    ;;
  singularity)
    if ! command -v singularity >/dev/null 2>&1; then
      echo "singularity profile selected but singularity not found in PATH"
      exit 1
    fi
    ;;
  conda)
    if ! command -v conda >/dev/null 2>&1; then
      echo "conda profile selected but conda not found in PATH"
      exit 1
    fi
    ;;
esac

echo "Using pgsc_calc runtime profile: $PROFILE"
echo "Using min overlap threshold: $MIN_OVERLAP"
if [[ "$RESUME_FLAG" == "1" ]]; then
  echo "Using Nextflow resume: enabled"
else
  echo "Using Nextflow resume: disabled"
fi
if [[ "$PROFILE" == "conda" ]]; then
  echo "Note: pgsc_calc v2.2.0 currently ships an unsatisfiable conda pin"
  echo "(pgscatalog.core=1.0.0 with pgscatalog.calc=0.3.0)."
  echo "If conda solve fails, patch ~/.nextflow/assets/pgscatalog/pgsc_calc/environments/pgscatalog_utils/environment.yml"
  echo "and change pgscatalog.calc=0.3.0 -> pgscatalog.calc=0.3.1, then rerun."
fi

mkdir -p "$OUTDIR"

# Collapse PGS IDs into a comma-separated list for the pipeline argument.
# Trim whitespace, ignore comments/blank lines, strip trailing inline comments,
# and strip CR characters to tolerate Windows-edited files.
PGS_IDS="$(
  awk '
    {
      sub(/\r$/, "")
      sub(/[[:space:]]+#.*$/, "", $0)
    }
    /^[[:space:]]*#/ {
      next
    }
    /^[[:space:]]*$/ {
      next
    }
    {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      print
    }
  ' "$PGS_IDS_FILE" | paste -sd, -
)"

if [[ -z "$PGS_IDS" ]]; then
  echo "No PGS IDs found in $PGS_IDS_FILE"
  exit 1
fi

echo "Running pgsc_calc with PGS IDs: $PGS_IDS"

EXTRA_NEXTFLOW_ARGS=()
mkdir -p "$WORKDIR"
OVERRIDE_CONFIG="$WORKDIR/prs_plink2_override.config"
PLINK2_EXTRA_ARGS="--new-id-max-allele-len 100 missing --split-par b38"
if [[ -f "$SEX_UPDATE_FILE" ]]; then
  # Require a clean "FID IID SEX" file to avoid hidden runtime rewrites.
  if [[ ! -s "$SEX_UPDATE_FILE" ]]; then
    echo "Sex update file is empty: $SEX_UPDATE_FILE"
    exit 1
  fi

  # Validate "FID IID SEX" with SEX in {0,1,2}; no comments/blank lines allowed.
  if ! awk 'NF != 3 || ($3 != 0 && $3 != 1 && $3 != 2) { exit 1 }' "$SEX_UPDATE_FILE"; then
    echo "Invalid format in $SEX_UPDATE_FILE"
    echo "Expected: FID IID SEX   (SEX must be 0, 1, or 2)"
    echo "No comments or blank lines are allowed in this file."
    exit 1
  fi

  sex_update_path="$(realpath "$SEX_UPDATE_FILE")"
  escaped_sex_update_path="$(printf "%s" "$sex_update_path" | sed "s/'/'\\\\''/g")"
  PLINK2_EXTRA_ARGS="$PLINK2_EXTRA_ARGS --update-sex '$escaped_sex_update_path'"
  echo "Using sex update file: $SEX_UPDATE_FILE"
fi
cat > "$OVERRIDE_CONFIG" <<EOF
process {
  withName: PLINK2_VCF {
    ext.args = "${PLINK2_EXTRA_ARGS}"
  }
}
EOF
EXTRA_NEXTFLOW_ARGS+=(-c "$OVERRIDE_CONFIG")
if [[ -n "$GENOTYPES_CACHE" ]]; then
  EXTRA_NEXTFLOW_ARGS+=(--genotypes_cache "$GENOTYPES_CACHE")
fi

nextflow run pgscatalog/pgsc_calc \
  -profile "$PROFILE" \
  --input "$SAMPLESHEET" \
  --pgs_id "$PGS_IDS" \
  --target_build GRCh38 \
  --min_overlap "$MIN_OVERLAP" \
  --max_memory "$MAX_MEMORY" \
  --outdir "$OUTDIR" \
  -work-dir "$WORKDIR" \
  "${EXTRA_NEXTFLOW_ARGS[@]}" \
  $([[ "$RESUME_FLAG" == "1" ]] && echo "-resume")

echo "Done. Results in: $OUTDIR"
