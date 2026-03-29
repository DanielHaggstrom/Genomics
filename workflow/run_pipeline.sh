#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PIPELINE_ENV="config/pipeline.env"
SECRETS_ENV="config/secrets.env"
STAGE="all"

usage() {
  cat <<EOF
Usage: bash workflow/run_pipeline.sh [--stage all|monogenic|prs]

Environment files:
  - config/pipeline.env (required; copy from config/pipeline.env.example)
  - config/secrets.env (optional; copy from config/secrets.env.example)
  - .env (optional; copy from .env.example)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage)
      if [[ $# -lt 2 || "$2" == -* ]]; then
        echo "--stage requires one of: all, monogenic, prs"
        usage
        exit 1
      fi
      STAGE="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ ! -f "$PIPELINE_ENV" ]]; then
  echo "Missing $PIPELINE_ENV"
  echo "Create it from config/pipeline.env.example"
  exit 1
fi

if [[ -f ".env" ]]; then
  set -a
  # shellcheck source=/dev/null
  source ".env"
  set +a
fi

if [[ -f "$SECRETS_ENV" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$SECRETS_ENV"
  set +a
fi

set -a
# shellcheck source=/dev/null
source "$PIPELINE_ENV"
set +a

SAMPLE_ID="${SAMPLE_ID:-}"
ANNOTATED_VCF_GZ="${ANNOTATED_VCF_GZ:-}"
RUN_MONOGENIC="${RUN_MONOGENIC:-1}"
RUN_PRS="${RUN_PRS:-1}"
MONOGENIC_MODE="${MONOGENIC_MODE:-moderate}"
MONOGENIC_MAX_AF="${MONOGENIC_MAX_AF:-0.05}"
PRS_AUTOBUILD_SAMPLESHEET="${PRS_AUTOBUILD_SAMPLESHEET:-0}"
PRS_SAMPLESHEET="${PRS_SAMPLESHEET:-config/prs_samplesheet.csv}"
PRS_PGS_IDS_FILE="${PRS_PGS_IDS_FILE:-config/prs_pgs_ids.txt}"
PGSC_PROFILE="${PGSC_PROFILE:-auto}"
PGSC_MAX_MEMORY="${PGSC_MAX_MEMORY:-14.GB}"
PGSC_MIN_OVERLAP="${PGSC_MIN_OVERLAP:-0.75}"
PGSC_RESUME="${PGSC_RESUME:-1}"

run_monogenic() {
  if [[ -z "$SAMPLE_ID" ]]; then
    echo "SAMPLE_ID is required in $PIPELINE_ENV for monogenic outputs"
    exit 1
  fi
  if [[ -z "$ANNOTATED_VCF_GZ" ]]; then
    echo "ANNOTATED_VCF_GZ is required for monogenic stage"
    exit 1
  fi
  if [[ ! -f "$ANNOTATED_VCF_GZ" ]]; then
    echo "Monogenic input not found: $ANNOTATED_VCF_GZ"
    exit 1
  fi

  local outdir="results/monogenic/${SAMPLE_ID}"
  mkdir -p "$outdir"

  python3 workflow/monogenic/extract_priority_variants.py \
    --input-vcf "$ANNOTATED_VCF_GZ" \
    --output-tsv "$outdir/priority_variants.tsv" \
    --summary-txt "$outdir/priority_variants.summary.txt"

  python3 workflow/monogenic/build_monogenic_shortlists.py \
    --input-tsv "$outdir/priority_variants.tsv" \
    --outdir "$outdir"

  python3 workflow/monogenic/triage_clinvar_enriched.py \
    --input-vcf "$ANNOTATED_VCF_GZ" \
    --output-tsv "$outdir/high_confidence_clinvar_shortlist.tsv" \
    --summary-md "$outdir/high_confidence_clinvar_summary.md" \
    --mode "$MONOGENIC_MODE" \
    --max-af "$MONOGENIC_MAX_AF"

  python3 workflow/monogenic/dedup_shortlist_by_variant.py \
    --input-tsv "$outdir/high_confidence_clinvar_shortlist.tsv" \
    --output-tsv "$outdir/high_confidence_clinvar_unique_variants.tsv"

  echo "Monogenic outputs: $outdir"
}

run_prs() {
  if [[ "$PRS_AUTOBUILD_SAMPLESHEET" == "1" ]]; then
    bash workflow/prs/prepare_prs_samplesheet.sh
  elif [[ ! -f "$PRS_SAMPLESHEET" ]]; then
    echo "Missing PRS samplesheet: $PRS_SAMPLESHEET"
    echo "Set PRS_AUTOBUILD_SAMPLESHEET=1 or create PRS_SAMPLESHEET."
    exit 1
  fi

  PGSC_SAMPLESHEET="$PRS_SAMPLESHEET" \
  PGSC_PROFILE="$PGSC_PROFILE" \
  PGSC_MAX_MEMORY="$PGSC_MAX_MEMORY" \
  PGSC_MIN_OVERLAP="$PGSC_MIN_OVERLAP" \
  PGSC_RESUME="$PGSC_RESUME" \
  PGSC_PGS_IDS_FILE="$PRS_PGS_IDS_FILE" \
  bash workflow/prs/run_pgsc_calc.sh
}

case "$STAGE" in
  all)
    [[ "$RUN_MONOGENIC" == "1" ]] && run_monogenic
    [[ "$RUN_PRS" == "1" ]] && run_prs
    ;;
  monogenic)
    run_monogenic
    ;;
  prs)
    run_prs
    ;;
  *)
    echo "Invalid --stage value: $STAGE"
    usage
    exit 1
    ;;
esac

echo "Pipeline run finished."
