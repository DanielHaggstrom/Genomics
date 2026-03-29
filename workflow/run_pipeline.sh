#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PIPELINE_ENV="config/pipeline.env"
SECRETS_ENV="config/secrets.env"
PATHS_ENV="config/paths.env"
STAGE="all"

usage() {
  cat <<EOF
Usage: bash workflow/run_pipeline.sh [--stage all|monogenic|prs]

Environment files:
  - config/pipeline.env (required; copy from config/pipeline.env.example)
  - config/secrets.env (optional; copy from config/secrets.env.example)
  - .env (optional; copy from .env.example)
  - config/paths.env (optional; used by PRS samplesheet auto-build and monogenic VEP auto-annotation)
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

if [[ -f "$PATHS_ENV" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$PATHS_ENV"
  set +a
fi

SAMPLE_ID="${SAMPLE_ID:-}"
ANNOTATED_VCF_GZ="${ANNOTATED_VCF_GZ:-}"
RUN_MONOGENIC="${RUN_MONOGENIC:-1}"
RUN_PRS="${RUN_PRS:-1}"
MONOGENIC_MODE="${MONOGENIC_MODE:-moderate}"
MONOGENIC_MAX_AF="${MONOGENIC_MAX_AF:-0.05}"
MONOGENIC_AUTORUN_VEP="${MONOGENIC_AUTORUN_VEP:-0}"
MONOGENIC_SOURCE_VCF_GZ="${MONOGENIC_SOURCE_VCF_GZ:-${VCF_GZ:-}}"
MONOGENIC_VEP_DATA_DIR="${MONOGENIC_VEP_DATA_DIR:-}"
MONOGENIC_VEP_FORKS="${MONOGENIC_VEP_FORKS:-8}"
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

  local annotated_vcf="$ANNOTATED_VCF_GZ"
  if [[ "$MONOGENIC_AUTORUN_VEP" == "1" ]]; then
    if [[ -z "$MONOGENIC_SOURCE_VCF_GZ" ]]; then
      echo "MONOGENIC_SOURCE_VCF_GZ is required when MONOGENIC_AUTORUN_VEP=1"
      echo "Set it directly or provide VCF_GZ via $PATHS_ENV."
      exit 1
    fi
    if [[ ! -f "$MONOGENIC_SOURCE_VCF_GZ" ]]; then
      echo "Monogenic source VCF not found: $MONOGENIC_SOURCE_VCF_GZ"
      exit 1
    fi
    if [[ -z "$annotated_vcf" ]]; then
      local source_dir
      source_dir="$(dirname "$MONOGENIC_SOURCE_VCF_GZ")"
      local source_base
      source_base="$(basename "$MONOGENIC_SOURCE_VCF_GZ" .vcf.gz)"
      annotated_vcf="${source_dir}/vep_data/output/${source_base}.vep.clinvar_enriched.vcf.gz"
    fi

    local vep_args=(
      --input-vcf "$MONOGENIC_SOURCE_VCF_GZ"
      --output-vcf "$annotated_vcf"
      --forks "$MONOGENIC_VEP_FORKS"
    )
    if [[ -n "$MONOGENIC_VEP_DATA_DIR" ]]; then
      vep_args+=(--vep-data-dir "$MONOGENIC_VEP_DATA_DIR")
    fi

    bash workflow/monogenic/annotate_with_vep.sh "${vep_args[@]}"
  fi

  if [[ -z "$annotated_vcf" ]]; then
    echo "ANNOTATED_VCF_GZ is required for monogenic stage"
    echo "Or set MONOGENIC_AUTORUN_VEP=1 with MONOGENIC_SOURCE_VCF_GZ."
    exit 1
  fi
  if [[ ! -f "$annotated_vcf" ]]; then
    echo "Monogenic input not found: $annotated_vcf"
    exit 1
  fi

  local outdir="results/monogenic/${SAMPLE_ID}"
  mkdir -p "$outdir"
  rm -f \
    "$outdir/shortlist_summary.md" \
    "$outdir/high_confidence_clinvar_summary.md" \
    "$outdir/high_confidence_clinvar_shortlist.tsv" \
    "$outdir/high_confidence_clinvar_unique_variants.tsv"

  python3 workflow/monogenic/extract_priority_variants.py \
    --input-vcf "$annotated_vcf" \
    --output-tsv "$outdir/priority_variants.tsv" \
    --summary-txt "$outdir/priority_variants.summary.txt"

  python3 workflow/monogenic/build_monogenic_shortlists.py \
    --input-tsv "$outdir/priority_variants.tsv" \
    --outdir "$outdir"

  python3 workflow/monogenic/triage_clinvar_enriched.py \
    --input-vcf "$annotated_vcf" \
    --output-tsv "$outdir/rare_clinvar_review_queue.tsv" \
    --mode "$MONOGENIC_MODE" \
    --max-af "$MONOGENIC_MAX_AF"

  python3 workflow/monogenic/dedup_shortlist_by_variant.py \
    --input-tsv "$outdir/rare_clinvar_review_queue.tsv" \
    --output-tsv "$outdir/rare_clinvar_review_unique_variants.tsv"

  python3 workflow/monogenic/render_monogenic_report.py \
    --sample-id "$SAMPLE_ID" \
    --output-md "$outdir/monogenic_report.md" \
    --priority-summary "$outdir/priority_variants.summary.txt" \
    --clinvar-unique-tsv "$outdir/clinvar_plp_unique_variants.tsv" \
    --lof-unique-tsv "$outdir/lof_unique_variants.tsv" \
    --review-unique-tsv "$outdir/rare_clinvar_review_unique_variants.tsv" \
    --mode "$MONOGENIC_MODE" \
    --max-af "$MONOGENIC_MAX_AF"

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
