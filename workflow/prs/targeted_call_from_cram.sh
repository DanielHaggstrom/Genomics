#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: bash workflow/prs/targeted_call_from_cram.sh --cram sample.cram --reference-fasta ref.fa --targets-vcf union.vcf.gz --outdir outdir --sample-id SAMPLE
EOF
}

CRAM=""
REFERENCE_FASTA=""
TARGETS_VCF=""
OUTDIR=""
SAMPLE_ID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cram) CRAM="${2:-}"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="${2:-}"; shift 2 ;;
    --targets-vcf) TARGETS_VCF="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --sample-id) SAMPLE_ID="${2:-}"; shift 2 ;;
    *) usage; exit 1 ;;
  esac
done

if [[ -z "$CRAM" || -z "$REFERENCE_FASTA" || -z "$TARGETS_VCF" || -z "$OUTDIR" || -z "$SAMPLE_ID" ]]; then
  usage
  exit 1
fi

mkdir -p "$OUTDIR"
RAW_VCF="$OUTDIR/${SAMPLE_ID}.prs_targeted.raw.vcf.gz"
NORMALIZED_VCF="$OUTDIR/${SAMPLE_ID}.prs_targeted.vcf.gz"

if [[ ! -f "$RAW_VCF" ]]; then
  bcftools mpileup \
    -Ou \
    -f "$REFERENCE_FASTA" \
    -T "$TARGETS_VCF" \
    -a FORMAT/DP,FORMAT/AD,FORMAT/GQ \
    "$CRAM" \
  | bcftools call \
      -Aim \
      --insert-missed \
      -C alleles \
      -T "$TARGETS_VCF" \
      -Oz \
      -o "$RAW_VCF"
fi

bcftools norm \
  -f "$REFERENCE_FASTA" \
  -m -both \
  -Oz \
  -o "$NORMALIZED_VCF" \
  "$RAW_VCF"

tabix -f -p vcf "$RAW_VCF"
tabix -f -p vcf "$NORMALIZED_VCF"

echo "$NORMALIZED_VCF"
