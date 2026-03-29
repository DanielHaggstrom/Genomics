#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: bash workflow/core/call_variants_from_cram.sh --cram sample.cram --reference-fasta ref.fa --outdir outdir --sample-id SAMPLE
EOF
}

CRAM=""
REFERENCE_FASTA=""
OUTDIR=""
SAMPLE_ID=""
GATK_IMAGE="${GATK_IMAGE:-broadinstitute/gatk:4.6.1.0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cram) CRAM="${2:-}"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --sample-id) SAMPLE_ID="${2:-}"; shift 2 ;;
    *) usage; exit 1 ;;
  esac
done

if [[ -z "$CRAM" || -z "$REFERENCE_FASTA" || -z "$OUTDIR" || -z "$SAMPLE_ID" ]]; then
  usage
  exit 1
fi

mkdir -p "$OUTDIR"

CRAM_DIR="$(cd "$(dirname "$CRAM")" && pwd)"
CRAM_BASENAME="$(basename "$CRAM")"
REF_DIR="$(cd "$(dirname "$REFERENCE_FASTA")" && pwd)"
REF_BASENAME="$(basename "$REFERENCE_FASTA")"
OUTDIR_ABS="$(cd "$OUTDIR" && pwd)"

GVCF="$OUTDIR_ABS/${SAMPLE_ID}.g.vcf.gz"
RAW_VCF="$OUTDIR_ABS/${SAMPLE_ID}.raw.vcf.gz"
FILTERED_VCF="$OUTDIR_ABS/${SAMPLE_ID}.hard-filtered.vcf.gz"

run_gatk() {
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$CRAM_DIR:/input" \
    -v "$REF_DIR:/ref" \
    -v "$OUTDIR_ABS:/out" \
    "$GATK_IMAGE" "$@"
}

if [[ ! -f "$GVCF" ]]; then
  run_gatk gatk HaplotypeCaller \
    -R "/ref/$REF_BASENAME" \
    -I "/input/$CRAM_BASENAME" \
    -O "/out/$(basename "$GVCF")" \
    -ERC GVCF
fi

if [[ ! -f "$RAW_VCF" ]]; then
  run_gatk gatk GenotypeGVCFs \
    -R "/ref/$REF_BASENAME" \
    -V "/out/$(basename "$GVCF")" \
    -O "/out/$(basename "$RAW_VCF")"
fi

TMP_SNP="$OUTDIR_ABS/${SAMPLE_ID}.snps.filtered.vcf.gz"
TMP_INDEL="$OUTDIR_ABS/${SAMPLE_ID}.indels.filtered.vcf.gz"

if [[ ! -f "$TMP_SNP" ]]; then
  run_gatk gatk SelectVariants \
    -R "/ref/$REF_BASENAME" \
    -V "/out/$(basename "$RAW_VCF")" \
    --select-type-to-include SNP \
    -O "/out/${SAMPLE_ID}.snps.vcf.gz"

  run_gatk gatk VariantFiltration \
    -R "/ref/$REF_BASENAME" \
    -V "/out/${SAMPLE_ID}.snps.vcf.gz" \
    -O "/out/${SAMPLE_ID}.snps.filtered.vcf.gz" \
    --filter-name "SNP_HARD_FILTER" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0"
fi

if [[ ! -f "$TMP_INDEL" ]]; then
  run_gatk gatk SelectVariants \
    -R "/ref/$REF_BASENAME" \
    -V "/out/$(basename "$RAW_VCF")" \
    --select-type-to-include INDEL \
    -O "/out/${SAMPLE_ID}.indels.vcf.gz"

  run_gatk gatk VariantFiltration \
    -R "/ref/$REF_BASENAME" \
    -V "/out/${SAMPLE_ID}.indels.vcf.gz" \
    -O "/out/${SAMPLE_ID}.indels.filtered.vcf.gz" \
    --filter-name "INDEL_HARD_FILTER" \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"
fi

bcftools concat -a -Oz -o "$FILTERED_VCF" "$TMP_SNP" "$TMP_INDEL"
tabix -f -p vcf "$FILTERED_VCF"

echo "$FILTERED_VCF"
