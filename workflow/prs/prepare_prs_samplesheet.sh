#!/usr/bin/env bash
set -euo pipefail

# Build a pgsc_calc samplesheet from config/paths.env.
# Output: config/prs_samplesheet.csv

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

if [[ -f "config/pipeline.env" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "config/pipeline.env"
  set +a
fi

if [[ ! -f "config/paths.env" ]]; then
  echo "Missing config/paths.env"
  exit 1
fi

# shellcheck source=/dev/null
source "config/paths.env"

if [[ -z "${VCF_GZ:-}" ]]; then
  echo "VCF_GZ is not set in config/paths.env"
  exit 1
fi

# pgsc_calc expects path_prefix without the .vcf.gz extension.
path_prefix="${VCF_GZ%.vcf.gz}"

if [[ "$path_prefix" == "$VCF_GZ" ]]; then
  echo "VCF_GZ should end with .vcf.gz (got: $VCF_GZ)"
  exit 1
fi

sample_name="${SAMPLE_ID:-}"
if [[ -z "$sample_name" ]]; then
  sample_name="$(basename "$path_prefix")"
fi

# pgsc_calc reserves the sampleset name "reference" and rejects underscores.
# Prefer SAMPLE_ID from config/pipeline.env when available and normalize any
# fallback name into a compact alphanumeric/hyphen form.
sample_name="$(echo "$sample_name" | sed -E 's/[^A-Za-z0-9-]+/-/g; s/^-+//; s/-+$//; s/-+/-/g')"

if [[ -z "$sample_name" ]]; then
  echo "Could not derive a valid sampleset name"
  exit 1
fi
if [[ "${sample_name,,}" == "reference" ]]; then
  echo "sampleset name 'reference' is reserved by pgsc_calc"
  exit 1
fi

cat > config/prs_samplesheet.csv <<EOF
sampleset,path_prefix,chrom,format,vcf_genotype_field
${sample_name},${path_prefix},,vcf,GT
EOF

echo "Wrote config/prs_samplesheet.csv"
