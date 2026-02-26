#!/usr/bin/env bash
set -euo pipefail

# Build a pgsc_calc samplesheet from config/paths.env.
# Output: config/prs_samplesheet.csv

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

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

sample_name="$(basename "$path_prefix")"
# Keep sample IDs simple/safe for downstream tools.
# Replace characters outside [A-Za-z0-9_-] with underscore.
sample_name="$(echo "$sample_name" | sed 's/[^A-Za-z0-9_-]/_/g')"

cat > config/prs_samplesheet.csv <<EOF
sampleset,path_prefix,chrom,format,vcf_genotype_field
${sample_name},${path_prefix},,vcf,GT
EOF

echo "Wrote config/prs_samplesheet.csv"
