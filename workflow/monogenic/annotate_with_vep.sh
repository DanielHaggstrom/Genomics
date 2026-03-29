#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

INPUT_VCF_GZ=""
OUTPUT_VCF_GZ=""
VEP_DATA_DIR=""
VEP_FORKS="${MONOGENIC_VEP_FORKS:-8}"
DOCKER_BIN="${MONOGENIC_DOCKER_BIN:-}"
VEP_IMAGE="${MONOGENIC_VEP_IMAGE:-ensemblorg/ensembl-vep:release_115.0}"
CLINVAR_URL="${MONOGENIC_CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz}"
CLINVAR_TBI_URL="${MONOGENIC_CLINVAR_TBI_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi}"
KEEP_CACHE_TARBALL="${MONOGENIC_KEEP_VEP_TARBALL:-0}"
VEP_CACHE_VERSION_OVERRIDE="${MONOGENIC_VEP_CACHE_VERSION:-}"

usage() {
  cat <<EOF
Usage: bash workflow/monogenic/annotate_with_vep.sh --input-vcf path/to/sample.vcf.gz --output-vcf path/to/output.vcf.gz [--vep-data-dir path/to/vep_data] [--forks 8]

Requirements:
  - Input VCF must be bgzip-compressed (.vcf.gz)
  - Docker Desktop / docker engine reachable via docker.exe or docker
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-vcf)
      INPUT_VCF_GZ="${2:-}"
      shift 2
      ;;
    --output-vcf)
      OUTPUT_VCF_GZ="${2:-}"
      shift 2
      ;;
    --vep-data-dir)
      VEP_DATA_DIR="${2:-}"
      shift 2
      ;;
    --forks)
      VEP_FORKS="${2:-}"
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

if [[ -z "$INPUT_VCF_GZ" || -z "$OUTPUT_VCF_GZ" ]]; then
  usage
  exit 1
fi

if [[ ! -f "$INPUT_VCF_GZ" ]]; then
  echo "Input VCF not found: $INPUT_VCF_GZ"
  exit 1
fi

if [[ "$INPUT_VCF_GZ" != *.vcf.gz ]]; then
  echo "Input VCF must end with .vcf.gz: $INPUT_VCF_GZ"
  exit 1
fi

if [[ -z "$VEP_DATA_DIR" ]]; then
  VEP_DATA_DIR="$(dirname "$INPUT_VCF_GZ")/vep_data"
fi

SOURCE_DIR="$(cd "$(dirname "$INPUT_VCF_GZ")" && pwd)"
mkdir -p "$(dirname "$OUTPUT_VCF_GZ")"
OUTPUT_DIR="$(cd "$(dirname "$OUTPUT_VCF_GZ")" && pwd)"
VEP_DATA_DIR="$(mkdir -p "$VEP_DATA_DIR" && cd "$VEP_DATA_DIR" && pwd)"

case "$OUTPUT_DIR" in
  "$SOURCE_DIR" | "$SOURCE_DIR"/*)
    ;;
  *)
    echo "Output VCF must live under the same directory as the source VCF."
    echo "Source directory: $SOURCE_DIR"
    echo "Output directory: $OUTPUT_DIR"
    exit 1
    ;;
esac

case "$VEP_DATA_DIR" in
  "$SOURCE_DIR" | "$SOURCE_DIR"/*)
    ;;
  *)
    echo "VEP data directory must live under the same directory as the source VCF."
    echo "Source directory: $SOURCE_DIR"
    echo "VEP data directory: $VEP_DATA_DIR"
    exit 1
    ;;
esac

if [[ -z "$DOCKER_BIN" ]]; then
  if command -v docker.exe >/dev/null 2>&1; then
    DOCKER_BIN="docker.exe"
  elif command -v docker >/dev/null 2>&1; then
    DOCKER_BIN="docker"
  else
    echo "Neither docker.exe nor docker is available in PATH."
    exit 1
  fi
fi

ensure_docker_ready() {
  if "$DOCKER_BIN" version --format '{{.Server.Version}}' >/dev/null 2>&1; then
    return 0
  fi

  if [[ "$DOCKER_BIN" == "docker.exe" ]] && command -v powershell.exe >/dev/null 2>&1; then
    powershell.exe -NoProfile -Command "& { \$p='C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe'; if (Test-Path \$p) { Start-Process -FilePath \$p } }" >/dev/null 2>&1 || true
    for _ in {1..24}; do
      sleep 5
      if "$DOCKER_BIN" version --format '{{.Server.Version}}' >/dev/null 2>&1; then
        return 0
      fi
    done
  fi

  echo "Docker is not ready. Start Docker Desktop and retry."
  exit 1
}

detect_cache_version() {
  if [[ -n "$VEP_CACHE_VERSION_OVERRIDE" ]]; then
    printf '%s\n' "$VEP_CACHE_VERSION_OVERRIDE"
    return 0
  fi

  local version_output
  version_output="$("$DOCKER_BIN" run --rm "$VEP_IMAGE" vep --help 2>&1 || true)"
  local cache_version
  cache_version="$(printf '%s\n' "$version_output" | sed -n 's/.*ensembl-vep[[:space:]]*:[[:space:]]*\([0-9][0-9]*\)\..*/\1/p' | head -n 1)"
  if [[ -z "$cache_version" ]]; then
    cache_version="115"
  fi
  printf '%s\n' "$cache_version"
}

download_if_missing() {
  local url="$1"
  local target="$2"
  if [[ -f "$target" && -s "$target" ]]; then
    return 0
  fi
  mkdir -p "$(dirname "$target")"
  curl -L --fail -o "$target" "$url"
}

ensure_input_index() {
  if [[ ! -f "${INPUT_VCF_GZ}.tbi" || "$INPUT_VCF_GZ" -nt "${INPUT_VCF_GZ}.tbi" ]]; then
    tabix -f -p vcf "$INPUT_VCF_GZ"
  fi
}

ensure_cache() {
  local cache_version="$1"
  local cache_root="$VEP_DATA_DIR/cache"
  local cache_info="$cache_root/homo_sapiens/${cache_version}_GRCh38/info.txt"
  local cache_tarball="$cache_root/tmp/homo_sapiens_vep_${cache_version}_GRCh38.tar.gz"
  local cache_url="https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache/homo_sapiens_vep_${cache_version}_GRCh38.tar.gz"

  if [[ -f "$cache_info" ]]; then
    echo "Using existing VEP cache: $cache_info"
    return 0
  fi

  mkdir -p "$cache_root/tmp"
  echo "Downloading VEP cache to $cache_tarball"
  curl -L -C - --fail -o "$cache_tarball" "$cache_url"

  if ! gzip -t "$cache_tarball"; then
    echo "Downloaded VEP cache is not a valid gzip archive: $cache_tarball"
    exit 1
  fi

  echo "Extracting VEP cache"
  tar -xzf "$cache_tarball" -C "$cache_root"

  if [[ ! -f "$cache_info" ]]; then
    echo "VEP cache extraction completed but info file is missing: $cache_info"
    exit 1
  fi

  if [[ "$KEEP_CACHE_TARBALL" != "1" ]]; then
    rm -f "$cache_tarball"
  fi
}

ensure_clinvar() {
  local clinvar_vcf="$VEP_DATA_DIR/custom/clinvar.vcf.gz"
  local clinvar_tbi="${clinvar_vcf}.tbi"
  download_if_missing "$CLINVAR_URL" "$clinvar_vcf"
  download_if_missing "$CLINVAR_TBI_URL" "$clinvar_tbi"
}

ensure_docker_ready
ensure_input_index
CACHE_VERSION="$(detect_cache_version)"
ensure_cache "$CACHE_VERSION"
ensure_clinvar

mkdir -p "$(dirname "$OUTPUT_VCF_GZ")"

if [[ -f "$OUTPUT_VCF_GZ" && -f "${OUTPUT_VCF_GZ}.tbi" ]]; then
  echo "Using existing annotated VCF: $OUTPUT_VCF_GZ"
  exit 0
fi

if command -v wslpath >/dev/null 2>&1; then
  MOUNT_DIR_WIN="$(wslpath -w "$SOURCE_DIR")"
else
  MOUNT_DIR_WIN="$SOURCE_DIR"
fi
INPUT_IN_CONTAINER="/data/$(basename "$INPUT_VCF_GZ")"
OUTPUT_RELATIVE="${OUTPUT_VCF_GZ#${SOURCE_DIR}/}"
OUTPUT_IN_CONTAINER="/data/${OUTPUT_RELATIVE}"
WARNING_IN_CONTAINER="${OUTPUT_IN_CONTAINER}_warnings.txt"
CACHE_IN_CONTAINER="/data/${VEP_DATA_DIR#${SOURCE_DIR}/}/cache"
CLINVAR_IN_CONTAINER="/data/${VEP_DATA_DIR#${SOURCE_DIR}/}/custom/clinvar.vcf.gz"

echo "Running VEP annotation"
"$DOCKER_BIN" run --rm -v "${MOUNT_DIR_WIN}:/data" "$VEP_IMAGE" \
  vep \
    --cache \
    --offline \
    --assembly GRCh38 \
    --species homo_sapiens \
    --dir_cache "$CACHE_IN_CONTAINER" \
    --input_file "$INPUT_IN_CONTAINER" \
    --output_file "$OUTPUT_IN_CONTAINER" \
    --warning_file "$WARNING_IN_CONTAINER" \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --fork "$VEP_FORKS" \
    --no_stats \
    --symbol \
    --canonical \
    --check_existing \
    --max_af \
    --custom "file=${CLINVAR_IN_CONTAINER},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN"

tabix -f -p vcf "$OUTPUT_VCF_GZ"
echo "Annotated VCF ready: $OUTPUT_VCF_GZ"
