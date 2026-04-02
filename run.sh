#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 ]]; then
  echo "Usage: bash run.sh manifests/<sample>.yaml [profile] [nextflow args...]"
  exit 1
fi

MANIFEST="$1"
shift
PROFILE="docker"
if [[ $# -gt 0 && "$1" != -* ]]; then
  PROFILE="$1"
  shift
fi
NEXTFLOW_ARGS=("$@")
ABS_MANIFEST="$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")"

PROJECT_CACHE_ROOT="$(
  PYTHONPATH="$ROOT_DIR" python3 -m workflow.core.manifest print-project-cache-root \
    --manifest "$ABS_MANIFEST"
)"

mkdir -p "$PROJECT_CACHE_ROOT/nextflow/home" "$PROJECT_CACHE_ROOT/nextflow/launcher" "$PROJECT_CACHE_ROOT/nextflow/work"
export NXF_HOME="$PROJECT_CACHE_ROOT/nextflow/home"

cd "$PROJECT_CACHE_ROOT/nextflow/launcher"

nextflow \
  -log "$PROJECT_CACHE_ROOT/nextflow/nextflow.log" \
  run "$ROOT_DIR" \
  -params-file "$ABS_MANIFEST" \
  -profile "$PROFILE" \
  -work-dir "$PROJECT_CACHE_ROOT/nextflow/work" \
  "${NEXTFLOW_ARGS[@]}"
