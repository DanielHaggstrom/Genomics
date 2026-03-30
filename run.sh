#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash run.sh manifests/<sample>.yaml [profile]"
  exit 1
fi

MANIFEST="$1"
PROFILE="${2:-docker}"
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
  -work-dir "$PROJECT_CACHE_ROOT/nextflow/work"
