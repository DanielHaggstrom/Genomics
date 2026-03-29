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

mkdir -p "$ROOT_DIR/output/logs" "$ROOT_DIR/output/nextflow_state/home" "$ROOT_DIR/output/nextflow_launcher"
export NXF_HOME="$ROOT_DIR/output/nextflow_state/home"

cd "$ROOT_DIR/output/nextflow_launcher"

nextflow -log "$ROOT_DIR/output/logs/nextflow.log" run "$ROOT_DIR" -params-file "$ABS_MANIFEST" -profile "$PROFILE"
