#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

DRY_RUN=0

usage() {
  cat <<EOF
Usage: bash workflow/utils/cleanup_runs.sh [--dry-run]

Deletes ignored, generated pipeline artifacts from the repository root:
  - results/
  - logs/
  - work/
  - .nextflow/
  - .nextflow.log*
  - trace.txt, timeline.html, report.html, dag.*
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
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

remove_target() {
  local target="$1"
  local resolved="$ROOT_DIR/${target#./}"

  case "$resolved" in
    "$ROOT_DIR" | "$ROOT_DIR"/*)
      ;;
    *)
      echo "Refusing to remove path outside repo root: $target"
      exit 1
      ;;
  esac

  if [[ ! -e "$target" && ! -L "$target" ]]; then
    return 0
  fi

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "Would remove: $target"
    return 0
  fi

  rm -rf -- "$target"
  echo "Removed: $target"
}

shopt -s nullglob
TARGETS=(
  ".nextflow"
  "logs"
  "results"
  "work"
  ".nextflow.log"*
  "trace.txt"
  "timeline.html"
  "report.html"
  "dag."*
)
shopt -u nullglob

for target in "${TARGETS[@]}"; do
  remove_target "$target"
done

if [[ "$DRY_RUN" == "1" ]]; then
  echo "Dry run complete."
else
  echo "Cleanup complete."
fi
