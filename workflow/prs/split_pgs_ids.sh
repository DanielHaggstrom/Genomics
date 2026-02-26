#!/usr/bin/env bash
set -euo pipefail

# Split a PGS ID list into smaller batch files.
#
# Usage:
#   bash workflow/prs/split_pgs_ids.sh [input_ids_file] [batch_size] [outdir]
#
# Defaults:
#   input_ids_file: config/prs_pgs_ids.txt
#   batch_size: 10
#   outdir: config/prs_batches

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

INPUT_FILE="${1:-config/prs_pgs_ids.txt}"
BATCH_SIZE="${2:-10}"
OUTDIR="${3:-config/prs_batches}"

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Missing input file: $INPUT_FILE"
  exit 1
fi
if ! [[ "$BATCH_SIZE" =~ ^[0-9]+$ ]] || [[ "$BATCH_SIZE" -lt 1 ]]; then
  echo "batch_size must be a positive integer (got: $BATCH_SIZE)"
  exit 1
fi

mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/pgs_ids_batch_*.txt

mapfile -t IDS < <(
  tr -d '\r' < "$INPUT_FILE" \
    | grep -v '^\s*#' \
    | grep -v '^\s*$'
)

TOTAL="${#IDS[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
  echo "No IDs found in $INPUT_FILE"
  exit 1
fi

batch_num=1
for ((i = 0; i < TOTAL; i += BATCH_SIZE)); do
  out_file="$OUTDIR/pgs_ids_batch_$(printf '%02d' "$batch_num").txt"
  : > "$out_file"
  for ((j = i; j < i + BATCH_SIZE && j < TOTAL; j++)); do
    echo "${IDS[j]}" >> "$out_file"
  done
  echo "Wrote $out_file"
  ((batch_num += 1))
done

echo "Done. Total IDs: $TOTAL"
echo "Batches: $((batch_num - 1))"
