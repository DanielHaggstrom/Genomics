#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

PATTERN='(AKIA[0-9A-Z]{16}|BEGIN (RSA|EC|OPENSSH|PGP) PRIVATE KEY|api[_-]?key\s*=|token\s*=|secret\s*=|password\s*=|/Users/|C:/Users/|/mnt/c/Users/)'
TRACKED_FILES=()
while IFS= read -r f; do
  case "$f" in
    *.example*|.env.example|workflow/utils/prepublish_check.sh)
      continue
      ;;
    *)
      TRACKED_FILES+=("$f")
      ;;
  esac
done < <(git ls-files)

if command -v rg >/dev/null 2>&1; then
  set +e
  MATCHES="$(rg -n --pcre2 -i "$PATTERN" "${TRACKED_FILES[@]}" || true)"
  set -e
else
  set +e
  MATCHES="$(grep -nEi "$PATTERN" "${TRACKED_FILES[@]}" 2>/dev/null || true)"
  set -e
fi

if [[ -n "$MATCHES" ]]; then
  echo "Potential secret/private-path leaks found in tracked files:"
  echo "$MATCHES"
  exit 1
fi

echo "Prepublish check passed: no obvious secret/private-path leaks in tracked files."
