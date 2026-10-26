#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  export PATH="$ROOT_DIR/.venv/bin:$PATH"
fi
export PYTHONPATH="$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}"
exec python3 -m workflow.run "$@"
