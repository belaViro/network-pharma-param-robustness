#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    echo "Error: Python >= 3.10 is required but was not found on PATH." >&2
    exit 1
fi

exec "$PYTHON_BIN" reproduce.py "$@"
