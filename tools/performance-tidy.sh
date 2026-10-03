#!/bin/sh
# Advisory performance audit; source files are never fixed automatically.
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "${PYTHON:-python3}" "$SCRIPT_DIR/performance_tidy.py" "$@"
