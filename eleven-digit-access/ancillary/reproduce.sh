#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
python_bin=${PYTHON:-python3}
exec "$python_bin" -B "$here/reproduce.py" "$@"
