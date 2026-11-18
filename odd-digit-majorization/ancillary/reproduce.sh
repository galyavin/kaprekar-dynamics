#!/usr/bin/env bash
set -eu
bundle_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
interpreter="${PYTHON:-$bundle_dir/.venv/bin/python}"
if [[ -z "${PYTHON:-}" && ! -x "$interpreter" ]]; then
  interpreter=python3
fi
exec "$interpreter" "$bundle_dir/run.py" "$@"
