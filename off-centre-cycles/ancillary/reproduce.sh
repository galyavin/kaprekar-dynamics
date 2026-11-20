#!/usr/bin/env bash
set -euo pipefail

bundle=$(cd "$(dirname "$0")" && pwd)
python=${PYTHON:-python3}
work=$(mktemp -d "${TMPDIR:-/tmp}/kaprekar-odd.XXXXXXXX")
trap 'rm -rf "$work"' EXIT

(cd "$bundle" && shasum -a 256 --quiet -c SHA256SUMS)
cp -R "$bundle/inputs/." "$work/"
cp "$bundle/check_primary.py" "$work/"
cd "$work"
export PYTHONDONTWRITEBYTECODE=1
"$python" -I -B -u "$bundle/replay.py" "$work" "${P9_COMPUTE_LOCK:-$work/compute.lock}"
