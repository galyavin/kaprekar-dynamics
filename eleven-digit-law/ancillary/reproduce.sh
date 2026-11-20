#!/usr/bin/env bash
set -euo pipefail

bundle=$(cd "$(dirname "$0")" && pwd)
python=${PYTHON:-python3}
work=$(mktemp -d "${TMPDIR:-/tmp}/kaprekar-p8.XXXXXX")
trap 'rm -rf "$work"' EXIT
export PYTHONDONTWRITEBYTECODE=1

(cd "$bundle/inputs" && shasum -a 256 --quiet -c "$bundle/SHA256SUMS")
cp -R "$bundle/inputs/." "$work/"
mkdir -p "$work/.task-monitor"
if [[ -n "${P8_COMPUTE_LOCK:-}" ]]; then
    ln -s "$P8_COMPUTE_LOCK" "$work/.task-monitor/compute.lock"
else
    touch "$work/.task-monitor/compute.lock"
fi
cp "$bundle/LARGE.sha256" "$work/"
cd "$work"
"$python" -u "$bundle/replay.py" "$@"
