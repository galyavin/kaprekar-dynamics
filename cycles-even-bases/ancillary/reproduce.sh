#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
python_bin=${PYTHON:-python3}
compiler=${CXX:-c++}

"$python_bin" - <<'PY'
import hashlib, json
from pathlib import Path
manifest = json.loads(Path("results/artifact_sha256.json").read_text())
for name, digest in manifest.items():
    if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
        raise SystemExit("SHA-256 mismatch in " + name)
print("Verified SHA-256 for %d files." % len(manifest))
PY

scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT
"$python_bin" rules.py > "$scratch/boundary.json"
cmp results/boundary.json "$scratch/boundary.json"
echo 'Passed five digit examples, the full timed boundary table and all critical-state rules at 120 bases.'
"$compiler" -O2 -std=c++17 cycles.cpp -o "$scratch/cycles"
"$scratch/cycles" > "$scratch/cycles.json"
cmp results/cycles.json "$scratch/cycles.json"
echo 'Passed cycle classification at all 500 even bases through 1000 and Table 1 through exponent 12.'
