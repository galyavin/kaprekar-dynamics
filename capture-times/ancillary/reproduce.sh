#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python_bin=${PYTHON:-python3}
export LC_ALL=C PYTHONHASHSEED=0
unset PYTHONOPTIMIZE

"$python_bin" -B - <<'PY'
from hashlib import sha256
from pathlib import Path

for line in Path("SHA256SUMS").read_text().splitlines():
    expected, filename = line.split("  ", 1)
    with Path(filename).open("rb") as artifact:
        digest = sha256()
        for block in iter(lambda: artifact.read(1 << 20), b""):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise SystemExit(f"SHA-256 mismatch: {filename}")
print("SHA-256 manifest verified")
PY

work_dir=$(mktemp -d)
trap 'rm -rf -- "$work_dir"' EXIT
printf 'Runtime: '
"$python_bin" --version

replay() {
    local label=$1 expected=$2
    shift 2
    printf 'Checking %s\n' "$label"
    "$python_bin" -B "$@" > "$work_dir/output"
    diff -u "$expected" "$work_dir/output"
    printf 'Verified %s\n' "$label"
}

replay "physical and order 3–13 certificates" recon14_physical/verification.log \
    recon14_physical/verify.py recon14_physical/kaprekar-off-centre-cycles.tex \
    recon14_physical/kaprekar-off-centre-cycles-ancillary.zip
replay "radius, entrance and timing controls" recon15/checks.txt recon15/check.py
replay "period-dividing-2r graph and entry controls" recon_carriers/verification.txt \
    recon_carriers/verify.py
replay "order-15 certificates" recon_root15/verification.json \
    recon_root15/replay.py recon_root15
printf '%s\n' 'All four reconstructions passed.'
