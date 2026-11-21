#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
unset PYTHONOPTIMIZE PYTHONPATH
export LC_ALL=C PYTHONHASHSEED=0
python3 --version
python3 -c 'import sys; assert sys.version_info >= (3, 9), "Python 3.9 or later required"'
shasum -a 256 -c SHA256SUMS

cross=inputs/kaprekar-investigation/cross_length
python3 -B "$cross/rc11_reconstruction/cones.py"
python3 -B "$cross/rc11_reconstruction/physical.py"
python3 -B "$cross/grc/check_instances.py"
python3 -B "$cross/grc/check_corner.py"
g8_output=$(mktemp)
trap 'rm -f "$g8_output"' EXIT
python3 -B "$cross/grc_reconstruction/instances.py" | tee "$g8_output"
grep -Fqx '[G8] exact centered exits [(1, (12, 12, 30, 42, 24)), (5, (36, 60, 78, 90, 48)), (9, (60, 108, 126, 138, 72)), (13, (84, 156, 174, 186, 96))]' "$g8_output"
python3 -B "$cross/grc_reconstruction/physical.py"
python3 -B recompute.py
python3 -B check_rows.py
python3 -B check_evidence.py
printf '%s\n' 'PASS: P10 exact certificates, independent checks and evidence totals.'
