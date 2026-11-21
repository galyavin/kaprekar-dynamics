#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "$0")"
if [[ $# -ne 0 ]]; then
    echo "Usage: bash reproduce.sh" >&2
    exit 2
fi

# In a shared project checkout, an optional existing lock serializes the replay.
if [[ -n ${P11_COMPUTE_LOCK:-} ]]; then
    exec python3 - "$P11_COMPUTE_LOCK" <<'PY'
import fcntl
import os
import subprocess
import sys

with open(sys.argv[1], "r+") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    environment = os.environ.copy()
    del environment["P11_COMPUTE_LOCK"]
    raise SystemExit(subprocess.call(["bash", "reproduce.sh"], env=environment))
PY
fi

export PYTHONDONTWRITEBYTECODE=1
export PYTHONHASHSEED=0
export LC_ALL=C
shasum -a 256 -c SHA256SUMS
python3 --version
"${CXX:-c++}" --version
"${CC:-cc}" --version

work=$(mktemp -d "${TMPDIR:-/tmp}/kaprekar-p11-replay.XXXXXXXX")
trap 'rm -rf "$work"' EXIT

check_output() {
    local expected=$1 actual
    shift
    actual="$work/$(basename "$expected")"
    "$@" > "$actual"
    diff -u "$expected" "$actual"
    echo "Reproduced $expected"
}

# The published Task 207 log starts with two machine-specific timing lines.
sed '1,2d' arith/census.txt > "$work/census.expected"
python3 arith/census.py replay arith/certificates.jsonl independent/census_7_1599.json > "$work/census.actual"
diff -u "$work/census.expected" "$work/census.actual"
python3 census_totals.py
check_output arith2/census-check.txt python3 arith2/check_census.py
check_output recon14_arith/verify.txt python3 recon14_arith/verify.py
check_output recon14_lowt/verify.txt python3 recon14_lowt/verify.py
check_output recon14_sieve/replay.txt python3 recon14_sieve/replay.py

"${CXX:-c++}" -std=c++17 -O3 recon14_arith/scan_strict.cpp -o "$work/scan_strict"
check_output recon14_arith/scan_strict.txt "$work/scan_strict" 1000000000

"${CC:-cc}" -std=c11 -O2 independent/idx2scan.c -o "$work/idx2scan"
sed '$d' independent/idx2_strict_1e9.txt > "$work/index2.expected"
"$work/idx2scan" 1000000000 0 > "$work/index2.actual"
diff -u "$work/index2.expected" "$work/index2.actual"
echo "Reproduced independent/idx2_strict_1e9.txt (without the recorded timing line)"
