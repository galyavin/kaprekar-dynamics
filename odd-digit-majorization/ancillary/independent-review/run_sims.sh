#!/usr/bin/env bash
set -euo pipefail
# Fresh-base simulations (run under the compute lock)
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
echo "=== route A0,B1,AB table: B=120010+beta (n=200000) and 290015+beta (n=100000), R=0.001B"
"$PY" sim_route.py 200000 0.001 120010 120011 120012 120013 120014
"$PY" sim_route.py 100000 0.001 290015 290016 290017 290018 290019
echo "=== uniform strings: two bases per class, same residue mod 5, different parity and mod 3"
"$PY" sim_uniform.py 160000 60010 60011 60012 60013 60014 60035 60036 60037 60038 60039
echo "=== exact max_N T_B"
for B in 40 41 42 43 44 45 46 47 48 49 50 64 97 128 150 173 199 223 251; do ./kap9 maxT $B; done
