#!/usr/bin/env bash
set -euo pipefail
# Exhaustive check of eq. (descent) at every base 101..160 and at 173,191,199,211,223 (run under compute lock)
cd "$(dirname "$0")"
( for B in $(seq 101 160) 173 191 199 211 223; do echo $B; done ) | xargs -P 6 -I{} ./kap9 descent {} | sort -t= -k2 -n
./kap9 descrand 2000000 11
./kap9 descrand 2000000 12
