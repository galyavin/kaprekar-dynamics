#!/usr/bin/env bash
set -euo pipefail
# Replicate the four family-B kernels with |z|>3: same bases with a fresh seed, and bases above 10^6.
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
"$PY" -B sim_targeted.py targeted_sameB.npz 777 12000 "A2/A:A0/AA|2|451067" "A0/A:A0/AA|4|451009" "B1/B:A0/A:A1/AA|3|451048" "B1/B:A0/A:A0/AA|3|451038" > logs_targeted.txt 2>&1
"$PY" -B sim_targeted.py targeted_bigB.npz 778 12000 "A2/A:A0/AA|2|1200042" "A0/A:A0/AA|4|1200019" "B1/B:A0/A:A1/AA|3|1200033" "B1/B:A0/A:A0/AA|3|1200043" >> logs_targeted.txt 2>&1
"$PY" -B -W ignore analyze_allroutes.py targeted_sameB.npz targeted_bigB.npz >> logs_targeted.txt 2>&1
mv allroutes_analysis.json targeted_analysis.json
