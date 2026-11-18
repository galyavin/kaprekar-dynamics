#!/usr/bin/env bash
set -euo pipefail
# All 165 kernels at fresh bases, two base families of opposite parity (run under compute lock).
# Raw samples are saved; analysis (with error-bounded atom grouping) is analyze_allroutes.py.
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
"$PY" -B sim_allroutes.py 3000 60 300000 allroutes_raw_300k.npz > logs_allroutes_300k.txt 2>&1
"$PY" -B sim_allroutes.py 3000 60 450000 allroutes_raw_450k.npz 5 > logs_allroutes_450k.txt 2>&1
"$PY" -B analyze_allroutes.py allroutes_raw_300k.npz allroutes_raw_450k.npz > logs_allroutes_analysis.txt 2>&1
mv allroutes_analysis.json allroutes_analysis_families.json
