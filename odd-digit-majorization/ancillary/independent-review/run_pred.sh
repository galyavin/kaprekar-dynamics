#!/usr/bin/env bash
set -euo pipefail
# Predicted global statistics from the supplement's kernels, 8 x 6000 samples (under compute lock)
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
pids=()
for s in 1 2 3 4 5 6 7 8; do
  "$PY" -B pred_global.py 6000 "$s" "pred_part$s.json" > "logs_pred_part$s.txt" 2>&1 &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
exit "$status"
