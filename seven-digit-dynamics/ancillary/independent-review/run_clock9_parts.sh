#!/bin/bash
# depth-nine claim, own encoding, split by the 16 source sign patterns; Yices and z3 per part.
cd "$(dirname "$0")"
mkdir -p logs
run_part() {
  k=$1
  r=$( { /usr/bin/time -p timeout 600 yices-smt2 smt2/clock9_part$k.smt2; } 2>&1 | tr '\n' ' ')
  echo "part $k yices: $r"
}
export -f run_part
seq 0 15 | xargs -P 8 -I{} bash -c 'run_part {}' | sort -n -k2 > logs/clock9_parts_yices.txt
cat logs/clock9_parts_yices.txt
