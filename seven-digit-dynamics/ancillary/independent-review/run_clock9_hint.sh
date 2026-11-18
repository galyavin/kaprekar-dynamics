#!/bin/bash
# Yices on the nine partitions that timed out, with implied ordering facts added (equivalent queries).
cd "$(dirname "$0")"
run_part() {
  k=$1
  r=$( { /usr/bin/time -p timeout 450 yices-smt2 smt2/clock9_hint_part$k.smt2; } 2>&1 | tr '\n' ' ')
  echo "part $k yices+implied-order: $r"
}
export -f run_part
printf "1\n2\n3\n4\n5\n8\n9\n12\n13\n" | xargs -P 9 -I{} bash -c 'run_part {}' | sort -n -k2 > logs/clock9_hint_yices.txt
cat logs/clock9_hint_yices.txt
