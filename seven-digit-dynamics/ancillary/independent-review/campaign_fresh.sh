#!/bin/bash
# Fresh-base campaign. Run in a working copy after compiling exhaust.c and traj.c.
set -eu
pids=()
wait_batch() {
  local pid failed=0
  for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
  pids=()
  return "$failed"
}
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
mkdir -p logs
date
# (a) exhaustive lemma checks at every base 18..300, and reproduction at the paper's bases
( ./exhaust 18 300 1 0 1 > logs/exhaust_18_300.txt ) & pids+=("$!")
( for B in 129 130 400 401; do ./exhaust $B $B 1 0 1; done > logs/exhaust_paperbases.txt ) & pids+=("$!")
# (b) fresh bases, exact weighted law over all B^7-B words
( for B in 500 501 502 503; do ./exhaust $B $B 1 1 1; done > logs/exhaust_500_503.txt ) & pids+=("$!")
( for B in 504 505 506 507; do ./exhaust $B $B 1 1 1; done > logs/exhaust_504_507.txt ) & pids+=("$!")
wait_batch
( ./exhaust 1000 1000 1 1 1 > logs/exhaust_1000.txt ) & pids+=("$!")
( ./exhaust 1001 1001 1 1 1 > logs/exhaust_1001.txt ) & pids+=("$!")
# (c) two clocks / parity weights: boxes of 512 residues at fresh large bases
for k in 0 1 2 3 4 5 6 7; do B=$((200000+k)); ( ./traj clocks $B 100 $((1000+k)) > logs/clocks_$B.txt ) & pids+=("$!"); done
wait_batch
for B in 1000003 1000004; do ( ./traj clocks $B 30 $B > logs/clocks_$B.txt ) & pids+=("$!"); done
# (d) population law by uniform random words at fresh bases
for B in 20011 20012 20013 20014; do ( ./traj pop $B 100000 $B > logs/pop_$B.txt ) & pids+=("$!"); done
for B in 100005 100006; do ( ./traj pop $B 20000 $B > logs/pop_$B.txt ) & pids+=("$!"); done
( ./traj limitstats 20000000 99 > logs/limitstats.txt 2>&1 ) & pids+=("$!")
( "$PY" prop_tests.py 12 300 300 > logs/prop_tests.txt 2>&1 ) & pids+=("$!")
wait_batch
date
