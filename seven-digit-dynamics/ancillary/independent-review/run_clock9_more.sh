#!/bin/bash
# depth-nine: Yices on the aux query, then z3 on the plain (no auxiliary bounds) query
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
echo "yices clock9 aux:"; /usr/bin/time -p timeout 900 yices-smt2 smt2/clock9_aux.smt2
echo "z3 clock9 plain:"; "$PY" smt_clock9.py plain - 900
