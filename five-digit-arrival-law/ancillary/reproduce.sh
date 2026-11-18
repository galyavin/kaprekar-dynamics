#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
python_bin="$(command -v "${PYTHON:-python3}")"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p results

check_lean() {
    (
        cd lean
        lake --version
        lean --version
        lake build
        lake env leanchecker --fresh FiveDigit
    ) 2>&1 | tee results/lean_replay.log
}

case "${1:-all}" in
    all)
        "$python_bin" run.py all
        "$python_bin" independent-review/z3_independent.py | tee results/z3_independent_rerun.log
        "$python_bin" headline_checks.py
        check_lean
        ;;
    headline) "$python_bin" headline_checks.py ;;
    lean) check_lean ;;
    review-full) "$python_bin" independent-review/replay.py ;;
    *) "$python_bin" run.py "$@" ;;
esac
