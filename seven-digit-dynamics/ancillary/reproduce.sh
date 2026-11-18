#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python_bin="${PYTHON:-python3}"
"$python_bin" - <<'PY'
from pathlib import Path
import hashlib,json,sys
if not __debug__:
    sys.exit('Disable optimized Python for certificate replay.')
root = Path.cwd()
for name, digest in json.loads((root/'results/artifact_sha256.json').read_text()).items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
        sys.exit('Checksum mismatch: ' + name)
PY
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
if ! "$python_bin" -c 'import cvc5, sympy, z3' >/dev/null 2>&1; then
  "$python_bin" -m venv "$work/venv"
  python_bin="$work/venv/bin/python"
  "$python_bin" -m pip install --disable-pip-version-check -r requirements.txt
fi
command -v c++ >/dev/null
command -v yices-smt2 >/dev/null
"$python_bin" -c 'import sys,cvc5,sympy,z3; print("Python",sys.version.split()[0],"Z3",z3.get_version_string(),"cvc5",cvc5.__version__,"SymPy",sympy.__version__)'
yices-smt2 --version
cp -R source "$work/"
for check in seven_clock_audit seven_entry_audit seven_first_access seven_zero_transience_audit seven_boundary_audit seven_facet_checks seven_return_checks seven_tables seven_global_checks seven_clock_checks
do
  echo "Running $check"
  "$python_bin" -u "$work/source/$check.py"
done
echo "Running independent transition model"
"$python_bin" -B independent-review/model.py
echo "All seven-digit checks passed."
