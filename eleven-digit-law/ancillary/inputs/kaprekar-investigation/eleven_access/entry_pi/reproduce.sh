#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "$0")" && pwd)
cd "$here"
for freeze in ./*-FREEZE.sha256; do
    shasum -a 256 -c "$freeze"
done
snapshot=$(mktemp -d "${TMPDIR:-/tmp}/task098-proof.XXXXXXXX")
trap 'rm -rf "$snapshot"' EXIT
cp catalogue.py regular_geometry.py verify.py ../global_work/{faces,lp,work}.py "$snapshot/"
mkdir "$snapshot/evidence"
cp ../global_work/evidence/energy.json.gz "$snapshot/evidence/"
cd "$snapshot"
export PYTHONDONTWRITEBYTECODE=1
python3 catalogue.py > catalogue.log
diff -u "$here/catalogue.json" catalogue.json
python3 verify.py
