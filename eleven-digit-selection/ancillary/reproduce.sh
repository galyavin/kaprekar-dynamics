#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bundle=$PWD
scratch=$(mktemp -d "${TMPDIR:-/tmp}/kaprekar-p7.XXXXXXXX")
trap 'rm -rf "$scratch"' EXIT
export PYTHONDONTWRITEBYTECODE=1
python_cmd=${PYTHON:-python3}
"$python_cmd" --version

"$python_cmd" - "$bundle" "$scratch" "${KAPREKAR_COMPUTE_LOCK:-$scratch/compute.lock}" <<'PY'
import fcntl
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

bundle, scratch, lock_path = map(Path, sys.argv[1:])
python = sys.executable
leaf_digest = '1c62f809a5019babd80d429ec50f5a743fb9531ea064fee519c05c1bd33df63b'
frozen_cohort_digest = 'f9a075479e3f4b0e6f07e0330805ced5e09cc6d156f900e67b3c538037c69e7d'


def run(label, command, expected=None, cwd=scratch, env=None):
    print(label, flush=True)
    start = time.monotonic()
    if expected is None:
        subprocess.run(command, cwd=cwd, env=env, check=True)
    else:
        output = subprocess.check_output(command, cwd=cwd, env=env)
        if output != (expected if isinstance(expected, bytes) else expected.read_bytes()):
            raise SystemExit(f'Output differs from the frozen record for {label}')
    print(f'{label} passed in {time.monotonic() - start:.2f} s', flush=True)


with lock_path.open('a') as lock:
    wait = time.monotonic()
    fcntl.flock(lock, fcntl.LOCK_EX)
    print(f'Compute lock acquired after {time.monotonic() - wait:.2f} s', flush=True)
    start = time.monotonic()
    primary = bundle / 'parity'
    audit = bundle / 'audit'
    cohort = primary / 'reconstruction.txt'
    if hashlib.sha256(cohort.read_bytes()).hexdigest() != frozen_cohort_digest:
        raise SystemExit('Frozen cohort digest mismatch')
    run('Cohort lattice', [python, primary / 'cohort.py'])
    run('Whole-word witness', [python, primary / 'exact_control.py'],
        primary / 'exact_control.txt')
    run('Independent guarded cohort replay', [python, audit / 'cohort_replay.py', cohort])
    dump = scratch / 'leaves.jsonl.gz'
    run('Even exclusion tree and leaf dump',
        [python, bundle / 'dump_tree.py', primary, dump], primary / 'exclusion.txt')
    digest = hashlib.sha256()
    with gzip.open(dump, 'rb') as leaves:
        for block in iter(lambda: leaves.read(1 << 20), b''):
            digest.update(block)
    if digest.hexdigest() != leaf_digest:
        raise SystemExit(f'Leaf digest mismatch: {digest.hexdigest()}')
    print(f'Uncompressed terminal leaf SHA-256 {leaf_digest}', flush=True)
    # The supplied audit takes a lock itself. The parent holds the shared
    # lock for the whole run, so give the audit a distinct private lock.
    audit_lock = scratch / 'audit.lock'
    audit_lock.touch()
    run('Independent full terminal-partition audit',
        [python, audit / 'audit.py', dump, audit_lock])
    mixture = bundle / 'mixture'
    cxx = os.environ.get('CXX', 'clang++')
    subprocess.run([cxx, '--version'], check=True)
    verifier = scratch / 'mixture_verify'
    run('Build six-residue route verifier',
        [cxx, '-O3', '-std=c++17', '-Wall', '-Wextra', mixture / 'verify.cpp', '-o', verifier])
    tree = mixture / 'tree.txt'
    header, *record_lines, footer = tree.read_bytes().splitlines(keepends=True)
    if not (header.startswith(b'# lock acquired after ') and footer.startswith(b'# exit 0, wall ')):
        raise SystemExit('Unexpected provenance in the frozen route record')
    records = b''.join(record_lines)
    run('Six-residue route verification', [verifier], records)
    tree_audit = subprocess.check_output([python, mixture / 'tree_check.py', tree], cwd=scratch)
    if not tree_audit.rstrip().endswith(b'ALL TREE SUMMARY CHECKS PASS'):
        raise SystemExit(tree_audit.decode())
    print('Exact lane clock and weight audit passed', flush=True)
    run('Independent mixture arithmetic', [python, mixture / 'mixture.py', 'tables'],
        mixture / 'tables.log')
    sample_dump = scratch / 'mix-leaves-0.txt.gz'
    dump_env = dict(os.environ, KAPREKAR_COMPUTE_LOCK=str(scratch / 'dump.lock'))
    beta_zero = records.split(b'beta 1 ', 1)[0]
    run('Independent leaf-dump smoke check for residue 0',
        [python, mixture / 'dump_leaves.py', '0', '0', sample_dump], beta_zero, env=dump_env)
    with gzip.open(sample_dump, 'rb') as sample:
        if sample.read() != (mixture / 'sample-leaves-0.txt').read_bytes():
            raise SystemExit('Leaf-dump sample differs from its frozen record')

    box = scratch / 'box'
    shutil.copytree(bundle / 'box', box)
    frozen = '0c857f4cc571cc24be00d260d0eda53ba4380bc3bb8d149f08c5c56631d1c893'
    if hashlib.sha256((box / 'frozen-test.md').read_bytes()).hexdigest() != frozen:
        raise SystemExit('Frozen second-box test digest mismatch')
    run('Box prefix for cohort 1', [python, 'channel.py'], box / 'endpoint_verification.json', cwd=box)
    run('Box prefix for cohorts 0 and 2', [python, 'certificate.py'], box / 'routes.json', cwd=box)
    rebuild_cohort_one = """import json
from pathlib import Path
import cohort1

records = json.loads(Path('endpoint_verification.json').read_text())['residues']
_, result = cohort1.enumerate_quotient(cohort1.build_sections(records))
print(json.dumps(result, default=str, indent=2, sort_keys=True))
"""
    run('Box affine matching for cohort 1', [python, '-c', rebuild_cohort_one],
        box / 'accepted_quotient.json', cwd=box)
    run('Box whole-population affine matching', [python, 'quotient.py'], box / 'quotient.json', cwd=box)
    run('Box residual routes and plateaus', [python, 'verify.py'], box / 'result.json', cwd=box)
    print(f'Checks through residual routes passed in {time.monotonic() - start:.2f} s', flush=True)
    run('Twelve-shift structural sections', [python, 'structural.py'], box / 'structure.json', cwd=box)
    run('Guarded coverage of the second box', [python, 'coverage.py'], box / 'coverage.json', cwd=box)
    run('Exact pullbacks and limit weights', [python, 'limit_weights.py', 'coverage.json'],
        box / 'limit_weights.json', cwd=box)
    run('Twelve-shift locality controls', [python, 'locality.py'], box / 'locality.json', cwd=box)
    run('Marked clock weights', [python, 'marked_weights.py'], box / 'marked-weights.json', cwd=box)
    run('Marked clock controls', [python, 'marked_clock.py'], box / 'marked-clock.json', cwd=box)
    run('Separate twelve-shift algebra', [python, bundle / 'independent-review/algebra.py'],
        bundle / 'independent-review/results-algebra.txt')
    print(f'All checks passed in {time.monotonic() - start:.2f} s under the lock', flush=True)
PY
