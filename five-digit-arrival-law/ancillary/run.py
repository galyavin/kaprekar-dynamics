"""Reproduce P2 checks from the bundled source snapshot."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
parser = argparse.ArgumentParser()
parser.add_argument('stage', nargs='?', default='all', choices=['all','audit','limit','joint','maximum','overshoot','numbers'])
args = parser.parse_args()
if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
    raise RuntimeError('Run without -O, -OO, or PYTHONOPTIMIZE.')
OUT.mkdir(exist_ok=True)
work = OUT / 'run'
work.mkdir(exist_ok=True)
manifest = {}
for source in sorted((HERE / 'source').iterdir()):
    if source.is_file():
        shutil.copy2(source, work / source.name)
        manifest[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
(OUT/'source-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
stages = {
    'audit': ['five_clock_audit.py','five_joint_audit.py','five_quantitative_audit.py'],
    'limit': ['five_limit_checks.py'],
    'joint': ['five_joint_checks.py'],
    'maximum': ['five_quantitative_checks.py'],
    'overshoot': ['p2_overshoot_checks.py'],
    'numbers': ['p2_numbers.py'],
}
selected = stages if args.stage == 'all' else {args.stage:stages[args.stage]}
import z3
versions = {'python':sys.version, 'z3':z3.get_version_string(),
            'compiler':subprocess.check_output([os.getenv('CXX','c++'),'--version'],text=True).splitlines()[0]}
(OUT/'versions.json').write_text(json.dumps(versions,indent=2)+'\n')
exe = work/'five_holdout'
if args.stage in ('all','limit','joint'):
    prefixes = ('five_limit_holdout_', 'five_joint_holdout_') if args.stage == 'all' else (('five_limit_holdout_',) if args.stage == 'limit' else ('five_joint_holdout_',))
    for prefix in prefixes:
        for old in work.glob(prefix+'*.json'):
            old.unlink()
    subprocess.run([os.getenv('CXX','c++'),'-O2','-std=c++17',str(work/'five_holdout.cpp'),'-o',str(exe)],check=True)
results = []
for group, scripts in selected.items():
    for name in scripts:
        command = [sys.executable, '-u', str(work/name)]
        if name in ('five_limit_checks.py','five_joint_checks.py'):
            command += ['--holdout-executable',str(exe)]
        start = time.monotonic()
        print('Running',name,flush=True)
        with (OUT/(name.removesuffix('.py')+'.log')).open('w') as log:
            proc = subprocess.Popen(command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,text=True)
            for line in proc.stdout:
                log.write(line); log.flush(); print(line,end='',flush=True)
            code=proc.wait()
        row={'file':name,'returncode':code,'seconds':round(time.monotonic()-start,3)}
        results.append(row)
        (OUT/('run-'+args.stage+'.json')).write_text(json.dumps(results,indent=2)+'\n')
        if code: raise SystemExit(code)
print('All selected checks passed.',flush=True)
