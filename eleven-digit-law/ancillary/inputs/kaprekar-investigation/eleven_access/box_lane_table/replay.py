"""Run the checker against a read-only input snapshot and the shared compute lock."""
import fcntl
from hashlib import sha256
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = ACCESS.parents[1]
FROZEN_SHA256 = '7fe097dfd975e3cb5ab5e90b4d6ff1b3bbab4e122ad066d8a7e99a16a43cc0ec'
assert sha256((HERE / 'frozen-test.md').read_bytes()).hexdigest() == FROZEN_SHA256
INPUTS = {
    'verify.py': HERE / 'verify.py',
    'certificate.py': ACCESS / 'box_mod12/certificate.py',
    'quotient.py': ACCESS / 'box_mod12/quotient.py',
    'control.py': ACCESS / 'box_mod12/control.py',
    'routes.json': ACCESS / 'box_mod12/routes.json',
    'quotient.json': ACCESS / 'box_mod12/quotient.json',
    'channel.py': ACCESS / 'channel_continuation/certify.py',
    'endpoint_verification.json': ACCESS / 'channel_continuation/endpoint_verification.json',
    'cohort1.py': ACCESS / 'cohort1_reconstruction/reconstruct.py',
}

with tempfile.TemporaryDirectory(prefix='task066-') as temp:
    temp = Path(temp)
    for name, path in INPUTS.items():
        shutil.copy2(path, temp / name)
    with (ROOT / '.task-monitor/compute.lock').open('rb') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        subprocess.run([sys.executable, '-B', 'verify.py', *sys.argv[1:]], cwd=temp, check=True)
