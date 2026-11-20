"""Run a Task 081 verification command with frozen read-only inputs."""
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

here = Path(__file__).resolve().parent
access = here.parent
root = access.parents[1]
inputs = {
    'verify.py': access / 'box_lane_table/verify.py',
    'certificate.py': access / 'box_mod12/certificate.py',
    'quotient.py': access / 'box_mod12/quotient.py',
    'control.py': access / 'box_mod12/control.py',
    'routes.json': access / 'box_mod12/routes.json',
    'quotient.json': access / 'box_mod12/quotient.json',
    'channel.py': access / 'channel_continuation/certify.py',
    'endpoint_verification.json': access / 'channel_continuation/endpoint_verification.json',
    'cohort1.py': access / 'cohort1_reconstruction/reconstruct.py',
    'word2.c': root / 'docs/drafts/research-direction/o2l-task056/word2.c',
}

if __name__ == '__main__':
    argv = sys.argv[1:]
    small = argv[0] == '--small'
    entry, *args = argv[small:]
    with tempfile.TemporaryDirectory(prefix='task081-') as name:
        temp = Path(name)
        for filename, path in inputs.items():
            shutil.copy2(path, temp / filename)
        for path in here.glob('*.py'):
            shutil.copy2(path, temp / path.name)
        command = [sys.executable, '-Bu', entry, *args]
        if small:
            subprocess.run(command, cwd=temp, check=True, timeout=120)
        else:
            with (root / '.task-monitor/compute.lock').open('rb') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                subprocess.run(command, cwd=temp, check=True)
