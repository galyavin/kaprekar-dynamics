#!/usr/bin/env python3
"""Run a shell command while holding the exclusive fcntl compute lock.
Usage: locked.py LOCKFILE LOGFILE -- cmd ...
"""
import fcntl
import subprocess
import sys
import time

lock, log = sys.argv[1], sys.argv[2]
cmd = sys.argv[sys.argv.index("--") + 1:]
with open(lock, "a") as f:
    t = time.time()
    fcntl.flock(f, fcntl.LOCK_EX)
    with open(log, "a") as out:
        out.write(f"# lock acquired after {time.time() - t:.1f}s; cmd: {' '.join(cmd)}\n")
        out.flush()
        t = time.time()
        r = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT)
        out.write(f"# exit {r.returncode}, wall {time.time() - t:.1f}s\n")
    fcntl.flock(f, fcntl.LOCK_UN)
sys.exit(r.returncode)
