"""Run a command while holding the shared exclusive compute lock."""
import fcntl, os, subprocess, sys, tempfile
LOCK = os.environ.get("COMPUTE_LOCK", os.path.join(tempfile.gettempdir(), "p2-compute.lock"))
with open(LOCK, "a") as fd:
    fcntl.flock(fd, fcntl.LOCK_EX)
    try:
        sys.exit(subprocess.call(sys.argv[1:]))
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
