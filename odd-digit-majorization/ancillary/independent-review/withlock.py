"""Run a command while holding the project's exclusive fcntl compute lock."""
import fcntl, os, subprocess, sys, time
LOCK = os.environ.get("KAPREKAR_COMPUTE_LOCK", os.path.join(os.path.dirname(__file__), "compute.lock"))
with open(LOCK, "a") as f:
    t = time.time()
    fcntl.flock(f, fcntl.LOCK_EX)
    sys.stderr.write(f"[lock acquired after {time.time()-t:.1f}s]\n")
    rc = subprocess.call(sys.argv[1:])
    fcntl.flock(f, fcntl.LOCK_UN)
sys.exit(rc)
