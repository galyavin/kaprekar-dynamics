"""Fresh-base campaign runner: runs the simulation binary on the listed bases in parallel and writes logs under an exclusive lock.
Usage: python3 campaign.py SIM_BINARY LOGDIR
"""
import fcntl
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# (bases, samples per base); residues mod 10 cover both parities and all classes mod 5
PLAN = (
    [(b, 20000) for b in range(3534, 3544)]
    + [(b, 20000) for b in range(12340, 12350)]
    + [(b, 5000) for b in range(123450, 123460)]
    + [(b, 1000) for b in range(1234560, 1234570)]
    + [(b, 100) for b in range(12345670, 12345680)]
    + [(b, 12) for b in (123456780, 123456781, 123456782, 123456783, 123456784)]
)


def run(job):
    sim, logdir, B, n = job
    out = os.path.join(logdir, f"sim_B{B}.log")
    t0 = time.time()
    # split big jobs into chunks with distinct seeds
    chunks = max(1, n // 2500)
    with open(out, "w") as f:
        for c in range(chunks):
            m = n // chunks + (1 if c < n % chunks else 0)
            subprocess.run([sim, str(B), str(m), str(1000 + c)], stdout=f, check=True)
    return B, n, time.time() - t0


def main():
    sim, logdir = sys.argv[1], sys.argv[2]
    os.makedirs(logdir, exist_ok=True)
    lock_path = os.environ.get("COMPUTE_LOCK") or os.path.join(logdir, "compute.lock")
    lf = open(lock_path, "a")
    t0 = time.time()
    print("waiting for compute lock", flush=True)
    fcntl.flock(lf, fcntl.LOCK_EX)
    print(f"lock acquired after {time.time() - t0:.0f}s", flush=True)
    try:
        jobs = [(sim, logdir, B, n) for B, n in PLAN]
        # largest first
        jobs.sort(key=lambda j: -j[2] * j[3] if j[2] > 10**6 else -j[2] * j[3] / 10)
        with ThreadPoolExecutor(max_workers=9) as ex:
            for B, n, dt in ex.map(run, jobs):
                print(f"done B={B} n={n} {dt:.1f}s", flush=True)
    finally:
        fcntl.flock(lf, fcntl.LOCK_UN)
        print(f"lock released; total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
