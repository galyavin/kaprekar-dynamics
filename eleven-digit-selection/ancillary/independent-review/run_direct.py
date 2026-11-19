"""Rerun the exact original-word sampling and direct subtraction checks."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
BOX_CASES = (
    (4200185, 5000, 100, 6221953),
    (4800233, 5000, 100, 747529),
    (4300175, 5000, 100, 890241),
    (5000159, 5000, 100, 533093),
    (87234569, 16, 16, 956643),
    (87234575, 16, 16, 643850),
)
U_CASES = ((2934720, 1000, 521959), (2934721, 1000, 714183))


def run(binary, arguments):
    start = time.monotonic()
    output = subprocess.check_output([binary, *map(str, arguments)], text=True)
    return output + f"wall_seconds {time.monotonic() - start:.2f}\n"


def main():
    with tempfile.TemporaryDirectory(prefix="kaprekar-direct-") as temp:
        binary = Path(temp) / "direct"
        compiler = os.environ.get("CXX", "clang++")
        subprocess.run([compiler, "-std=c++17", "-O3", "-Wall", "-Wextra", HERE / "direct.cpp",
                        "-o", binary], check=True)
        lock_path = Path(os.environ.get("KAPREKAR_COMPUTE_LOCK", str(Path(temp) / "compute.lock")))
        with lock_path.open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            start = time.monotonic()
            print("Compute lock acquired", flush=True)
            cases = [("box", *case) for case in BOX_CASES] + [("u", *case) for case in U_CASES]
            for mode, base, *_ in cases:
                print(f"{mode} B={base}, mod 12={base % 12}", flush=True)
            print(subprocess.check_output([sys.executable, HERE / "constants.py", binary], text=True), flush=True)
            workers = min(len(cases), os.cpu_count() or 1)
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for output in pool.map(lambda args: run(binary, args), cases):
                    print(output, flush=True)
            print(f"total_wall_seconds {time.monotonic() - start:.2f}", flush=True)


if __name__ == "__main__":
    main()
