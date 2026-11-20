"""Run the frozen experiments in a private snapshot under the compute lock."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import fcntl
import shutil
import subprocess
import tempfile


def run(directory, executable, base, size, seed, mode):
    path = directory / f"{executable}-{mode}-{base}.tsv"
    with path.open("w") as out:
        subprocess.run([f"./{executable}", str(base), str(size), str(seed), mode],
                       cwd=directory, stdout=out, check=True)
    return path.name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--followup", action="store_true")
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    repo = source.parents[2]
    directory = Path(tempfile.mkdtemp(prefix="task097-holdout-")).resolve()
    print(directory, flush=True)
    inputs = ("PREDICTIONS.md", "FREEZE.sha256", "graph.cpp", "actual.cpp")
    for name in inputs:
        shutil.copyfile(source / name, directory / name)
    subprocess.run(["shasum", "-a", "256", "-c", "FREEZE.sha256"],
                   cwd=directory, check=True)
    for binary in ("graph", "actual"):
        subprocess.run(["clang++", "-O3", "-std=c++17", "-Wall", "-Wextra",
                        f"{binary}.cpp", "-o", binary], cwd=directory, check=True)
    jobs = []
    if args.followup:
        for name in ("FOLLOWUP.md", "FOLLOWUP.sha256"):
            shutil.copyfile(source / name, directory / name)
        subprocess.run(["shasum", "-a", "256", "-c", "FOLLOWUP.sha256"],
                       cwd=directory, check=True)
        for base in (48001, 48007):
            jobs.extend([
                ("actual", base, 5000, 970042031, "uniform-pi"),
                ("graph", base, 5000, 970042039, "uniform-pi"),
            ])
    else:
        for start in (12000, 48000):
            for r in range(12):
                base = start+r
                jobs.extend([
                    ("actual", base, 8000, 970041003, "original"),
                    ("graph", base, 3000, 970041007, "entry"),
                    ("actual", base, 1000, 970041013, "uniform-pi"),
                    ("graph", base, 1000, 970041021, "uniform-pi"),
                ])
        for base in (1200005, 1200011):
            jobs.append(("actual", base, 24, 970041031, "box"))
    with (repo / ".task-monitor" / "compute.lock").open("r+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        print("acquired compute lock", flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run, directory, *job) for job in jobs]
            for future in futures:
                print(future.result(), flush=True)


if __name__ == "__main__":
    main()

