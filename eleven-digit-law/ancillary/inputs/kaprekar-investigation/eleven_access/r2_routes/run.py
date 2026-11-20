"""Run a fresh locked diagnostic in a private snapshot and print its summary."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    source = Path(__file__).resolve().parent
    repo = source.parents[2]
    subprocess.run(["shasum", "-a", "256", "-c", "FREEZE.sha256"],
                   cwd=source, stdout=sys.stderr, check=True)
    with tempfile.TemporaryDirectory(prefix="r2-class3-") as folder:
        snapshot = Path(folder)
        for name in ("diagnostic.cpp", "check_physical.cpp", "summarize.py"):
            shutil.copy2(source/name, snapshot/name)
        shutil.copy2(source.parent/"stage_graph/graph.cpp", snapshot/"graph.cpp")
        for name in ("diagnostic", "check_physical"):
            subprocess.run(["clang++", "-std=c++17", "-O3", "-Wall", "-Wextra", "-Wno-return-type",
                            f"{name}.cpp", "-o", name], cwd=snapshot, check=True)

        def shard(begin):
            with (snapshot/f"rows-{begin:05}.txt").open("wb") as output:
                subprocess.run(["./diagnostic", str(begin), "7500"],
                               cwd=snapshot, stdout=output, check=True)
            print(f"finished vectors {begin} to {begin+7499}", file=sys.stderr, flush=True)

        with (repo/".task-monitor/compute.lock").open("rb") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            subprocess.run(["./check_physical"], cwd=snapshot, stdout=sys.stderr, check=True)
            with ThreadPoolExecutor(max_workers=8) as pool:
                list(pool.map(shard, range(0, 60_000, 7500)))
        subprocess.run([sys.executable, "-B", "summarize.py", folder], cwd=snapshot, check=True)


if __name__ == "__main__":
    main()
