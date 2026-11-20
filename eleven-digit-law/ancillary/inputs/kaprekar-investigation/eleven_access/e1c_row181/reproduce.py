"""Replay Task 139 and the new row 181 bound in locked private snapshots."""

import argparse
import fcntl
from pathlib import Path
from shutil import copyfile
from subprocess import run
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = HERE.parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound-only", action="store_true", help="skip the Task 139 replay")
    args = parser.parse_args()
    if not args.bound_only:
        run([sys.executable, "-B", ACCESS / "e1b_guarded/reproduce.py"], cwd=ROOT, check=True)
    inputs = [
        HERE / "mass_bound.py",
        ACCESS / "entry_pi/catalogue.json",
        ACCESS / "e1_basins/catalogue_checks.py",
        *(ACCESS / name for name in (
            "en11b/entrance.py", "en11b/clipping.py", "en11b/density.cpp",
            "e1b_guarded/anchor.py", "e1b_guarded/guarded.py", "e1b_guarded/mass.py",
            "stage_graph/pi_geometry.py", "sb_exhaustion/geometry.py",
        )),
        ROOT / "docs/drafts/research-direction/astra-strategy-12/one_step.py",
    ]
    with TemporaryDirectory(prefix="e1c-") as directory:
        snapshot = Path(directory)
        for source in inputs:
            target = snapshot / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            copyfile(source, target)
        here = snapshot / HERE.relative_to(ROOT)
        density = here / "density"
        with (ROOT / ".task-monitor/compute.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            run(["clang++", "-O2", "-std=c++17", here.parent / "en11b/density.cpp", "-o", density],
                cwd=here, check=True)
            run([sys.executable, "-B", "mass_bound.py", density], cwd=here, check=True)


if __name__ == "__main__":
    main()
