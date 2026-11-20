"""Run Task 139 from a private snapshot under the shared compute lock."""

import argparse
import fcntl
from pathlib import Path
from shutil import copyfile
from subprocess import run
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def invoke(cwd, script, *args):
    run((sys.executable, "-B", script, *args), cwd=cwd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", type=Path, metavar="OUTPUT",
                        help="write a fresh gzip certificate after a private run")
    args = parser.parse_args()
    inputs = [
        (HERE.parent / "e1_basins/catalogue_checks.py", "e1_basins/catalogue_checks.py"),
        (HERE.parent / "global_work/evidence/energy.json.gz", "energy.json.gz"),
        (HERE.parent / "entry_pi/catalogue.json", "entry_pi/catalogue.json"),
        (ROOT / "docs/drafts/research-direction/astra-strategy-13/basins.out", "basins.out"),
        *((HERE.parent / p, p) for p in (
            "en11b/entrance.py", "en11b/clipping.py", "stage_graph/pi_geometry.py",
            "sb_exhaustion/geometry.py",
        )),
        (ROOT / "docs/drafts/research-direction/astra-strategy-12/one_step.py",
         "docs/drafts/research-direction/astra-strategy-12/one_step.py"),
        *((p, f"e1b_guarded/{p.name}") for p in HERE.glob("*.py")),
    ]
    inputs = [(source, target if target.startswith("docs/") else f"kaprekar-investigation/eleven_access/{target}")
              for source, target in inputs]
    if args.search is None:
        inputs.append((HERE / "certificate.json.gz",
                       "kaprekar-investigation/eleven_access/e1b_guarded/certificate.json.gz"))
    with TemporaryDirectory(prefix="e1b-") as directory:
        temp = Path(directory)
        for source, target in inputs:
            path = temp / target
            path.parent.mkdir(parents=True, exist_ok=True)
            copyfile(source, path)
        with (ROOT / ".task-monitor/compute.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            access = temp / "kaprekar-investigation/eleven_access"
            invoke(access, "e1_basins/catalogue_checks.py", "--catalogue", "entry_pi/catalogue.json",
                   "--energy", "energy.json.gz", "--astra", "basins.out", "--diagnose")
            target = access / "e1b_guarded"
            if args.search is None:
                invoke(target, "sweep.py", "--replay", "certificate.json.gz")
            else:
                invoke(target, "sweep.py", "--output", "certificate.json.gz")
            invoke(target, "row181.py", "certificate.json.gz")
            invoke(target, "mass.py")
        if args.search is not None:
            copyfile(target / "certificate.json.gz", args.search)


if __name__ == "__main__":
    main()
