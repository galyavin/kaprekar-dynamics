"""Replay the P8 certificates from a hash-checked private source tree."""
import argparse
from contextlib import contextmanager
import fcntl
import json
from math import isclose
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path.cwd()
ACCESS = ROOT / "kaprekar-investigation/eleven_access"
DRAFT = ROOT / "docs/drafts/research-direction"
AUDIT = DRAFT / "astra-theorem-restricted"
BOX_LEAVES = DRAFT / "replay-task071/results/leaves.jsonl"


def matches_panel_summary(actual, expected):
    """Compare row hashes and counts exactly, allowing bounded float rounding."""
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, float):
        return isclose(actual, expected, rel_tol=0, abs_tol=1e-13)
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            matches_panel_summary(actual[key], value) for key, value in expected.items()
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            matches_panel_summary(left, right) for left, right in zip(actual, expected)
        )
    return actual == expected


def run(command, *, cwd=ROOT, expected=None, output=None, statistical=False):
    if expected is not None:
        data = subprocess.check_output(command, cwd=cwd)
        reference = expected.read_bytes()
        if statistical:
            matches = matches_panel_summary(json.loads(data), json.loads(reference))
        elif expected.suffix == ".json":
            matches = json.loads(data) == json.loads(reference)
        else:
            matches = data == reference
        if not matches:
            raise RuntimeError(f"Replay differs from {expected.relative_to(ROOT)}")
        return
    subprocess.run(command, cwd=cwd, stdout=output, check=True)


def python(script, *args, expected=None, statistical=False):
    run([sys.executable, "-B", script, *args], expected=expected, statistical=statistical)


@contextmanager
def compute_lock():
    with (ROOT / ".task-monitor/compute.lock").open("rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def build(source, executable, *flags):
    run(["clang++", "-std=c++17", "-O3", *flags, source, "-o", executable])


def regenerate_box_leaves():
    files = {
        "verify.py": "box_lane_table/verify.py",
        "certificate.py": "box_mod12/certificate.py",
        "quotient.py": "box_mod12/quotient.py",
        "control.py": "box_mod12/control.py",
        "routes.json": "box_mod12/routes.json",
        "quotient.json": "box_mod12/quotient.json",
        "channel.py": "channel_continuation/certify.py",
        "endpoint_verification.json": "channel_continuation/endpoint_verification.json",
        "cohort1.py": "cohort1_reconstruction/reconstruct.py",
    }
    BOX_LEAVES.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="p8-box-") as directory:
        work = Path(directory)
        for name, source in files.items():
            shutil.copyfile(ACCESS / source, work / name)
        shutil.copyfile(BOX_LEAVES.parents[1] / "leafdump_driver.py", work / "leafdump_driver.py")
        with compute_lock(), BOX_LEAVES.open("wb") as output:
            run([sys.executable, "-B", "leafdump_driver.py"], cwd=work, output=output)
    run(["shasum", "-a", "256", "-c", "LARGE.sha256"])


def compare_records(script, archive, regex, *args):
    """Compare order-independent output from the parallel traversals."""
    fresh = subprocess.check_output([sys.executable, "-Bu", script, *args], cwd=ROOT, text=True)
    pattern = re.compile(regex)
    def records(text):
        selected = []
        for line in text.splitlines():
            match = pattern.search(line)
            if match:
                selected.append(line[match.start():])
        return sorted(selected)

    new = records(fresh)
    old = records(archive.read_text())
    if not old or new != old:
        raise RuntimeError(f"Closure records differ from {archive.relative_to(ROOT)}")
    for line in new:
        print(line)


def local_controls():
    regenerate_box_leaves()
    stage = ACCESS / "stage_graph"
    box = ACCESS / "box_structure"
    route = ACCESS / "r2_routes"
    period = ACCESS / "class_compat"
    period_check = ACCESS / "class_compat_check"
    independent = ACCESS / "restricted_theorem_check"
    print("Finite-route reader and original-lattice controls", flush=True)
    python(stage / "exact_cells.py", expected=stage / "cell-lotteries.json")
    python(box / "marked_weights.py", expected=box / "marked-weights.json")
    python(box / "snapshot.py", "--small", "marked_clock.py", expected=box / "marked-clock.json")
    python(route / "control_weights.py", expected=route / "control-weights.json")
    python(route / "box_route.py", expected=route / "box-route.txt")

    python(independent / "check_arithmetic.py", expected=independent / "arithmetic.log")
    with TemporaryDirectory(prefix="p8-reader-") as directory:
        reader = Path(directory) / "reader"
        build(independent / "reader.cpp", reader, "-ftrapv")
        python(independent / "check_controls.py", reader, expected=independent / "controls.log")
    python(period / "premises.py", expected=period / "premises.json")
    python(period / "period66_witness.py", expected=period / "period66-witness.out")
    python(period_check / "witness.py", expected=period_check / "witness.json")
    python(period_check / "local_controls.py", expected=period_check / "local-controls.json")
    python(AUDIT / "local_period.py")
    python(AUDIT / "terminal_witness.py")

    print("Exact cone comparisons and selection controls", flush=True)
    selection = ACCESS / "s1_selection"
    python(selection / "exact.py", expected=selection / "exact.json")
    python(selection / "larger_cell.py", expected=selection / "larger-cell.json")


def mass_controls(quick):
    basin = ACCESS / "e1_basins"
    event = ACCESS / "e1c_row181"
    entry = ACCESS / "en11b"
    destinations = ACCESS / "destinations"
    with compute_lock():
        run(["bash", ACCESS / "entry_pi/reproduce.sh"])
    for script, args in [("verify.py", [ACCESS / "entry_pi/catalogue.json"]),
                         ("survivor.py", []), ("singular_seam.py", [])]:
        python(destinations / script, *args)
    python(AUDIT / "basins.py")
    python(AUDIT / "mass.py")
    if quick:
        return

    print("Exact basin certificates, entry integration and three-step event", flush=True)
    with TemporaryDirectory(prefix="p8-mass-") as directory:
        reader = Path(directory) / "density"
        build(entry / "density.cpp", reader)
        with compute_lock():
            python(entry / "entrance.py", "--depth", "4", "--evaluator", reader,
                   expected=entry / "entrance.json")
    python(basin / "reproduce.py", "--recompute-mass")
    python(ACCESS / "e1b_guarded/reproduce.py", expected=ACCESS / "e1b_guarded/verification.log")
    archive = (event / "verification.log").read_bytes()
    first_cell = re.search(rb"^\(\d, \d, \d\) \d: ", archive, re.MULTILINE)
    if first_cell is None:
        raise RuntimeError("Three-step event archive is missing its cell records")
    expected = event / "expected-event.out"
    expected.write_bytes(archive[first_cell.start():])
    python(event / "reproduce.py", "--bound-only", expected=expected)
    python(AUDIT / "event_mass.py", "--all", expected=AUDIT / "event-mass.json")


def rank_controls():
    stage = ACCESS / "z1_stages"
    continuation = ACCESS / "z1b_rank_one"
    print("Six-root inherited ingress and primary rank-one graph", flush=True)
    compare_records(stage / "cover.py", stage / "rank-three.out", r"^closed ",
                    "--rank", "3", "--jobs", "6", "--limit", "40000")
    compare_records(continuation / "verify.py", continuation / "verification.out",
                    r"rank-one closed ")
    print("Six-root ordinary-difference reconstruction", flush=True)
    independent = ACCESS / "restricted_theorem_check"
    combined = independent / "expected-coverage.out"
    combined.write_text((independent / "coverage.log").read_text()
                        + (independent / "coverage-odd.log").read_text())
    compare_records(independent / "cover_exact.py", combined,
                    r"^line closed ", "--rank", "2")

    print("First-coordinate decrease for the four compared lifts", flush=True)
    anchor = ACCESS / "s1c_pointwise"
    fresh = subprocess.check_output([sys.executable, "-Bu", anchor / "anchor_bound.py",
                                     "--roots", "3", "5"], text=True)
    prefix = "anchor-clock-summary "
    row = next(line[len(prefix):] for line in fresh.splitlines() if line.startswith(prefix))
    assert json.loads(row) == json.loads((anchor / "anchor-bound.json").read_text())
    print(row)


def period_controls():
    primary = ACCESS / "class_compat"
    second = ACCESS / "class_compat_check"
    print("Centered period inventories and exceptional terminal sections", flush=True)
    python(primary / "stride_inventory.py", expected=primary / "stride-inventory.json")
    python(primary / "prime_reader.py", expected=primary / "prime-reader.json")
    python(second / "ingress.py", expected=second / "ingress.json")
    python(second / "ingress.py", "--roots", "0", "2", "4", "--canonical",
           "--index-11-at-rank-three", expected=second / "canonical.json")
    python(AUDIT / "period_audit.py", expected=AUDIT / "period-inventory.json")


def physical_tests():
    print("Prospective physical diagnostics (empirical coverage)", flush=True)
    route = ACCESS / "r2_routes"
    selection = ACCESS / "s1_selection"
    python(route / "run.py", expected=route / "diagnostic-results.json")
    with TemporaryDirectory(prefix="p8-holdout-") as directory:
        python(selection / "run.py", "holdout", directory)
        for test in ("holdout", "local"):
            python(selection / "summarize.py", test, directory,
                   expected=selection / f"{test}.json", statistical=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "quick", "local", "mass", "rank", "period", "physical"),
                        default="all")
    part = parser.parse_args().part
    run(["clang++", "--version"])
    run(["yices-smt2", "--version"])
    import sympy
    print("Python", sys.version.split()[0], "| SymPy", sympy.__version__, flush=True)
    if part in ("all", "quick", "local"):
        local_controls()
    if part in ("all", "quick", "mass"):
        mass_controls(part == "quick")
    if part in ("all", "rank"):
        rank_controls()
    if part in ("all", "period"):
        period_controls()
    if part in ("all", "physical"):
        physical_tests()
    print(f"P8 reproduction passed ({part}).")


if __name__ == "__main__":
    main()
