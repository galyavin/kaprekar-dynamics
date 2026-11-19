"""Verify the supplied evidence and run the proof checks in a private copy."""

import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from time import monotonic


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def command(*args, cwd=None, capture=False):
    return subprocess.run(args, cwd=cwd, check=True, text=True,
                          stdout=subprocess.PIPE if capture else None).stdout


def python(directory, script, *args, capture=False):
    return command(sys.executable, "-B", script, *map(str, args),
                   cwd=directory, capture=capture)


def compressed(path):
    return gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()


def verify_manifest(here):
    expected = json.loads((here / "results" / "artifact_sha256.json").read_text())
    actual = {p.relative_to(here).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(here.rglob("*")) if p.is_file()
              and p.name != "artifact_sha256.json"}
    require(actual == expected, "The supplied file inventory or SHA-256 values differ.")
    print(f"Verified {len(actual)} input files.", flush=True)


def check(state, result, z3_cli, yices, workers):
    import sympy
    import z3

    from importlib.metadata import version
    require(version("z3-solver") == "5.4.1.0" and sympy.__version__ == "1.14.0"
            and version("mpmath") == "1.3.0", "Use the pinned Python dependencies.")
    z3_version = command(z3_cli, "--version", capture=True).splitlines()[0]
    yices_version = command(yices, "--version", capture=True).splitlines()[0]
    require(re.search(r"\b5\.4\.1\b", z3_version) and yices_version == "Yices 2.7.0",
            "Z3 5.4.1 and Yices 2.7.0 are required.")
    print(sys.version.split()[0], z3_version, yices_version, flush=True)
    original = state / "original"
    local = state / "local"
    generator = state / "generator"
    independent = state / "independent"
    centered = state / "centered"

    # Discover and check the local obligations before the long exhaustions.
    radius = result / "radius-queries"
    python(local, "radius.py", original / "energy.json.gz", radius)
    radius_queries = sorted(radius.glob("*.smt2"))
    require(len(radius_queries) == 6, "Expected six original radius predicates.")
    for query in radius_queries:
        for solver in (z3_cli, yices):
            actual = command(solver, str(query), capture=True)
            require(actual.strip() == "unsat", f"Radius obligation: {query.name}: {actual}")
    print("Six radius predicates: both solvers UNSAT.", flush=True)
    radius_independent = local / "independent_radius"
    output = python(radius_independent, "yices_words.py", capture=True)
    require(output.count("-> unsat") == 5, output)
    print(output, end="")
    output = python(radius_independent, "yices_ray.py", 10, 2, capture=True)
    require(output.strip() == "unsat", output)
    print("Independent ten-step ray predicate: UNSAT.", flush=True)
    python(independent, "ray_bound.py", z3_cli, yices)
    python(local, "seam_certificate.py", "--output", result / "seam-queries")
    python(local, "seam_independent.py")
    python(local, "cycles.py")
    python(local, "basin_certificate.py", "--output", result / "four-cycle-basin.json")
    python(independent, "basin_reverse.py")
    python(centered, "eleven_two_cone_audit.py")
    python(centered, "eleven_two_cone_checks.py")

    frozen = independent / "frozen"
    regenerated = result / "independent"
    python(independent, "exhaust.py", regenerated, "--frozen-certificates",
           frozen / "certificates.jsonl.gz")
    python(independent, "audit.py", regenerated)
    for supplied in sorted(frozen.iterdir()):
        require(compressed(supplied) == compressed(regenerated / supplied.name),
                f"Exact independent check differs at {supplied.name}.")
    print("All five independent records reproduced.", flush=True)
    python(independent, "compare.py", regenerated, original, result / "comparison.json")
    python(independent, "replay.py", original, result / "replay.json", z3_cli, yices)

    generated = result / "generator"
    generated.mkdir()
    generator_lock = generated / "compute.lock"
    generator_lock.touch()
    for scope in ("positive", "zero"):
        python(generator, "tree.py", 5, scope, generated / scope, generator_lock)
    python(generator, "energy.py", generated / "positive.json", generated / "zero.json",
           generated / "energy.json")
    python(generator, "lp.py", 5, generated / "positive.json", generated / "work.json",
           generated / "energy.json")
    for stem in ("positive", "zero", "energy", "work"):
        plain = generated / f"{stem}.json"
        (generated / f"{stem}.json.gz").write_bytes(gzip.compress(plain.read_bytes(), mtime=0))
    python(independent, "compare.py", regenerated, generated,
           result / "generator-comparison.json")

    python(centered, "eleven_centered_clock_checks.py")
    python(centered, "eleven_centered_clock_audit.py", "--workers", workers)
    coverage = json.loads((centered / "eleven_centered_clock_audit_results.json").read_text())
    require(sum(len(case["records"]) for case in coverage["cases"]) == 97465,
            "Wrong number of centered coverage nodes.")
    require(coverage["status"] == "complete: all physical source cases and nodes accepted",
            "Centered coverage is incomplete.")
    print("97,465 centered coverage nodes accepted.", flush=True)


def main():
    require(__debug__, "Run with Python assertions enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4, help="parallel centered-coverage batches")
    args = parser.parse_args()
    require(args.workers >= 1, "Use at least one worker.")
    here = Path(__file__).resolve().parent
    verify_manifest(here)
    z3_cli = (os.environ.get("Z3") or shutil.which("z3")
              or str(Path(sys.executable).parent / "z3"))
    yices = shutil.which("yices-smt2")
    require(yices is not None, "Add yices-smt2 2.7.0 to PATH.")
    work = Path(tempfile.mkdtemp(prefix="kaprekar-eleven-"))
    state, result = work / "input", work / "results"
    shutil.copytree(here, state)
    result.mkdir()
    print(f"Private reproduction directory: {work}", flush=True)
    start = monotonic()
    lock_path = Path(os.environ.get("COMPUTE_LOCK", result / "compute.lock"))
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        check(state, result, z3_cli, yices, args.workers)
    print(f"All reproduction checks passed in {monotonic() - start:.1f} seconds.", flush=True)


if __name__ == "__main__":
    main()
