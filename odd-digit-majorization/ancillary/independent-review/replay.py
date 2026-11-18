"""Rebuild and check the independent computations in a disposable working copy."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ARCHIVE = Path(__file__).resolve().parent


def execute(work, name, args):
    start = time.monotonic()
    print(f"Starting {name}", flush=True)
    log = work / f"rerun_{name}.txt"
    with log.open("w") as out:
        result = subprocess.run(args, cwd=work, env={**os.environ, "PYTHON": sys.executable},
                                stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        print(log.read_text()[-4000:], file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, args)
    print(f"{name} {time.monotonic()-start:.1f}s", flush=True)
    return log.read_text()


def require(ok, name):
    if not ok:
        raise RuntimeError(f"Unexpected result for {name}")


def canonical(log):
    return [line for line in log.splitlines() if not line.startswith("[lock acquired")]


def same_summary(actual, expected):
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            same_summary(actual[key], value) for key, value in expected.items())
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            same_summary(a, b) for a, b in zip(actual, expected))
    if isinstance(expected, float):
        return math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10)
    return actual == expected


def replay(work):
    execute(work, "compile", [os.environ.get("CC", "cc"), "-O2", "kap9.c", "-o", "kap9"])
    scripts = sorted(p.stem for p in work.glob("t*.py"))
    scripts += ["z_entry", "z_clock", "z_reseed", "z_vacuity", "z_vacuity_iv"]

    def check(name):
        output = execute(work, name, [sys.executable, "-B", name + ".py"])
        if name.startswith("t"):
            expected = (ARCHIVE / f"logs_{name}.txt").read_text()
            require(canonical(output) == canonical(expected), name)
        elif name == "z_vacuity":
            unexpected = [line for line in output.splitlines() if "**UNEXPECTED**" in line]
            require(len(unexpected) == 1 and "min D >= 2 eta" in unexpected[0]
                    and "unsat" in unexpected[0], name)
            require(output.count("hypotheses rank") == 8 and output.count("j=4: unsat") == 2, name)
        elif name == "z_vacuity_iv":
            require(len(re.findall(r"^\(iv\).*\bsat\b", output, re.M)) == 3, name)
        else:
            expect = {"z_entry": (10, 1), "z_clock": (28, 0), "z_reseed": (76, 2)}
            verdicts = re.findall(r"\b(unsat|sat)\s+OK\s+", output)
            require((verdicts.count("unsat"), verdicts.count("sat")) == expect[name]
                    and "**FAIL**" not in output, name)

    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(check, scripts))

    for name in ("descent", "sims", "pred", "allroutes", "targeted"):
        output = execute(work, name, ["bash", "run_" + name + ".sh"])
        if name in ("descent", "sims"):
            expected = (ARCHIVE / f"logs_{name}.txt").read_text()
            require(canonical(output) == canonical(expected), name)

    for name in ("compare_uniform", "summarize_route"):
        output = execute(work, name, [sys.executable, "-B", name + ".py"])
        expected = (ARCHIVE / f"logs_{name}.txt").read_text()
        require(canonical(output) == canonical(expected), name)

    for original in ARCHIVE.glob("*.json"):
        actual = json.loads((work / original.name).read_text())
        expected = json.loads(original.read_text())
        require(same_summary(actual, expected), original.name)
    print("Independent checks reproduced the archived results.", flush=True)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="kaprekar-review-") as tmp:
        root = Path(tmp)
        shutil.copytree(ARCHIVE.parent / "source", root / "source")
        work = root / "independent-review"
        work.mkdir()
        for pattern in ("*.py", "*.c", "*.sh"):
            for src in ARCHIVE.glob(pattern):
                shutil.copy2(src, work / src.name)
        replay(work)
