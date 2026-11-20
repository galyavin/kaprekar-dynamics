"""Replay the primary and independently reconstructed exact certificates."""
import fcntl
import os
from pathlib import Path
import subprocess
import sys


def run(work, lock_path):
    evidence = work / "kaprekar-investigation/cross_length"
    compiler = os.environ.get("CXX", "c++")
    env = dict(os.environ)
    env.pop("PYTHONOPTIMIZE", None)

    def command(args, **kw):
        return subprocess.run(args, check=True, env=env, **kw)

    def python(path, *args, **kw):
        return command([sys.executable, "-E", "-B", "-u", path, *args], cwd=work, **kw)

    def check_output(path, expected, *args):
        actual = python(path, *args, stdout=subprocess.PIPE).stdout
        if actual != expected.read_bytes():
            sys.stdout.buffer.write(actual)
            raise SystemExit(f"Changed output for {path.name}; expected {expected.name}.")
        print(f"{expected.name} reproduced", flush=True)

    with lock_path.open("a+") as lock:
        print("Waiting for exclusive compute lock.", flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
        command([sys.executable, "--version"])
        command([compiler, "--version"])

        print("Primary exact certificates", flush=True)
        python(work / "check_primary.py")
        python(evidence / "lc_classification/verify_root.py", "--literal")
        python(evidence / "lc2_order12/verify.py", "--literal")
        itinerary = evidence / "lc2_order12/itinerary.txt"
        check_output(evidence / "lc2_order12/reconstruct.py",
                     evidence / "lc2_order12/certificate.jsonl", itinerary)
        python(evidence / "lc3_all_orders/verify.py")
        python(evidence / "lc_classification/verify_25.py")
        python(evidence / "tl3_survival/verify.py")
        python(evidence / "e13b_passage/marked_cell.py", evidence / "e13_selection")

        print("Order-eight weights and finite cube experiments", flush=True)
        z1 = evidence / "z1_components"
        binary = work / "roots"
        command([compiler, "-std=c++17", "-O3", "-Wall", "-Wextra", "-Werror",
                 z1 / "roots.cpp", "-o", binary])
        check_output(z1 / "weights.py", z1 / "weights.json", binary)
        check_output(z1 / "physical.py", z1 / "physical.txt", evidence / "lc_classification")
        check_output(z1 / "mixture_panel.py", z1 / "mixture-panel.json",
                     binary, z1 / "weights.json")
        check_output(z1 / "parity.py", z1 / "parity.tsv", binary)

        print("Independent reconstruction", flush=True)
        independent = evidence / "survival_reconstruction"
        chunks = []
        for stage in ("strict", "wall", "class0", "cube"):
            result = python(independent / f"{stage}.py", work, stdout=subprocess.PIPE).stdout
            sys.stdout.buffer.write(result)
            sys.stdout.flush()
            chunks.append(result)
        if b"".join(chunks) != (independent / "reconstruction.log").read_bytes():
            raise SystemExit("Independent reconstruction differs from the archived log.")
        print("All primary checks passed. Independent reconstruction reproduced byte for byte.")


if __name__ == "__main__":
    run(Path(sys.argv[1]), Path(sys.argv[2]))
