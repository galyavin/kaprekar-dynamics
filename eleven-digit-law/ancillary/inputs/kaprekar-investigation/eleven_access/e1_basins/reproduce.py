"""Replay the exact certificates from a private task snapshot."""
import argparse
import fcntl
from fractions import Fraction as Q
import json
from math import comb, factorial
from pathlib import Path
from shutil import copyfile
from subprocess import run
import sys
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def invoke(*args, **kwargs):
    return run(args, check=True, text=True, capture_output=True, **kwargs).stdout


def fixed_decimal(value, digits, up=False):
    unit = 10**digits
    integer, remainder = divmod(value.numerator*unit,value.denominator)
    integer += bool(up and remainder)
    return f"{integer//unit}.{integer%unit:0{digits}d}"


def ledger(output, archive):
    first, *lines = output.splitlines()
    _, grid, _, through, _, ambiguous, _, outside = first.split()
    grid, through, ambiguous, outside = map(int,(grid,through,ambiguous,outside))
    entries = [tuple(map(int,line.split())) for line in lines]
    assert [t for t,_,_ in entries] == list(range(1,through+1))
    assert ambiguous + outside + sum(count for _,count,_ in entries) == comb(grid,5)
    combined = sum(n for t,_,n in entries if t >= 5)
    extra = Q(factorial(11)*combined,(2*grid)**5*(24*grid)**6)
    baseline = Q(json.loads(archive.read_text())[4]["by_step_each_class"])
    bound = baseline + extra
    assert extra == Q(16984163369908523775525,604462909807314587353088)
    print("Original mass of certified boxes entering at steps 5-19 (exact):", extra)
    print("Increase in proved lower bound, percentage points (rounded down):",
          fixed_decimal(100*extra,10))
    print("m(D) in percent (proved lower bound):",fixed_decimal(100*bound,10))
    print("m(D^c) in percent (proved upper bound):",fixed_decimal(100*(1-bound),10,up=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute-mass", action="store_true",
                        help="compile and rerun the full native scan under the compute lock")
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--solver", default="yices-smt2")
    args = parser.parse_args()

    inputs = (
        HERE/"catalogue_checks.py",
        HERE/"ray_check.py",
        HERE/"folded_mass.cpp",
        HERE/"mass.out",
        HERE.parent/"entry_pi/catalogue.json",
        HERE.parent/"global_work/evidence/energy.json.gz",
        HERE.parent/"en11b/entrance.json",
        ROOT/"docs/drafts/research-direction/astra-strategy-13/basins.out",
    )
    with TemporaryDirectory(prefix="e1-basins-") as tmp:
        temp = Path(tmp)
        for file in inputs:
            copyfile(file,temp/file.name)
        catalogue = invoke(sys.executable,"-B","catalogue_checks.py",
                            "--catalogue","catalogue.json","--energy","energy.json.gz",
                            "--astra","basins.out","--diagnose",cwd=temp)
        print(catalogue,end="")
        print(invoke(sys.executable,"-B","ray_check.py","--solver",args.solver,cwd=temp),end="")
        output = (temp/"mass.out").read_text()
        if args.recompute_mass:
            with (ROOT/".task-monitor/compute.lock").open("a") as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                invoke(args.compiler,"-std=c++17","-O3","-Wall","-Wextra","-Werror",
                       "folded_mass.cpp","-o","folded_mass",cwd=temp)
                fresh = invoke("./folded_mass","128","19",cwd=temp)
            assert fresh == output
            print("Exact folded-box scan reproduced.")
        ledger(output,temp/"entrance.json")


if __name__ == "__main__":
    main()
