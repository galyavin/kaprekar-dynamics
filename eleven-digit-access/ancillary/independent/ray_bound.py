"""Check the sharp ten-step maximum-norm expansion of the ray return."""
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys
import tempfile


def apply(x):
    v, c, d, z = x
    return (-d-z, max(v, -c-d), min(v, -c-d), -v-c)


def build_query():
    smt = ["(set-logic QF_LRA)"]
    for t in range(11):
        for key in "vcdz":
            smt.append(f"(declare-fun {key}{t} () Real)")
    smt.append("(assert (>= c0 d0))")
    for name in "vcdz":
        smt.append(f"(assert (and (>= {name}0 (- 1)) (<= {name}0 1)))")
    faces = [f"(= {name}0 {end})" for name in "vcdz" for end in ("1", "(- 1)")]
    smt.append("(assert (or " + " ".join(faces) + "))")
    for t in range(10):
        v, c, d, z = [f"{k}{t}" for k in "vcdz"]
        alt = f"(- (+ {c} {d}))"
        test = f"(>= {v} {alt})"
        terms = [f"(- (+ {d} {z}))", f"(ite {test} {v} {alt})",
                 f"(ite {test} {alt} {v})", f"(- (+ {v} {c}))"]
        for k, term in zip("vcdz", terms):
            smt.append(f"(assert (= {k}{t+1} {term}))")
    for name in "vcdz":
        val = f"(* 77 {name}10)"
        smt.append(f"(assert (and (> {val} (- 295)) (< {val} 295)))")
    return "\n".join(smt + ["(check-sat)", ""])


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    x = tuple(Q(n, 77) for n in (-53, 12, -15, 77))
    assert x[1] >= x[2] and max(map(abs, x)) == 1
    for _ in range(10):
        x = apply(x)
    assert max(map(abs, x)) == Q(295, 77)
    assert x == tuple(Q(n, 77) for n in (295, 256, 57, -295))
    print("Attaining image:", *x, flush=True)
    with tempfile.TemporaryDirectory() as directory:
        query = Path(directory) / "strict.smt2"
        query.write_text(build_query())
        for solver in sys.argv[1:] or ["yices-smt2", "z3"]:
            result = subprocess.run([solver, str(query)], check=True, text=True,
                                    capture_output=True)
            assert result.stdout.strip() == "unsat", result
            print(solver, result.stdout.strip(), flush=True)
