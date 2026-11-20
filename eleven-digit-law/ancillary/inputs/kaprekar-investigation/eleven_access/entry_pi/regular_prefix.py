"""Exact finite capture or repetition on the regular rotation seams."""
import argparse
import json
from pathlib import Path
from subprocess import run

import z3

from regular_geometry import active_at_boundary, dot


def rule(word, turns):
    active = active_at_boundary(word)
    solver = z3.SolverFor("QF_LRA")
    v = tuple(z3.Real(f"v_{j}") for j in range(6))
    solver.add(sum(v) == 0)
    comparisons = [tuple(dot(form, v) for _, _, form, _ in active)]
    for step in range(1, turns * len(active) + 1):
        k = (step - 1) % len(active)
        value = comparisons[-1][k]
        correction = z3.Real(f"min_{step}")
        solver.add(correction == z3.If(value < 0, -value, 0))
        next_values = tuple(z3.Real(f"q_{step}_{j}") for j in range(len(active)))
        for new, old, (_, _, form, _) in zip(next_values, comparisons[-1], active):
            solver.add(new == old + dot(form, active[k][3]) * correction)
        comparisons.append(next_values)
    sections = comparisons[::len(active)]
    # A convergent tail outside the seams cannot already be captured or
    # have locally periodic normal coordinates.
    for section in sections:
        solver.add(z3.Or(*(q < 0 for q in section)))
    for i, a in enumerate(sections):
        for b in sections[i + 1:]:
            solver.add(z3.Or(*(x != y for x, y in zip(a, b))))
    return solver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("min_dim", type=int)
    parser.add_argument("max_dim", type=int)
    parser.add_argument("turns", type=int)
    args = parser.parse_args()
    rows = json.loads(Path("catalogue.json").read_text())
    for row in sorted(rows, key=lambda r: (-r["dimension"], r["word"])):
        word = row["word"]
        if (len(word) != 5 or any("-" in phase for phase in word)
                or not args.min_dim <= row["dimension"] <= args.max_dim):
            continue
        if not row["paired"] or row["dimension"] == 5:
            continue
        solver = rule(word, args.turns)
        solver.set(timeout=20000)
        z3_status = solver.check()
        print(row["dimension"], " ".join(word), z3_status, flush=True)
        text = solver.to_smt2()
        Path("prefix-" + "_".join(word) + ".smt2").write_text(text)
        if z3_status == z3.sat:
            model = solver.model()
            print("  v:", [str(model.eval(z3.Real(f"v_{j}"))) for j in range(6)], flush=True)
        if z3_status != z3.unsat:
            break
        yices = run(["yices-smt2"], input="(set-logic QF_LRA)\n" + text,
                    capture_output=True, text=True, check=True, timeout=20)
        if yices.stdout.strip() != "unsat":
            raise RuntimeError(yices.stdout + yices.stderr)
        print("  Yices agrees", flush=True)


if __name__ == "__main__":
    main()
