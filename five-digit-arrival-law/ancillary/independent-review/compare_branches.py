"""Comparison of the project's 26-branch transition encoding with the
reviewer's sorting-network encoding (run only AFTER the reviewer's own checks).

1. The project's branch guards partition the valid domain (exhaustive and
   pairwise disjoint) for every B >= 2 -- a property the project's audits rely
   on but do not check in the solver: a gap would silently drop orbits from
   the disjunctive transition relation and make the six-step and (36) checks
   vacuously UNSAT.
2. On each branch, the branch's affine image equals the reviewer's
   sorting-network image (so the two transition relations coincide).

The project modules are imported read-only from a scratch copy so that no
bytecode is written into the project tree.
usage: python compare_branches.py PATH_TO_SCRATCH_COPY_OF_SOURCE
"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import z3                                   # noqa: E402
from five_excursions import BRANCHES        # noqa: E402  (project)
import z3_independent as R                  # noqa: E402  (reviewer)


def ev(row, u, v, B):
    return row[0] * u + row[1] * v + row[2] * B + row[3]


B, u, v = z3.Ints("B u v")
valid = z3.And(B >= 2, 1 <= u, u <= B - 1, 0 <= v, v <= u)
G = {lab: z3.And([ev(g, u, v, B) >= 0 for g in w.guards]) for lab, w in BRANCHES.items()}
print("project branches:", len(BRANCHES))

s = z3.Solver()
s.add(valid, z3.Not(z3.Or(list(G.values()))))
r = s.check()
print(f"[{'PASS' if r == z3.unsat else 'FAIL'}] branches cover every valid pair, all B>=2: {r}")

labs = list(G)
overl = 0
for i in range(len(labs)):
    for j in range(i + 1, len(labs)):
        s = z3.Solver()
        s.add(G[labs[i]], G[labs[j]])
        if s.check() != z3.unsat:
            overl += 1
            print("   overlap", labs[i], labs[j], s.model())
print(f"[{'PASS' if overl == 0 else 'FAIL'}] branches pairwise disjoint ({overl} overlaps)")

bad = 0
for lab, w in BRANCHES.items():
    s = z3.Solver()
    s.add(B >= 2, G[lab])
    u2, v2 = R.step(s, B, u, v, "cmp")
    s.add(z3.Or(u2 != ev(w.matrix[0], u, v, B), v2 != ev(w.matrix[1], u, v, B)))
    if s.check() != z3.unsat:
        bad += 1
        print("   mismatch on branch", lab, s.model())
print(f"[{'PASS' if bad == 0 else 'FAIL'}] each branch image equals the sorting-network "
      f"image ({bad} mismatches)")
