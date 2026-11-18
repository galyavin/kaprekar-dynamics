"""Unbounded residue and rare-predecessor checks for the joint-law proof."""
if not __debug__:
    raise RuntimeError(
        'Proof audits require assertions; disable -O/-OO and PYTHONOPTIMIZE.')

import json
from pathlib import Path
import z3

from five_clock import margin_step
from five_excursions import BRANCHES


def audit():
    B, u, v = z3.Ints('B u v')
    def ev(row):
        return row[0] * u + row[1] * v + row[2] * B + row[3]
    def rare(x, y):
        return z3.Or(z3.And(2*x == B, 2*y == B-2),
                     z3.And(2*x == B+2, y == 0),
                     z3.And(2*x == B, y == 0))
    def core(x, y):
        return z3.And(1 <= x, x < B, 0 <= y, y <= x,
                      2*y <= B-1, 2*(x-y) <= B-3,
                      x+y-B >= -2, 2*x-y-B >= -2)
    checks = {}
    for label, word in BRANCHES.items():
        un, vn = [ev(row) for row in word.matrix[:2]]
        guards = [ev(g) >= 0 for g in word.guards]
        s = z3.Solver()
        s.add(B >= 30, B % 2 == 0, *guards,
              rare(un, vn), z3.Not(rare(u, v)))
        assert s.check() == z3.unsat, ('outside rare predecessor', label)
        r = z3.Solver()
        r.add(B >= 30, B % 3 == 1, core(u, v), *guards,
              ((u+v) % 3 == 0) != ((un+vn) % 3 == 0))
        assert r.check() == z3.unsat, ('core label invariance', label)
        checks[label] = {'outside_rare_predecessor': 'unsat',
                         'core_residue_change': 'unsat'}
    # Verify the finite terminal box, including both directions of the label
    # criterion. Completeness and entry into this box are already proved.
    terminal = 0
    for a in range(-2, 10):
        for b in range(-2, 10):
            q, path = (a, b), {}
            while q not in path:
                path[q] = len(path)
                q = margin_step(q)
            period = len(path)-path[q]
            residue = (a+b) % 3
            expected = 1 if residue == 0 else 4 if residue == 2 else (
                4 if a % 3 == b % 3 == 2 else 3)
            assert period == expected
            terminal += 1
    # Four strict leading cells. The corresponding integer maps differ only
    # by constant offsets; determinant ±1 gives a bijection modulo 3.
    matrices = [((0,-1),(-1,-1)), ((0,-1),(1,-1)),
                ((1,0),(1,-1)), ((1,0),(1,1))]
    dets = [a*d-b*c for ((a,b),(c,d)) in matrices]
    assert set(dets) <= {-1, 1}
    for m in matrices:
        images = {(sum(m[0][i]*q[i] for i in range(2)) % 3,
                   sum(m[1][i]*q[i] for i in range(2)) % 3)
                  for q in ((a,b) for a in range(3) for b in range(3))}
        assert len(images) == 9
    result = {'solver': z3.get_version_string(),
              'domain': 'all integer B >= 30; parity/residue restrictions as stated',
              'branch_checks': checks, 'terminal_margin_starts': terminal,
              'strict_leading_cell_determinants': dets,
              'rare_basin_input_weight': '40*B^2-80*B-10',
              'rare_actual_histogram': {'0': '2', '1': '35*B^2-75*B-12',
                                        '2': '5*B^2-5*B'},
              'scope': 'Encoded residue/predecessor implications and terminal checks; '
                       'weighted equidistribution and limit proof are analytic, not SMT or Lean.'}
    Path(__file__).with_name('five_joint_audit_results.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(f'All-base checks passed on {len(checks)} sort branches; '
          f'{terminal} terminal starts; four residue bijections.')


if __name__ == '__main__':
    audit()
