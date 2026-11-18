"""All-base strip-clock estimate used in the sharp second-order maximum."""
if not __debug__:
    raise RuntimeError(
        'Proof audits require assertions; disable -O/-OO and PYTHONOPTIMIZE.')

import json
from pathlib import Path
import z3
from five_excursions import BRANCHES


def audit():
    B, u0, v0 = z3.Ints('B u0 v0')
    u, v = u0, v0
    solver = z3.Solver()
    solver.set(timeout=45000)
    solver.add(B >= 30, 1 <= u0, u0 < B, 0 <= v0, v0 <= u0,
               2*u0 >= B-1, 2*v0 <= B+1, 2*(u0-v0) <= B+2,
               4*u0 >= 3*B)
    def H(u, v):
        return z3.And(1 <= u, u < B, 0 <= v, v <= u,
                      2*v <= B-1, 2*(u-v) <= B-3,
                      u+v-B >= -2, 2*u-v-B >= -2)
    def ev(row, u, v):
        return row[0]*u+row[1]*v+row[2]*B+row[3]
    results = []
    for k in range(7):
        solver.push()
        solver.add(H(u, v), z3.Or(2*(u+v-B) > 2*u0-B+40,
                                  2*(2*u-v-B) > 2*u0-B+40))
        result = solver.check()
        assert result == z3.unsat, ('strip clock', k, result)
        results.append(str(result))
        solver.pop()
        solver.add(z3.Not(H(u, v)))
        if k < 6:
            un, vn = z3.Ints(f'u{k+1} v{k+1}')
            solver.add(z3.Or(*[
                z3.And(*[ev(g, u, v) >= 0 for g in word.guards],
                       un == ev(word.matrix[0], u, v),
                       vn == ev(word.matrix[1], u, v))
                for word in BRANCHES.values()]))
            u, v = un, vn
    data = {'solver': z3.get_version_string(),
            'domain': 'all integer B >= 30; q0 in J_B; 4*u0 >= 3*B',
            'first_core_entry_times_checked': list(range(7)),
            'results': results,
            'bound': 'max(a_entry,b_entry) <= u0-B/2+20',
            'scope': 'Additional strip-clock implication only. Entry within six '
                     'steps and rare-basin exclusion are proved in the existing audits.'}
    Path(__file__).with_name('five_quantitative_audit_results.json').write_text(
        json.dumps(data, indent=2)+'\n')
    print('Seven first-entry implications UNSAT, all integer B >= 30.')


if __name__ == '__main__':
    audit()
