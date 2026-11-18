"""Exact, unbounded integer checks for the five-digit clock proof.

Requires z3-solver. UNSAT checks range over all integer bases B >= 30;
these are not finite-base sampling. The human proof is five_digit_limit.md.
The strip-entry lemma is a computer-assisted finite affine case check.
"""
if not __debug__:
    raise RuntimeError(
        'Proof audits require assertions; disable -O/-OO and PYTHONOPTIMIZE.')

import json
from fractions import Fraction
from pathlib import Path

import z3

from five_clock import margin_step
from five_excursions import BRANCHES


def audit():
    B = z3.Int('B')
    def ev(row, u, v):
        return row[0]*u+row[1]*v+row[2]*B+row[3]
    def valid(u, v):
        return z3.And(1 <= u, u < B, 0 <= v, v <= u)
    def H(u, v):
        return z3.And(valid(u, v), 2*v <= B-1, 2*(u-v) <= B-3,
                      u+v-B >= -2, 2*u-v-B >= -2)
    def T2(u, v):
        return z3.Or(z3.And(2*u == B, 2*v == B-2),
                     z3.And(2*u == B+2, v == 0))
    def unsat(solver, description):
        result = solver.check()
        if result != z3.unsat:
            witness = str(solver.model()) if result == z3.sat else solver.reason_unknown()
            raise AssertionError((description, str(result), witness))
    u, v = z3.Ints('u v')
    domain = z3.And(B >= 30, valid(u, v))
    guards = {label: z3.And([ev(row, u, v) >= 0 for row in word.guards])
              for label, word in BRANCHES.items()}
    assert len(guards) == 26

    def coverage_query(included):
        solver = z3.Solver()
        solver.add(domain, z3.Not(z3.Or([guards[label] for label in included])))
        return solver

    unsat(coverage_query(guards), 'uncovered valid state')
    overlap = z3.Solver()
    overlap.add(domain, z3.Sum([z3.If(guard, 1, 0) for guard in guards.values()]) >= 2)
    unsat(overlap, 'valid state covered by multiple branches')
    print('26 branch guards, all B >= 30: uncovered state UNSAT; overlapping guards UNSAT.')

    omitted = 'Z0'
    gap = coverage_query(label for label in guards if label != omitted)
    status = gap.check()
    if status != z3.sat:
        raise AssertionError(('missing branch negative control', omitted, str(status)))
    model = gap.model()
    witness = {str(x): model.eval(x).as_long() for x in (B, u, v)}
    assert BRANCHES[omitted].contains(witness['B'], (witness['u'], witness['v']))
    print('Without', omitted, 'an uncovered state is SAT:', witness)

    active = []
    for label, w in BRANCHES.items():
        s = z3.Solver()
        s.add(B >= 30, H(u, v), *[ev(g, u, v) >= 0 for g in w.guards])
        result = s.check()
        if result == z3.unsat:
            continue
        assert result == z3.sat
        active.append(label)
        un, vn = [ev(row, u, v) for row in w.matrix[:2]]
        a, b = u+v-B, 2*u-v-B
        an = z3.If(a >= 1, b, z3.If(a >= -1, b, b+1))
        bn = z3.If(a >= 1, a-3, z3.If(a >= -1, -2*a, 3))
        s.push()
        s.add(z3.Or(un+vn-B != an, 2*un-vn-B != bn))
        unsat(s, ('core map', label))
        s.pop()
        s.add(z3.Not(H(un, vn)))
        unsat(s, ('core invariance', label))
    print('Core map and forward invariance: all B >= 30; branches', active)

    for label, w in BRANCHES.items():
        s = z3.Solver()
        s.add(B >= 30, v >= 1, *[ev(g, u, v) >= 0 for g in w.guards])
        un, vn = [ev(row, u, v) for row in w.matrix[:2]]
        s.add(z3.Or(z3.Not(valid(un, vn)), 2*un < B-1, 2*(un-vn) > B+2))
        unsat(s, ('first-step bounds', label))
    print('First nonzero-coordinate step: universal strip bounds verified.')

    # The entire finite strip-entry statement, including zero/equality cases.
    # Fresh coordinates avoid substitutions that could grow exponentially.
    s = z3.Solver()
    s.set(timeout=45000)
    u, v = z3.Ints('u0 v0')
    s.add(B >= 30, valid(u, v), 2*u >= B-1, 2*v <= B+1,
          2*(u-v) <= B+2)
    for k in range(7):
        s.add(z3.Not(H(u, v)), z3.Not(T2(u, v)))
        if k == 6:
            unsat(s, 'strip avoids H and T2 for six steps')
            break
        un, vn = z3.Ints(f'u{k+1} v{k+1}')
        s.add(z3.Or(*[
            z3.And(*[ev(g, u, v) >= 0 for g in w.guards],
                   un == ev(w.matrix[0], u, v), vn == ev(w.matrix[1], u, v))
            for w in BRANCHES.values()]))
        u, v = un, vn
    print('Strip entry into H or T2 within six steps: UNSAT for all B >= 30.')

    # Independent finite terminal graph of the base-free margin map.
    cycles, errors, tails = set(), [], []
    for a in range(-2, 10):
        for b in range(-2, 10):
            q, path, index = (a, b), [], {}
            while q not in index:
                index[q] = len(path)
                path.append(q)
                q = margin_step(q)
                assert len(path) < 100
            h = index[q]
            cycle = path[h:]
            j = cycle.index(min(cycle))
            cycles.add(tuple(cycle[j:]+cycle[:j]))
            tails.append(h)
            errors.append(Fraction(h)-Fraction(2, 3)*max(a, b))
            assert max(map(max, path)) <= 10
            assert max(margin_step(margin_step((a, b)))) <= 9
    assert max(tails) == 10 and min(errors) == -3 and max(errors) == Fraction(16, 3)
    assert len(cycles) == 4
    print('144 terminal starts: max tail 10; clock discrepancy [-3, 16/3].')

    # Exact two-step decrement outside the terminal box.
    a, b = z3.Ints('a b')
    def f(a, b):
        return (z3.If(a >= 1, b, z3.If(a >= -1, b, b+1)),
                z3.If(a >= 1, a-3, z3.If(a >= -1, -2*a, 3)))
    maximum = lambda a, b: z3.If(a >= b, a, b)
    minimum = lambda a, b: z3.If(a <= b, a, b)
    a2, b2 = f(*f(a, b))
    s = z3.Solver()
    s.add(a >= -2, b >= -2, maximum(a, b) >= 10)
    s.add(maximum(a2, b2) != maximum(a, b)-3+z3.If(minimum(a, b) == -2, 1, 0))
    unsat(s, 'bulk two-step decrement')
    print('Bulk maximum drops by 3, or by 2 at the sole -2 correction.')
    result = {
        'solver': z3.get_version_string(), 'base_domain': 'all integer B >= 30',
        'branch_count': len(guards), 'uncovered_state': 'unsat', 'overlapping_guards': 'unsat',
        'coverage_negative_control': {'omitted': omitted, 'uncovered_state': 'sat',
                                      'witness': witness},
        'core_branches': active, 'core_map_and_invariance': True,
        'first_step_strip_bounds': True, 'strip_entry_steps': 6,
        'strip_entry_unbounded_smt': 'unsat', 'terminal_starts': 144,
        'terminal_cycles_in_margin_coordinates': sorted(cycles),
        'terminal_max_tail': 10, 'terminal_error_range': [str(min(errors)), str(max(errors))],
        'bulk_decrement_verified': True,
        'pair_clock_error_bound_on_core': ['-3', '6'],
        'actual_time_upper_bound': 'B/3 + ceil(sqrt(B)) + 20',
        'scope': 'SMT certifies the stated affine implications; analytic and counting arguments are in the proof note. Not a Lean formalization.',
    }
    Path(__file__).with_name('five_clock_audit_results.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    audit()
