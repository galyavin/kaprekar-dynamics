"""Independent real-linear seam checks. Run in the research Python with Yices on PATH."""
from fractions import Fraction as F
from itertools import combinations
import subprocess

import z3


def phases(a):
    d = (1 - 3*a)/2
    return [(1-a, 1-2*a, 2*d, d, a), (1-a, 1-2*a, 2*a+d, a+d, a)]


def raw(q, base, shift):
    return [base-shift, *q[:-1], q[-1]-shift,
            *(base-x-(shift if i else 0) for i, x in enumerate(q))]


def diff(digits):
    return [digits[i]-digits[-1-i] for i in range(5)]


def pairwise_sort(keys, entries):
    # Rank from pairwise comparisons with index tie-breaking. No comparator network.
    out = []
    for level in sorted(set(keys), reverse=True):
        group = [value for key, value in zip(keys, entries) if key == level]
        if len(group) == 1:
            out.extend(group)
            continue
        rank = [z3.Sum([z3.If(y >= x if j < i else y > x, 1, 0)
                        for j, y in enumerate(group) if j != i])
                for i, x in enumerate(group)]
        out.extend(z3.Sum([z3.If(r == j, x, 0) for r, x in zip(rank, group)])
                   for j in range(len(group)))
    return out


def actual(solver, q, base, leading):
    k, v = raw(leading, 1, 0), raw(q, base, 1)
    solver.add(base > q[0], q[-1] > 0, *[a >= b for a, b in zip(q, q[1:])])
    solver.add(*[v[i] > v[j] for i in range(11) for j in range(11) if k[i] > k[j]])
    return diff(pairwise_sort(k, v))


def max2(a, b):
    return z3.If(a >= b, a, b)


def paired(q, base):
    d = raw(q, base, 1)
    sorted_digits = [d[0]]
    for i in range(1, 6):
        a, b = d[i], d[-i]
        sorted_digits.extend([max2(a, b), a+b-max2(a, b)])
    return diff(sorted_digits)


def coords(q, base):
    x = [6*y-(5-i)*base for i, y in enumerate(q)]
    l = -x[0]
    return l, (x[1]+2*l)/6, (x[2]+3*l)/6, (x[3]+4*l)/6, (x[4]-l)/6


def query(name, phase, ret):
    p = phases(F(1, 6) if name == 'lower' else F(2, 11))
    B, *q = z3.Reals('B q1 q2 q3 q4 q5')
    slv = z3.Solver()
    first = actual(slv, q, B, p[phase])
    first_expected = paired(q, B)
    if phase == 1:
        first_expected[3] += max2(0, q[2]+q[3]-B+1)
    if ret:
        l, t_p, s, w, r = coords(q, B)
        j, k = max2(-r-1, 0)-max2(r-1, 0), max2(0, r+1)
        a, c = max2(-1, -r), max2(0, -l+2*s+1)
        t = s-t_p-c
        d = max2(a-1, t)
        e = max2(a, t-1)-max2(-t-1, -a)
        f = -max2(-a, 1-t)+j+1
        slv.add(w-l+t_p <= -1, l-2*w-2*k >= -1, w-j-t_p-1 <= 0)
        mid = z3.Reals('m1 m2 m3 m4 m5')
        slv.add(*[x == y for x, y in zip(mid, first)])
        second = actual(slv, mid, B, p[1])
        correction = max2(0, mid[2]+mid[3]-B+1)
        got = [*coords(second, B), second[3]-second[4]-q[3]+q[4],
               second[2]-second[3]-q[2]+q[3]]
        want = [l-6*d, e-2*d, f-3*d, w+j-t_p-4*d+correction,
                -t_p-k+d, r+j+k+correction, -c-max2(t-a-1, 0)-correction]
    else:
        got, want = first, first_expected
    # Check feasibility on each branch before asserting a discrepancy.
    assert slv.check() == z3.sat
    slv.add(z3.Or(*[x != y for x, y in zip(got, want)]))
    # No integer restriction: the queries include all physical offsets and ties.
    assert slv.check() == z3.unsat
    model = '(set-logic QF_LIRA)\n'+slv.to_smt2()
    # The ranking formula counts up to a fixed finite group with integer indicators.
    result = subprocess.run(['yices-smt2', '--timeout=15'],
                            input=model, capture_output=True, text=True, check=True)
    assert result.stdout.strip() == 'unsat', result
    return 'guard SAT in Z3, discrepancy UNSAT in Z3 and Yices'


def main():
    # Every open-edge strict comparison follows from affine endpoint signs.
    for phase in range(2):
        keys = [raw(phases(x)[phase], 1, 0) for x in (F(1, 6), F(2, 11), F(1, 5))]
        for i, j in combinations(range(11), 2):
            a, b, c = (v[i]-v[j] for v in keys)
            assert a*c >= 0 and ((b == 0) == (a == c == 0))
    print('Z3', z3.get_version_string(), subprocess.check_output(
        ['yices-smt2', '--version'], text=True).splitlines()[0])
    for site in ['open', 'lower']:
        for phase, ret in [(0, False), (1, False), (0, True)]:
            print(site, 'return' if ret else 'phase '+str(phase), query(site, phase, ret))


if __name__ == '__main__':
    main()
