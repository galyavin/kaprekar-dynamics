"""Dual-solver checks for exact Two-edge local sorting and return formulas."""

import argparse
import hashlib
from itertools import combinations
from pathlib import Path
import subprocess

import z3

from preimages import borrow, step, whole_step


# Each representative specifies only strict order between leading groups.
# The open edge has one such order. The lower endpoint adds the central tie.
PHASES = {
    'open': (120, (98,76,54,27,22), (98,76,71,49,22)),
    'lower': (12, (10,8,6,3,2), (10,8,7,5,2)),
}


def largest(a, b):
    return z3.If(a >= b, a, b)


def smallest(a, b):
    return z3.If(a <= b, a, b)


def positive(x):
    return largest(0, x)


def descending(entries):
    result = []
    for entry in entries:
        for i, value in enumerate(result):
            result[i], entry = largest(value, entry), smallest(value, entry)
        result.append(entry)
    return result


def borrowed(q, base):
    reflection = tuple(base-x-(i > 0) for i, x in enumerate(q))
    direct = (*q[:-1], q[-1]-1)
    return (base-1, *direct, *reflection)


def offsets(q, base):
    return tuple(6*x-(5-i)*base for i, x in enumerate(q))


def fold(digits):
    return tuple(digits[i]-digits[10-i] for i in range(5))


def local_map(solver, q, base, macro, scale):
    digits = borrowed(q, base)
    weights = (scale, *macro, *(scale-x for x in macro))
    groups = {w: [] for w in weights}
    for w, digit in zip(weights, digits):
        groups[w].append(digit)
    for i, j in combinations(range(11), 2):
        if weights[i] != weights[j]:
            greater, smaller = (i, j) if weights[i] > weights[j] else (j, i)
            solver.add(digits[greater] > digits[smaller])
    solver.add(q[0] < base, q[-1] > 0, *[a >= b for a,b in zip(q,q[1:])])
    return fold([d for w in sorted(groups, reverse=True) for d in descending(groups[w])])


def paired(q, base):
    digits = borrowed(q, base)
    local = [digits[0]]
    for i in range(5):
        local.extend(descending((digits[i+1], digits[10-i])))
    return fold(local)


def assume_distinct_values(solver, exact, proposed):
    solver.add(z3.Or(*(x != y for x, y in zip(exact, proposed))))


def swap_query(name, phase, base):
    scale, canon, other = PHASES[name]
    q = z3.Ints('q1 q2 q3 q4 q5')
    solver = z3.SolverFor('QF_LIA')
    raw = local_map(solver, q, base, (canon, other)[phase], scale)
    expected = list(paired(q, base))
    if phase:
        expected[3] += positive(q[2]+q[3]-base+1)
    assume_distinct_values(solver, raw, expected)
    return solver


def return_query(name, base):
    scale, canon, other = PHASES[name]
    q = z3.Ints('q1 q2 q3 q4 q5')
    middle = z3.Ints('v1 v2 v3 v4 v5')
    a, w, p, s, r = z3.Ints('a w p s r')
    solver = z3.SolverFor('QF_LIA')
    source = (-a, -2*a+6*p, -3*a+6*s, 6*w-4*a, a+6*r)
    solver.add(*(x == y for x,y in zip(offsets(q, base), source)))
    first = local_map(solver, q, base, canon, scale)
    solver.add(*(x == y for x,y in zip(middle, first)))
    raw = offsets(local_map(solver, middle, base, other, scale), base)

    A = largest(-1, -r)
    J = positive(-r-1)-positive(r-1)
    K = positive(r+1)
    C = positive(-a+2*s+1)
    t = s-p-C
    D = largest(A-1, t)
    E = largest(A, t-1)-largest(-t-1, -A)
    F = smallest(A, t-1)+J+1
    delta = positive(middle[2]+middle[3]-base+1)
    solver.add(w-a+p <= -1, a-2*w-2*K >= -1, w-J-p-1 <= 0)
    longitudinal = a-6*D
    expected = (
        -longitudinal,
        -2*longitudinal+6*(E-2*D),
        -3*longitudinal+6*(F-3*D),
        6*(w+J-p-4*D+delta)-4*longitudinal,
        longitudinal+6*(-p-K+D),
    )
    minimum_change = raw[3]-raw[4]-source[3]+source[4]
    maximum_change = raw[2]-raw[3]-source[2]+source[3]
    solver.add(z3.Or(
        *(x != y for x,y in zip(raw, expected)),
        minimum_change != 6*(r+J+K+delta),
        maximum_change != -6*(C+positive(t-A-1)+delta),
    ))
    return solver


def check_query(name, solver, output):
    solver.set(timeout=30000)
    transcript = '(set-logic QF_LIA)\n'+solver.to_smt2()
    result = str(solver.check())
    check = subprocess.run(['yices-smt2', '--timeout=30'], input=transcript,
                           capture_output=True, text=True, check=True, timeout=35)
    assert result == 'unsat' == check.stdout.strip(), (name,result,check.stdout)
    if output:
        (output/f'{name}.smt2').write_text(transcript)
    print(name, result, check.stdout.strip(), hashlib.sha256(transcript.encode()).hexdigest())


def guards(q, base):
    pairs = [(q[i]-(i == len(q)-1), base-q[-1-i]-(i != len(q)-1))
             for i in range(len(q))]
    return tuple(min(a)-max(b) for a,b in zip(pairs,pairs[1:]))


def full_word_controls():
    for t in (2, 3, 20, 100, 1003):
        base = 12*t
        q = (10*t,8*t,6*t,3*t-1,2*t)
        first = whole_step(q, base)
        assert first == step(q, base) == (10*t,8*t,7*t,5*t+1,2*t+1)
        assert guards(q, base) == (t-1,2*t,2*t,t-1)
        assert guards(first, base) == (2*t-2,-2,-2,2*t-1)
        assert max(0, first[2]+first[3]-base+1) == 2
        second = whole_step(first, base)
        raw = borrow(first, base)
        intended = [raw[0]]
        for i in range(5):
            intended.extend(sorted((raw[i+1],raw[10-i]),reverse=True))
        predicted = list(fold(intended))
        predicted[3] += 2
        assert second == tuple(predicted)

    # Internal five-digit switches in the known positive-pairing family.
    def five(q, base):
        word = (*q, 0, 0, 0)
        high = low = 0
        for left,right in zip(word,reversed(word)):
            high = high*base + left
            low = low*base + right
        n = high-low
        digits = []
        for _ in word:
            n, digit = divmod(n, base)
            digits.append(digit)
        assert n == 0
        digits.sort(reverse=True)
        return (digits[0]-digits[4],digits[1]-digits[3])

    for t in (20, 100):
        B = 60*t
        q = (42*t,18*t-1)
        signs = []
        margins = []
        for _ in range(12):
            assert min(guards(q, B)) > 0 and q[-1] > 0
            margins.append((sum(q)-B, 2*q[0]-q[1]-B))
            signs.append((sum(q)-B > 0)-(sum(q)-B < 0))
            q = five(q, B)
        M = 6*t+1
        assert margins[:5] == [(-1,M), (M,2), (2,M-3), (M-3,-1), (-1,M-6)]
        assert signs[0] == -1 and signs[1] == 1
        assert signs.count(-1) > 1 and signs.count(1) > 1
    print('Independent full-word lower-endpoint and five-digit controls passed.')


def main():
    if not __debug__:
        raise RuntimeError('Assertions are required for this certificate.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
    print('Z3', z3.get_version_string())
    subprocess.run(['yices-smt2','--version'], check=True)
    base = z3.Int('B')
    for name in PHASES:
        for phase in range(2):
            check_query(f'{name}-phase-{phase}', swap_query(name, phase, base), args.output)
        check_query(f'{name}-return', return_query(name, base), args.output)
    full_word_controls()


if __name__ == '__main__':
    main()
