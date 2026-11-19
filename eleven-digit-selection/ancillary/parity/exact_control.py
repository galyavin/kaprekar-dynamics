"""Check the reconstructed routes with whole-word subtraction at guarded endpoints."""
from fractions import Fraction as Q
from math import lcm
from affine_routes import (INITIAL_X, PARTNER, VARIABLES, ZERO, macros,
                           normalized_interval, repeat)
from cohort import gaps_form
from primitive import CENTER, step, trace

VERSION = (5, 4, 3, 2, 1)
SCALE = 3841
REPLAY_RETURNS = 30
CLOCK = (Q(901, 30000), Q(22211, 864000))
OFFSET = (Q(39), Q(209, 24))
PERIODS = (12, 24)


def physical_step(x, base):
    q = tuple((value + base * v) // 6
              for value, v in zip(x, VERSION))
    assert all(6 * p == value + base * v
               for p, value, v in zip(q, x, VERSION))
    assert base > q[0] >= q[1] >= q[2] >= q[3] >= q[4] > 0
    word = sum(p * (base ** (10 - i) - base ** i)
               for i, p in enumerate(q))
    digits = []
    for _ in range(11):
        word, d = divmod(word, base)
        digits.append(d)
    assert not word
    digits.sort(reverse=True)
    return tuple(6 * (digits[i] - digits[10 - i]) - base * v
                 for i, v in enumerate(VERSION))


def plain(state, x):
    values = tuple(f.evaluate(x) for f in state)
    assert all(v.denominator == 1 for v in values)
    return tuple(map(int, values))


def certify_partner(odd):
    base = 180000 * SCALE + odd
    state = tuple(v + odd * d for v, d in zip(VARIABLES, PARTNER))
    x = tuple(SCALE * v for v in INITIAL_X)
    clock = ZERO
    checked = 0
    projectors = []
    while any(any(f.row) for f in state):
        word, delta, div, num = macros(state)
        if not div:
            point = plain(state, x)
            state = step(state, word[0])[0]
            assert physical_step(point, base) == plain(state, x)
            checked += 1
            clock += 1
            continue
        remainder = num.evaluate(INITIAL_X) % div
        count = (num - remainder) / div + 1
        assert count.evaluate(x).denominator == 1
        n = int(count.evaluate(x))
        assert n > REPLAY_RETURNS
        projectors.append(count)
        point = plain(state, x)
        for j in (0, n - 1):
            here = tuple(v + j * d for v, d in zip(point, delta))
            for branch in word:
                target, code = step(here)
                assert code == branch
                assert target == physical_step(here, base)
                here = target
                checked += 1
        state, clock = repeat(state, clock, (word, delta, div, num), remainder)
        assert here == plain(state, x)
    endpoint = plain(state, x)
    start = tuple(v - REPLAY_RETURNS * d for v, d in zip(endpoint, delta))
    gaps = (-start[0], *(a - b for a, b in zip(start, start[1:])), start[-1])
    assert -min(gaps) > 504
    states, words, mu = trace(start)
    period = len(states) - mu
    assert period == PERIODS[odd]
    for i, a in enumerate(states):
        b = states[i + 1] if i < len(states) - 1 else states[mu]
        assert physical_step(a, base) == b
        checked += 1
    exact = clock - REPLAY_RETURNS * len(word) + mu
    assert exact.constant == OFFSET[odd]
    assert exact.evaluate(INITIAL_X) - exact.constant == 180000 * CLOCK[odd]
    assert exact.evaluate(x) == 180000 * SCALE * CLOCK[odd] + OFFSET[odd]
    lo, hi = normalized_interval(exact.row)
    print('arm', odd, 'base', base, 'endpoint', endpoint, 'terminal replay', mu,
          'period', period, 'whole-word checks', checked)
    print('exact difference clock form', exact, 'normalized range', lo, hi,
          'validated difference clock', exact.evaluate(x))
    coefficients = []
    for count in projectors:
        *chart, base_row = gaps_form(count)
        coefficients.extend(c + base_row for c in chart)
        coefficients.append(base_row)
    return lcm(6, *(c.denominator for c in coefficients))


if __name__ == '__main__':
    module = lcm(*(certify_partner(odd) for odd in range(2)))
    assert module == 3840
    print('symmetric six-gap coefficient modulus', module)
