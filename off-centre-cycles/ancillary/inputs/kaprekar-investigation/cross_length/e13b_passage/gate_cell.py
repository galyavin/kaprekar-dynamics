#!/usr/bin/env python3
"""Exact two-return gate certificate on a full source lattice in each class."""

from fractions import Fraction
from math import factorial

from premise import KICK, PERMUTATION, add, advance, centered, differences


DRIFT = (-42, -35, -28, -21, -14, -7)
POST_DRIFT = (-84, -70, -7, 7, -28, -14)
WORD = (12, 0, 0, 0, 0, 0) * 2
ZERO = (0,) * 6
D = 32
BASE_SLOPE = 140000
MAX_BASE_RESIDUAL = 7 * BASE_SLOPE - 1
GATE_DIRECTION = (3450, 2595, 613, -613, -544, -531)
GATE_COLUMNS = tuple(tuple(7 * (i == j) - 7 * (j == 2 and i == 3) for i in range(6))
                     for j in (0, 1, 2, 4, 5)) + (ZERO,)


def gaps(x):
    return tuple(a - b for a, b in zip((0, *x), (*x, 0)))


def slots(x, mask, shift):
    values = [-shift]
    for i in range(6):
        left = x[i] - shift * (i == 5)
        right = -x[5 - i] - shift * (i < 5)
        values.extend((right, left) if mask >> i & 1 else (left, right))
    return values


def update(x, mask, shift):
    values = slots(x, mask, shift)
    return tuple(values[i] - values[12 - i] for i in range(6))


def guards(base, x, mask, shift):
    values = slots(x, mask, shift)
    return (tuple(base + gap - shift * (i in (0, 6)) for i, gap in enumerate(gaps(x)))
            + tuple(base + values[2*i] - values[2*i + 1] for i in range(6))
            + tuple(values[2*i + 1] - values[2*i + 2] - shift * (mask >> i & 1) for i in range(6)))


def inverse_five(x, shift):
    g = gaps(x)
    for _ in range(5):
        g = add(g, KICK, -shift)
        g = tuple(g[PERMUTATION.index(i)] for i in range(7))
    return tuple(-sum(g[:i]) for i in range(1, 7))


def check(offset, slope, perturbations, start_n=0, strict_slope=False):
    worst_slope = slope - Fraction(sum(map(abs, perturbations)), D)
    assert worst_slope >= 0 and offset + start_n * worst_slope >= 0
    if strict_slope:
        assert worst_slope > 0


def check_forms(form, base, family, start_n=0, strict_slope=False):
    f0 = form(base, family[0], 7)
    fn = form(BASE_SLOPE, family[1], 0)
    fz = [form(0, vector, 0) for vector in family[2:]]
    for offset, slope, *perturbations in zip(f0, fn, *fz):
        check(offset, slope, perturbations, start_n, strict_slope)


def proof_word(base, family, word, off_anchor_gate=False, strict_slope=False):
    for mask in word:
        check_forms(lambda b, x, s: guards(b, x, mask, s), base, family, strict_slope=strict_slope)
        if off_anchor_gate:
            check_forms(lambda b, x, s: tuple(b + 64*g for g in gaps(x)[1:]), base, family)
        next_family = tuple(update(x, mask, 7 * (i == 0)) for i, x in enumerate(family))
        assert next_family[0][0] == family[0][0] - 7
        assert all(new[0] == old[0] for new, old in zip(next_family[1:], family[1:]))
        family = next_family
    return family


def shift_family(family, drift, constant, slope, last=0):
    return (add(family[0], drift, constant), add(family[1], drift, slope),
            *family[2:-1], add(family[-1], drift, last))


def prove(base, exit_seed):
    assert all((x - i*base) % 7 == 0 for i, x in enumerate(exit_seed, 1))
    assert exit_seed[2] + exit_seed[3] == -14
    assert all(x[2] + x[3] == 0 for x in (GATE_DIRECTION, *GATE_COLUMNS))
    exit_family = (exit_seed, GATE_DIRECTION, *GATE_COLUMNS)
    before_drift = tuple(inverse_five(x, 7 * (i == 0)) for i, x in enumerate(exit_family))
    source = shift_family(before_drift, DRIFT, -3177, -3177, -1)

    # All complete chamber copies and the five phases at the first failure.
    assert proof_word(base, source, (0,) * 6, strict_slope=True) == shift_family(source, DRIFT, 1, 0)
    last_full = shift_family(source, DRIFT, 3176, 3177, 1)
    proof_word(base, last_full, (0,) * 6)
    assert proof_word(base, before_drift, (0,) * 5) == exit_family
    assert guards(base, exit_seed, 0, 7)[15:17] == (-7, -7)

    # Both ends of the post-exit word. Affinity interpolates all copies.
    assert proof_word(base, exit_family, WORD, True) == shift_family(exit_family, POST_DRIFT, 1, 0)
    last_post = shift_family(exit_family, POST_DRIFT, 15, 16)
    arrived = proof_word(base, last_post, WORD, True)
    assert arrived == shift_family(exit_family, POST_DRIFT, 16, 16)
    check_forms(lambda b, x, s: tuple(b + 64*g for g in gaps(x)), base, arrived)

    # The exact first normalized crossing is within the post-exit word.
    check_forms(lambda b, x, s: (64*x[0] - b - (MAX_BASE_RESIDUAL + 1) * (s == 7),),
                base, exit_family, start_n=21)

    # Each physical source gap is at least B/500. Replacing the nominal
    # base by the actual base only increases this guard.
    check_forms(lambda b, x, s: tuple(500*(b+g) - 7*b for g in gaps(x)), base, source)

    # The zero-perturbation finite control uses whole-word subtraction.
    q = differences(base, gaps(source[0]))
    t_exit = 6 * 3177 + 5
    t_hit = None
    for t in range(t_exit + 193):
        g = centered(base, q)
        if t == t_exit:
            assert g == gaps(exit_seed)
        if t_hit is None and 64 * -min(g) <= base:
            t_hit = t
        if t == t_exit + 192:
            assert g == gaps(arrived[0])
        q = advance(base, q)
    assert t_hit == t_exit + (64 * exit_seed[0] - base + 447) // 448
    print(f'class {base % 7}: source X={source[0]}, exit X={exit_seed}, '
          f'first-exit step {t_exit}, first-gate step {t_hit}, post-word drift {POST_DRIFT}')


if __name__ == '__main__':
    prove(140000, (3451, 2597, 609, -623, -546, -532))
    prove(140006, (3450, 2595, 606, -620, -544, -531))
    print(f'Both classes: every B >= 4000000, original mass to first normalized gate >= '
          f'{Fraction(factorial(13), 500 * 1000**6 * (280000 * D)**6)} (inclusive).')
