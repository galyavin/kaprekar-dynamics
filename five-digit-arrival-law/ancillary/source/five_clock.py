"""Five-digit core clock and the all-base limiting arrival-time observable.

The proof is in five_digit_limit.md. Core errors are bounded uniformly;
limit_time is the limiting clock, not an exact finite-base arrival routine.
"""
from fractions import Fraction
from math import ceil

from five_digit import valid


def core_region(base, pair):
    u, v = pair
    return (base >= 30 and valid(base, pair) and 2*v <= base-1
            and 2*(u-v) <= base-3 and u+v-base >= -2
            and 2*u-v-base >= -2)


def margins(base, pair):
    u, v = pair
    return u+v-base, 2*u-v-base


def margin_step(pair):
    a, b = pair
    if a < -2 or b < -2:
        raise ValueError('Core margins must be at least -2.')
    if a >= 1:
        return b, a-3
    if a >= -1:
        return b, -2*a
    return b+1, 3


def core_clock(base, pair):
    """Approximate pair arrival h: -3 <= h-core_clock <= 6 on core_region."""
    if not core_region(base, pair):
        raise ValueError('Pair is outside the certified core.')
    return Fraction(2, 3)*max(margins(base, pair))


def first_continuum_step(x, y):
    """Leading digit-sort map on the open initial triangle."""
    if not 0 < y < x < 1:
        raise ValueError('Require 0 < y < x < 1.')
    if x < Fraction(1, 2):
        return 1-y, 1-x-y
    if x+y < 1:
        return 1-y, x-y
    if y < Fraction(1, 2):
        return x, x-y
    return x, x+y-1


def limit_time(x, y):
    """Explicit limit clock; rational inputs give an exact rational answer.

    One initial fold, a counted shear, then a range invariant on the strip.
    Values on branch equalities do not affect the probability distribution.
    """
    u, v = first_continuum_step(x, y)
    gap = 1-u
    count = max(0, ceil((v-Fraction(1, 2))/gap))
    v -= count*gap
    a, b = u+v-1, 2*u-v-1
    return Fraction(2, 3)*(max(0, a, b)-min(0, a, b))
