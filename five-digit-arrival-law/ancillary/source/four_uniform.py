"""Uniform four-slot arrival approximation as the odd part of the base grows.

The proof is in arithmetic_scaling.md. All probabilities and bounds are exact
Fractions. The generic pair/actual clocks reject the exceptional congruences;
the distribution theorem includes those inputs through a total-variation bound.
"""
from fractions import Fraction as Q


def base_parts(base):
    if not isinstance(base, int) or base < 2:
        raise ValueError('Require integer base >= 2.')
    odd, exponent = base, 0
    while odd % 2 == 0:
        odd //= 2
        exponent += 1
    return odd, exponent


def _valuation_time(z, n):
    if z == 0:
        return 1  # v_2(0)=infinity in the auxiliary clock.
    a = (z & -z).bit_length() - 1
    return 0 if a == n else 1 if a > n else n - a + 1


def regular_pair(base, pair):
    m, _ = base_parts(base)
    u, v = pair
    return (m >= 5 and 0 <= v <= u < base and u > 0
            and all(z % m for z in (u, v, u-v, u+v)))


def surrogate_time(base, pair):
    """Valuation prediction h+1; not exact on exceptional pairs or cycles."""
    _, n = base_parts(base)
    u, v = pair
    if not 0 <= v <= u < base:
        raise ValueError('Require 0 <= v <= u < B.')
    return 1 + max(_valuation_time(u, n), _valuation_time(v, n))


def generic_pair_arrival(base, pair):
    """Exact first pair-cycle arrival outside the invariant exceptional set."""
    if not regular_pair(base, pair):
        raise ValueError('The pair meets an excluded odd-core congruence.')
    return surrogate_time(base, pair) - 1


def generic_actual_arrival(base, digits):
    """Exact actual-input time on regular pairs, including cyclic inputs."""
    if len(digits) != 4 or any(not isinstance(d, int) or not 0 <= d < base for d in digits):
        raise ValueError('Require four base-B digits.')
    a, b, c, d = sorted(digits, reverse=True)
    pair = a-d, b-c
    h = generic_pair_arrival(base, pair)
    if h:
        return h+1
    # K(u,v)=(B^3-1)u+(B^2-B)v is injective on pair states.
    number = 0
    for digit in digits:
        number = base*number+digit
    u, remainder = divmod(number, base**3-1)
    v, rem = divmod(remainder, base**2-base)
    _, n = base_parts(base)
    scale = 1 << n
    on_cycle = (rem == 0 and 0 < v < u < base
                and u % (2*scale) == v % (2*scale) == scale)
    return 0 if on_cycle else 1


def limiting_law(exponent):
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError('Require a nonnegative integer exponent.')
    n = exponent
    return {1: Q(1, 4**(n+1)),
            **{t: Q(3, 4**(n+3-t)) for t in range(2, n+3)}}


def error_bound(base):
    m, _ = base_parts(base)
    if m < 5:
        raise ValueError('The uniform theorem requires odd part >= 5.')
    return min(Q(1), Q(26, m))


def proof_error_budget(base):
    """Unclipped bound before the conservative simplification to 26/m."""
    m, _ = base_parts(base)
    if m < 5:
        raise ValueError('The uniform theorem requires odd part >= 5.')
    conditioning = Q(1, base**3)
    return (Q(16, m) + Q(9, m)/(1-conditioning) + conditioning
            + Q(1, base**2)/(1-conditioning))


def random_base_probability(time):
    """Limiting law: base uniform in {2,...,X}, then input uniform in base."""
    if not isinstance(time, int) or time < 0:
        raise ValueError('Require a nonnegative integer time.')
    if time == 0:
        return Q(0)
    return Q(1, 7) if time == 1 else Q(3, 7*2**(time-2))
