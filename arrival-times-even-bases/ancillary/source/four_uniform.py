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


def limiting_law(exponent):
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError('Require a nonnegative integer exponent.')
    n = exponent
    return {1: Q(1, 4**(n+1)),
            **{t: Q(3, 4**(n+3-t)) for t in range(2, n+3)}}


def error_bound(base):
    m, _ = base_parts(base)
    if m < 5:
        raise ValueError('This check uses odd parts >= 5.')
    return min(Q(1), Q(26, m))


def proof_error_budget(base):
    """Unclipped bound before the conservative simplification to 26/m."""
    m, _ = base_parts(base)
    if m < 5:
        raise ValueError('This check uses odd parts >= 5.')
    conditioning = Q(1, base**3)
    return (Q(16, m) + Q(9, m)/(1-conditioning) + conditioning
            + Q(1, base**2)/(1-conditioning))
