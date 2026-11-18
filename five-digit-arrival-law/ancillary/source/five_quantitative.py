"""Sharp-order maximum witnesses and a conservative finite-base CDF bound."""
from fractions import Fraction as Q
from math import isqrt


def maximum_witness(base):
    """(pair, high-shear length) for the matching sqrt(B)/3 lower bound."""
    if not isinstance(base, int) or base < 30:
        raise ValueError('Require integer B >= 30.')
    v_end = (base - 1) // 2
    capacity = base - 3 - v_end
    n = (1 + isqrt(1 + 8 * capacity)) // 2
    assert n * (n - 1) // 2 <= capacity < n * (n + 1) // 2
    return (base - 1, v_end + n * (n - 1) // 2), n


def finite_cdf_error_bound(base, cutoff=None):
    """Uniform Kolmogorov bound for actual five-slot input times.

    Deliberately conservative. The bound is often trivial at practical bases;
    the theorem supplies an effective O((log B)^-2) rate, not a fitted rate.
    """
    if not isinstance(base, int) or base < 30:
        raise ValueError('Require integer B >= 30.')
    if cutoff is None:
        cutoff = max(1, (base.bit_length()-1)//2)
    if not isinstance(cutoff, int) or cutoff < 1:
        raise ValueError('Require a positive integer cutoff.')
    eps = Q(1, 2*(cutoff+1))
    omitted = 20*eps**2 - 40*eps**3 + 32*eps**5
    geometric = Q(100000*(cutoff+4)**2 *
                  (2**(cutoff+5)+cutoff+10), base)
    conditioning = Q(1, base**4-1)
    return min(Q(1), omitted+geometric+conditioning)
