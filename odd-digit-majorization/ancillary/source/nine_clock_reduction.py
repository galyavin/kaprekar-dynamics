"""Dispersion, residual tube, zero-return cover and sparse-exception time bound."""
from math import comb
from odd_gap import gaps, valid


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def dispersion(base, q):
    g = sorted(gaps(base, q))
    return sum((2*i-4)*a for i, a in enumerate(g))


def terminal_residual(base, q):
    g = gaps(base, q)
    m = min(g)
    return base-4*m-min(max(g), 2*m)


def return_cover(base, q):
    """A sufficient cover of zero states that ever return to zero.

    Membership does not assert that the state returns. For base >= 40 the
    cover contains at most 308*base states, including ranks one and two.
    """
    if not valid(base, q) or len(q) != 4 or q[3] != 0:
        raise ValueError('Expected a nonconstant zero-inner nine-slot state.')
    u, v, w, _ = q
    if v == 0:
        return True
    if w == 0:
        return (abs(u-v) <= 3 or abs(u+v-base) <= 3
                or abs(2*u-base) <= 4 or abs(2*v-base) <= 4
                or base-u <= 3)
    return (
        (v == w and (abs(u-v) <= 2 or abs(u+v-base) <= 2
                     or abs(2*v-base) <= 3))
        or (u == v and (abs(u+w-base) <= 2 or abs(2*u-base) <= 3))
        or (abs(2*u-base) <= 3 and abs(2*v-base) <= 3)
        or (base-u <= 3 and (abs(2*v-base) <= 3 or v == w
                            or abs(v+w-base) <= 3))
        or (abs(4*u-3*base) <= 8
            and any(abs(4*x-j*base) <= 12
                    for x in (v, w) for j in (1, 2, 3)))
    )


def zero_preparation_bound(base):
    if base < 40:
        raise ValueError('The proved zero preparation bound requires B >= 40.')
    return 1232*base+7


def terminal_tube_count_bound(base, threshold):
    if base < 2 or threshold < 0 or not isinstance(threshold, int):
        raise ValueError('Expected integer B >= 2 and threshold >= 0.')
    return 5*(base+1)*comb(threshold+4, 4)


def conditional_pair_clock_bound(base, horizon, threshold):
    """Conditional on the descent lemma for the supplied horizon and threshold.

    Outside residual <= threshold, assume some positive iterate of index
    1..horizon has smaller dispersion. For horizon=9 and threshold=100,
    the tree certificate is reconstructed by nine_clock_audit.py.
    """
    if not isinstance(horizon, int) or horizon < 1:
        raise ValueError('Expected a positive integer horizon.')
    return (zero_preparation_bound(base)+4*horizon*base
            +(8*horizon+1)*terminal_tube_count_bound(base, threshold)
            +horizon)
