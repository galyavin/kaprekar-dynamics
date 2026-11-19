"""Definitions for the physical-lattice centered clock.

The finite-horizon claim requires the complete independent certificate
audit. Actual finite-base conclusions additionally require paired guards.
"""
from math import comb
from odd_base_covariance import centered_step, centered, in_region

HORIZON = 32
THRESHOLD = 120
GROWTH = 12
ABSORBING = THRESHOLD+HORIZON*GROWTH
STATE_BOUND = comb(ABSORBING+5, 5)
V = (5, 4, 3, 2, 1)


def require(ok, label):
    if not ok: raise ArithmeticError(label)


def gaps(x):
    return (-x[0], *(a-b for a, b in zip(x, x[1:])), x[-1])


def residual(x):
    return -min(gaps(x))


def base_class(x):
    if len(x) != 5 or any(type(v) is not int for v in x):
        raise ValueError('Five integral centered coordinates required')
    beta = x[0] % 6
    if any((v-i*beta) % 6 for i, v in enumerate(x, 1)):
        raise ValueError('The physical centered lattice is required')
    return beta


def preparation_bound(x):
    base_class(x)
    return HORIZON*max(0, (residual(x)-THRESHOLD+5)//6)


def lyapunov(x):
    """Finite-lookahead one-step rank, proved after full clock acceptance."""
    base_class(x)
    best = HORIZON*residual(x)
    for j in range(1, HORIZON):
        x = centered_step(x)
        best = min(best, HORIZON*residual(x)+6*j)
    return best


def central_source(B, q):
    return (type(B) is int and B >= 3534 and len(q) == 5 and
            all(type(v) is int for v in q) and in_region(B, q) and
            min(B-q[0], *[a-b for a, b in zip(q, q[1:])], q[-1]) >= (B+7)//7+64)


def actual_bound(B, q):
    if not central_source(B, q):
        raise ValueError('Source is outside the declared central region')
    return preparation_bound(centered(B, q))+STATE_BOUND


def long_family(B, j=0):
    """A central physical source with a certified eight-step return block."""
    if type(B) is not int or B < 12000 or type(j) is not int or not 0 <= j <= B//6000:
        raise ValueError('Require B>=12000 and integer 0<=j<=floor(B/6000)')
    from eleven_two_cycle_clock import pair
    t = B//6000
    a = 18*t+(5*B)%6-18*j
    d = -108*t+(-2*B)%6+54*j
    return pair(B, a, d, (-4, -10, 5))
