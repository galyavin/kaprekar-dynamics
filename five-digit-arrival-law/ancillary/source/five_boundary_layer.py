"""Exact high excursions and the conditional five-slot square-root law.

Uses only the analytic core theorem from five_digit_limit.md. The exact exit
routine checks its domain and returns a certified clock interval. It makes no
claim of global entry for arbitrary pairs. See arithmetic_scaling.md.
"""
from fractions import Fraction as Q
from math import isqrt, sqrt

from five_clock import core_region, core_clock
from five_digit import transition, valid


def exit_data(base, pair):
    """Exact off-diagonal high segment, followed by at most one endpoint step.

    Raises ValueError outside the checked domain. The returned pair-time
    interval is justified only when the reached core state is transient.
    """
    if base < 30 or not valid(base, pair):
        raise ValueError('Require a valid pair in base >= 30.')
    u, v = pair
    r, threshold = base-u, (base+1)//2
    if r < 2 or v <= threshold or u <= v:
        raise ValueError('Require r>=2 and an off-diagonal high start.')
    # D(j)=j(j+2r-3)/2 is increasing for j>=0.
    target, coefficient = v-threshold, 2*r-3
    n = max(0, (isqrt(coefficient**2+8*target)-coefficient)//2)
    def drift(j):
        return j*(j+coefficient)//2
    if drift(n) < target:
        n += 1
    if not drift(n-1) < target <= drift(n):
        raise ArithmeticError('Invalid quadratic-exit inversion.')
    # u_j-v_j=(u-v)+j(r-2)+j(j-1)/2 >= u-v > 0.
    q = u-n, v-drift(n)
    R, overshoot = r+n, threshold-q[1]
    if not 0 <= overshoot <= R-3:
        raise ArithmeticError('Overshoot outside its exact integer block.')
    extra = 0
    if overshoot == 0:
        q = transition(base, q)
        extra = 1
    if not core_region(base, q):
        raise ValueError('High exit does not reach the certified core directly.')
    phi = core_clock(base, q)
    if phi-3 <= 0:
        raise ValueError('Endpoint is not certified transient by the core clock.')
    center = n+Q(2, 3)*(Q(base, 2)-R-min(overshoot, R-overshoot))
    return {'high_steps': n, 'endpoint_steps': extra, 'core_pair': q,
            'R': R, 'overshoot': overshoot, 'phase': Q(overshoot, R-2),
            'pair_lower': n+extra+phi-3,
            'pair_upper': n+extra+phi+6,
            'actual_lower': n+extra+phi-2,
            'actual_upper': n+extra+phi+7,
            'center_clock': center}


def kernel_interval(lam, z):
    if not lam > 0 or not Q(1, 2) < z < 1:
        raise ValueError('Require lambda>0 and 1/2<z<1.')
    s = sqrt(lam*lam+2*z-1)
    return -lam, s/3-lam


def kernel_cdf(time, lam, z):
    lo, hi = kernel_interval(lam, z)
    return max(0.0, min(1.0, (time-lo)/(hi-lo)))


def band_mass_constant(a, b, c, d):
    """Exact M when the four endpoints are Fraction-compatible rationals."""
    a, b, c, d = map(Q, (a, b, c, d))
    if not 0 < a < b or not Q(1, 2) < c < d < 1:
        raise ValueError('Require 0<a<b and 1/2<c<d<1.')
    return 60*(b*b-a*a)*((d*d-c*c)/2-(d**3-c**3)/3)


def uniform_band_threshold(a, b, c, d):
    """Conservative explicit B threshold for the compact-band guards.

    Require sqrt(B)>=max(16/a,16(b+2)/min(c-1/2,1-d),32).
    The expression is returned squared and rounded upward, with exact rationals.
    """
    a, b, c, d = map(Q, (a, b, c, d))
    band_mass_constant(a, b, c, d)
    epsilon = min(c-Q(1, 2), 1-d)
    bound = max(Q(32), 16/a, 16*(b+2)/epsilon)**2
    return (bound.numerator+bound.denominator-1)//bound.denominator
