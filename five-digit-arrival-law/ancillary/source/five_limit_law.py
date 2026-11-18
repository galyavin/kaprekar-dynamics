"""Certified rational enclosures of the five-digit limiting probability law.

These are limit-law integration bounds, not finite-base approximation bounds.
Only standard-library exact arithmetic is used. See five_digit_limit.md.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations
from math import ceil, comb

ZERO = (Q(0), Q(0), Q(0))
ONE = (Q(0), Q(0), Q(1))
X = (Q(1), Q(0), Q(0))
Y = (Q(0), Q(1), Q(0))
RHO = {(1, 1): Q(120), (0, 2): Q(-120),
       (2, 1): Q(-120), (1, 2): Q(120)}


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, c):
    return tuple(x*c for x in a)


def sub(a, b):
    return add(a, scale(b, -1))


def value(form, point):
    return form[0]*point[0]+form[1]*point[1]+form[2]


def clip(poly, form):
    """Intersect a counterclockwise rational polygon with form >= 0."""
    if not poly:
        return ()
    out = []
    prev = poly[-1]
    vp = value(form, prev)
    for point in poly:
        v = value(form, point)
        if (vp < 0) != (v < 0):
            ratio = vp/(vp-v)
            out.append(tuple(p+ratio*(q-p) for p, q in zip(prev, point)))
        if v >= 0:
            out.append(point)
        prev, vp = point, v
    cleaned = []
    for p in out:
        if not cleaned or p != cleaned[-1]:
            cleaned.append(p)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return tuple(cleaned) if len(cleaned) >= 3 else ()


def integrate(poly, polynomial):
    """Polynomial area integral by Green's theorem, with exact edge integrals."""
    if not poly:
        return Q(0)
    result = Q(0)
    for start, end in zip(poly, poly[1:]+poly[:1]):
        x, y = start
        dx, dy = end[0]-x, end[1]-y
        if not dy:
            continue
        for (p, q), coefficient in polynomial.items():
            moment = sum(Q(comb(p+1, i)*comb(q, j), i+j+1)
                         * x**(p+1-i)*dx**i*y**(q-j)*dy**j
                         for i in range(p+2) for j in range(q+1))
            result += coefficient*dy*moment/(p+1)
    return result


def multiply_linear(polynomial, form):
    out = {}
    for (p, q), c in polynomial.items():
        for key, d in (((p+1, q), form[0]), ((p, q+1), form[1]),
                       ((p, q), form[2])):
            out[key] = out.get(key, Q(0))+c*d
    return {k: v for k, v in out.items() if v}


def missing_mass(cutoff):
    if not isinstance(cutoff, int) or cutoff < 1:
        raise ValueError('cutoff must be a positive integer.')
    e = Q(1, 2*(cutoff+1))
    return 20*e**2-40*e**3+32*e**5


@lru_cache(maxsize=4)
def integration_cells(cutoff):
    """Initial fold / shear cells with min(y,1-x) >= 1/(2*(cutoff+1))."""
    missing_mass(cutoff)  # Validate before constructing cells.
    e = Q(1, 2*(cutoff+1))
    poly = ((Q(0), Q(0)), (Q(1), Q(0)), (Q(1), Q(1)))
    poly = clip(clip(poly, (0, 1, -e)), (-1, 0, 1-e))
    half = scale(ONE, Q(1, 2))
    cases = (
        ((sub(half, X),), sub(ONE, Y), sub(sub(ONE, X), Y)),
        ((sub(X, half), sub(ONE, add(X, Y))), sub(ONE, Y), sub(X, Y)),
        ((sub(add(X, Y), ONE), sub(half, Y)), X, sub(X, Y)),
        ((sub(Y, half),), X, sub(add(X, Y), ONE)),
    )
    cells = []
    for guards, u, v in cases:
        initial = poly
        for guard in guards:
            initial = clip(initial, guard)
        if not initial:
            continue
        gap = sub(ONE, u)
        for n in range(cutoff+1):
            cell = initial
            if n:
                cell = clip(cell, sub(sub(v, scale(gap, n-1)), half))
            vn = sub(v, scale(gap, n))
            cell = clip(cell, sub(half, vn))
            if cell and integrate(cell, {(0, 0): Q(1)}) > 0:
                cells.append((cell, sub(add(u, vn), ONE),
                              sub(sub(scale(u, 2), vn), ONE)))
    return tuple(cells)


def cdf_bounds(t, cutoff=100):
    """Return [lower, upper] enclosing P(L <= t), both exact Fractions."""
    t = Q(t)
    if t < 0:
        return Q(0), Q(0)
    if t >= Q(1, 3):
        return Q(1), Q(1)
    if not t:
        return Q(0), Q(0)
    radius = scale(ONE, 3*t/2)
    lower = Q(0)
    for poly, a, b in integration_cells(cutoff):
        for form in (a, b, sub(a, b)):
            poly = clip(clip(poly, sub(radius, form)), add(radius, form))
            if not poly:
                break
        lower += integrate(poly, RHO)
    e = Q(1, 2*(cutoff+1))
    # On the omitted region, L >= 1/3-min(y,1-x) > 1/3-e.
    if t <= Q(1, 3)-e:
        return lower, lower
    return lower, min(Q(1), lower+missing_mass(cutoff))


def moment_bounds(order=1, cutoff=100):
    """Enclose E[L**order]. This does not assert finite-base moment convergence."""
    if not isinstance(order, int) or order < 1:
        raise ValueError('order must be a positive integer.')
    lower = Q(0)
    for poly, a, b in integration_cells(cutoff):
        for lo, mid, hi in permutations((ZERO, a, b)):
            part = clip(clip(poly, sub(mid, lo)), sub(hi, mid))
            if not part:
                continue
            clock = scale(sub(hi, lo), Q(2, 3))
            integrand = RHO
            for _ in range(order):
                integrand = multiply_linear(integrand, clock)
            lower += integrate(part, integrand)
    e = Q(1, 2*(cutoff+1))
    def omitted_moment(c):
        return sum(Q(comb(order, j), 3**(order-j))*(-c)**j
                   * (40*e**(j+2)/(j+2)-120*e**(j+3)/(j+3)
                      +160*e**(j+5)/(j+5)) for j in range(order+1))
    return lower+omitted_moment(Q(1)), lower+omitted_moment(Q(2, 3))


def cdf(t, max_cutoff=10000):
    """Exact rational F(t) for rational t, with a finite adaptive cell sum."""
    t = Q(t)
    if t <= 0:
        return Q(0)
    if t >= Q(1, 3):
        return Q(1)
    cutoff = max(1, ceil(1/(2*(Q(1, 3)-t)))-1)
    if max_cutoff is not None and cutoff > max_cutoff:
        raise ValueError(f'Exact CDF needs cutoff {cutoff}; use cdf_bounds or raise max_cutoff.')
    lo, hi = cdf_bounds(t, cutoff)
    assert lo == hi
    return lo
