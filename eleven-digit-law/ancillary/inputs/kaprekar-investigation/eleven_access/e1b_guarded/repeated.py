"""Exact certificates for strict periodic sorting words in the local map."""

from fractions import Fraction as Q
from functools import lru_cache

import sympy as sp

from guarded import I, N, difference, multiply, pull, row_times
from catalogue_checks import power


X = sp.Symbol("x")


@lru_cache(maxsize=4096)
def polynomial_data(W):
    coeff = tuple(map(int, sp.Matrix(W).charpoly(X).all_coeffs()))
    poly = sp.Poly.from_list(coeff, X)
    # adj(xI-W) = sum_j x^(N-1-j) C_j by Cayley-Hamilton.
    blocks = [I]
    for c in coeff[1:N]:
        current = multiply(W, blocks[-1])
        blocks.append(tuple(tuple(a + c * b for a, b in zip(row, unit))
                            for row, unit in zip(current, I)))
    return poly, tuple(blocks)


def polynomial_interval(coeff, low, high):
    """Enclose a polynomial on the rational interval by interval Horner."""
    lo = hi = 0
    for c in coeff:
        products = (lo * low, lo * high, hi * low, hi * high)
        lo, hi = min(products) + c, max(products) + c
    return lo, hi


def adjugate_interval(guards, blocks, low, high):
    """Return a signed adjugate column strictly in the guard cone, if certified."""
    # Work with columns of g*adj(xI-W), one polynomial per guard and column.
    products = [tuple(row_times(g, C) for C in blocks) for g in guards]
    for j in range(N):
        intervals = [polynomial_interval([row[j] for row in item], low, high)
                     for item in products]
        if all(left > 0 for left, _ in intervals):
            return j, 1
        if all(right < 0 for _, right in intervals):
            return j, -1
    return None


def repeatable(cones, W, guards):
    """Return exact strict-ray certificates. Intervals isolate real eigenvalues."""
    if power(W, 120) == I:
        return ()
    # Any fixed eigenray must at least allow the same strict word twice.
    if not cones.feasible((*guards, *pull(guards, W))):
        return ()
    poly, blocks = polynomial_data(W)
    result = []
    # Irreducible factors allow refinement even for repeated spectrum.
    for root_factor, _ in poly.factor_list()[1]:
        roots = root_factor.intervals(eps=sp.Rational(1, 10**12), inf=0)
        for (a, b), multiplicity in roots:
            assert multiplicity == 1
            a, b = Q(a), Q(b)
            if a == b:
                if a == 1:
                    conditions = tuple(difference(row, unit) for row, unit in zip(W, I))
                    if cones.check(guards, conditions):
                        result.append({"type": "fixed", "eigenvalue": "1"})
                continue
            # Root at 0 or 1 is rational, so an irrational root is interior.
            if b <= 0:
                continue
            if a <= 0 or a <= 1 <= b:
                continue
            column = adjugate_interval(guards, blocks, a, b)
            if column is not None:
                j, sign = column
                result.append({"type": "contracting" if b < 1 else "expanding",
                               "interval": [str(a), str(b)], "column": j, "sign": sign})
    return tuple(result)


def invariant_repeat(cones, W, guards, cap=5):
    """Certify a nonempty open cone trapped in the word by rational implication."""
    if power(W, 120) == I:
        return None
    domain = guards
    for depth in range(1, cap + 1):
        next_guards = pull(domain, W)
        if cones.contains(next_guards, domain):
            return {"type": "invariant cone", "cone": domain, "depth": depth}
        refined = (*domain, *next_guards)
        if not cones.feasible(refined):
            return None
        domain = cones.facets(refined)
    return None


def replay_ray(cones, W, guards, ray):
    assert power(W, 120) != I
    if ray["type"] == "invariant cone":
        cone = tuple(map(tuple, ray["cone"]))
        assert cones.feasible(cone) and cones.contains(guards, cone)
        assert cones.contains(pull(cone, W), cone)
    elif ray["type"] == "fixed":
        assert ray["eigenvalue"] == "1"
        equations = tuple(difference(row, unit) for row, unit in zip(W, I))
        assert cones.check(guards, equations)
    else:
        a, b = map(Q, ray["interval"])
        assert 0 < a < b and (b < 1 if ray["type"] == "contracting" else a > 1)
        poly, blocks = polynomial_data(W)
        assert poly.count_roots(sp.Rational(a), sp.Rational(b)) == 1
        assert adjugate_interval(guards, blocks, a, b) == (ray["column"], ray["sign"])
