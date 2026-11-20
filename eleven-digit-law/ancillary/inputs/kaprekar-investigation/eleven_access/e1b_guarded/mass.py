"""Exact original-weight bound from the remaining leading q1 anchors."""

from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "en11b"))
from entrance import IDENTITY, ONE, Poly, measure, minus, quadrature

from guarded import dot, load_rows, phases_for
from catalogue_checks import family


def anchor_weight(cutoff):
    """Task 130's six-cell decomposition for q1(K^2 q) <= cutoff."""
    rules = quadrature()
    root = Poly.root()
    q = tuple(tuple(int(i > j) for i in range(6)) for j in range(5))
    half = (Q(1, 2),) * 6
    floor = (1 - cutoff,) * 6
    anchor = root.restrict(minus(IDENTITY[i], floor) for i in (0, 5))
    weights = []
    for m in range(6):
        domain = anchor.restrict(minus(row, half) if j < m else minus(half, row)
                                 for j, row in enumerate(q))
        if domain is None:
            weights.append(Q(0))
            continue
        lows = tuple(minus(ONE, row) if j < m else row for j, row in enumerate(q))
        accepted = domain.restrict(tuple(cutoff - a - b for a, b in zip(u, v))
                                   for u, v in combinations(lows, 2))
        weights.append(measure(accepted, rules, direct=True) if accepted is not None else Q(0))
    assert tuple(weights) == tuple(reversed(weights))
    return tuple(weights)


def decimal_down(value, places):
    scale = 10 ** places
    numerator = value.numerator * scale // value.denominator
    integer, frac = divmod(numerator, scale)
    return f"{integer}.{frac:0{places}d}"


def main():
    row = load_rows()[181]
    _, tangents = phases_for(row)
    _, _, _, points = family(row)
    cutoff = Q(7, 8)
    assert all(q[0] == cutoff for q in points)
    assert all(t[0] == 0 for T in tangents for t in T)

    weights = anchor_weight(Q(5, 6))
    assert sum(weights) == Q(5124467, 80621568)
    print("Task 130 q1(K²) ≤ 5/6 control (exact):", sum(weights))
    weights = anchor_weight(cutoff)
    print("Six sorting-split weights at 7/8 (exact):", [str(w) for w in weights])
    upper = sum(weights)
    print("Original weight of q1(K²) ≤ 7/8 (exact):", upper)
    bound = 1 - upper
    print("m(D) lower bound (exact):", bound)
    print("m(D) lower bound in percent (rounded down):", decimal_down(100*bound, 10))


if __name__ == "__main__":
    main()
