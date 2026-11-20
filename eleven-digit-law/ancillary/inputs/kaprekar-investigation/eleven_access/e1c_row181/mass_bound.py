"""Exact original weight of a three-step necessary event for row 181."""

from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from subprocess import PIPE, Popen
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / name) for name in ("en11b", "e1b_guarded")]

from anchor import GAPS, gap_values
from entrance import IDENTITY, Poly, image_domain, integrate, minus, quadrature
from guarded import dot, family, load_rows, phases_for
from mass import decimal_down


def verify_limit_gaps():
    row = load_rows()[181]
    _, tangents = phases_for(row)
    _, _, _, points = family(row)
    for q, (tangent,) in zip(points, tangents):
        fixed = sorted(gap for gap, normal in zip(gap_values(q), GAPS)
                       if dot(normal, tangent) == 0)
        assert fixed == [Q(1, 8)] * 3 + [Q(1, 4)]
        assert q[0] == Q(7, 8) and tangent[0] == 0


def event_cells():
    indices = tuple(range(6))
    base = Poly.root().restrict(minus(IDENTITY[i], (Q(1, 8),) * 6) for i in (0, 5))
    base = image_domain(base, IDENTITY)
    for smallest in combinations(indices, 3):
        other = tuple(i for i in indices if i not in smallest)
        small_cut = tuple(Q(3, 8) - int(i in smallest) for i in indices)
        small = base.restrict((small_cut, *(minus(IDENTITY[j], IDENTITY[i])
                                           for i in smallest for j in other)))
        if small is None:
            continue
        for largest in other:
            cell = small.restrict((minus(IDENTITY[largest], (Q(1, 4),) * 6),
                                   *(minus(IDENTITY[largest], IDENTITY[i])
                                     for i in other if i != largest)))
            if cell is not None:
                yield smallest, largest, cell


def compute(evaluator):
    rules = quadrature()
    upper = Q(0)
    for small, largest, cell in event_cells():
        weight, branches, distinct = integrate(cell, 3, rules, evaluator=evaluator)
        upper += weight
        print(f"{small} {largest}: {weight} ({branches} paths, {distinct} domains)", flush=True)
    bound = 1 - Q(5124467, 80621568) - upper
    assert 0 < upper < Q(1, 100)
    units = 100 * upper * 10**10
    rounded_up = Q(-(-units.numerator // units.denominator), 10**10)
    print("Row 181 upper mass (exact necessary event):", upper)
    print("Row 181 upper mass in percent (rounded up):", decimal_down(rounded_up, 10))
    print("m(D) lower bound (exact):", bound)
    print("m(D) lower bound in percent (rounded down):", decimal_down(100 * bound, 10))


if __name__ == "__main__":
    verify_limit_gaps()
    if len(sys.argv) == 1:
        compute(None)
    else:
        with Popen([str(Path(sys.argv[1]).resolve())], stdin=PIPE, stdout=PIPE, text=True) as evaluator:
            compute(evaluator)
            evaluator.stdin.close()
            assert evaluator.wait() == 0
