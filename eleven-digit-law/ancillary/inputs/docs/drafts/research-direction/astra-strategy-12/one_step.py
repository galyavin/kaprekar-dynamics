"""Exact original mass entering Task 097's regular chamber by one leading step.

Run with python3 -B one_step.py. The archived geometry is read-only.
"""
from fractions import Fraction
from itertools import product
from math import factorial, lcm, prod
from pathlib import Path
import json
import sys

GRAPH = Path(__file__).resolve().parents[4] / "kaprekar-investigation/eleven_access/stage_graph"
sys.dont_write_bytecode = True
sys.path.insert(0, str(GRAPH))
from pi_geometry import Poly, det, minus, within_interlevel


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def quadrature():
    # Degree-seven simplex formula. The degree-six identity actually used below
    # is verified on every monomial, up to permutation of barycentric variables.
    rules = []
    for i in range(4):
        size = 12 - 2 * i
        weight = Fraction((-1)**i * size**7, 64 * factorial(i) * factorial(12-i))
        nodes = [tuple(2*x+1 for x in exponents) for exponents in compositions(3-i, 6)]
        rules.append((size, weight, nodes))
    for powers in {tuple(sorted(p)) for p in compositions(6, 6)}:
        integral = sum(weight * Fraction(sum(prod(x**p for x, p in zip(node, powers))
                                               for node in nodes), size**6)
                       for size, weight, nodes in rules)
        assert integral == Fraction(prod(factorial(p) for p in powers), factorial(11))
    return rules


MASKS = list(product((False, True), repeat=5))


def lows_from_gaps(g):
    # q'_j = sum(g[j:]) (one-based j), low_j = 1 - q'_j - low_{j-1}.
    return (g[0], g[1], g[0]+g[2], g[1]+g[3], g[0]+g[2]+g[4])


def inverse_density(g, unit):
    low = lows_from_gaps(g)
    assert 0 <= low[0] <= low[1] <= low[2] <= low[3] <= low[4] and 2*low[4] <= unit
    total = 0
    for mask in MASKS:
        high = [unit-v for v, chosen in zip(low, mask) if chosen]
        rest = [v for v, chosen in zip(low, mask) if not chosen]
        q = [unit, *high, *reversed(rest), 0]
        total += prod(a-b for a,b in zip(q, q[1:]))
    return total


def mass(poly, rules):
    direct = Fraction(0)
    after_one = Fraction(0)
    for simplex in poly.simplices():
        vertices = [poly.vertices[i] for i in simplex]
        determinant = abs(det(vertices))
        scale = lcm(*(x.denominator for v in vertices for x in v))
        columns = tuple(zip(*(tuple(int(scale*x) for x in v) for v in vertices)))
        for size, weight, nodes in rules:
            unit = scale*size
            direct_integer_sum = after_integer_sum = 0
            for node in nodes:
                g = tuple(sum(x*b for x,b in zip(column,node)) for column in columns)
                direct_integer_sum += prod(g)
                after_integer_sum += inverse_density(g, unit)
            coeff = factorial(11) * determinant * weight / unit**6
            direct += coeff * direct_integer_sum
            after_one += coeff * after_integer_sum
    return direct, after_one


def report(label, x):
    print(f"{label}: {x} = {float(x):.12g}")


def main():
    rules = quadrature()
    faces, inter = within_interlevel()
    pi = Poly.root().restrict([*sorted(inter), *sorted(faces)])
    # Linear inequalities at vertices ensure the same 32 affine inverse
    # branches and nonnegative gaps throughout the closed polytope.
    for vertex in pi.vertices:
        inverse_density(vertex, 1)
    direct, by_one = mass(pi, rules)
    previous = json.loads((GRAPH / "pi-geometry.json").read_text())
    assert direct == Fraction(previous["pi_original_direct_weight"])
    report("original direct chamber mass", direct)
    report("original chamber mass by one step", by_one)
    report("new chamber mass at one step", by_one - direct)
    parts = {}
    for row, (phase, level, _) in sorted(faces.items()):
        winner = pi.restrict(minus(other,row) for other in sorted(faces) if other != row)
        if winner is not None:
            d, b = mass(winner, rules)
            parts[level] = parts.get(level, 0) + b
            old = next(r for r in previous["rows"] if tuple(r["gap_row"]) == row)
            assert d == Fraction(old["original_direct_pi_mass"])
    assert set(parts) == {2, 3} and sum(parts.values()) == by_one
    for level in sorted(parts):
        report(f"by-one source, first winner level {level}", parts[level])
    report("even-only level-3 period-five second return, original weight", parts[3]/3)


if __name__ == "__main__":
    main()
