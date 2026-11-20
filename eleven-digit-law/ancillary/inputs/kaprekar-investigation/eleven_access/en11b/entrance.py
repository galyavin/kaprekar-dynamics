"""Exact pullbacks and original population weights for fixed leading prefixes."""
from collections import Counter
from fractions import Fraction as Q
import json
from math import factorial, lcm, prod
from pathlib import Path
import argparse
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ACCESS = HERE.parent
sys.path[:0] = [
    str(ACCESS / "stage_graph"),
    str(ROOT / "docs/drafts/research-direction/astra-strategy-12"),
]
from one_step import MASKS, lows_from_gaps, quadrature
from pi_geometry import det, dot, minus, within_interlevel
from clipping import Poly, integer_vertices
from geometry import basis, integer_row

IDENTITY = tuple(basis(i) for i in range(6))
ZERO = (0,) * 6
ONE = (1,) * 6


def rows(function):
    return tuple(zip(*(function(vertex) for vertex in IDENTITY)))


LOWS = rows(lows_from_gaps)
IMAGE_CUTS = (LOWS[0], *(minus(b, a) for a, b in zip(LOWS, LOWS[1:])),
              minus(ONE, tuple(2 * v for v in LOWS[-1])))


def inverse_rows(mask):
    highs = [minus(ONE, low) for low, chosen in zip(LOWS, mask) if chosen]
    lows = [low for low, chosen in zip(LOWS, mask) if not chosen]
    points = (ONE, *highs, *reversed(lows), ZERO)
    return tuple(minus(a, b) for a, b in zip(points, points[1:]))


INVERSES = tuple(inverse_rows(mask) for mask in MASKS)
PI_PREDECESSOR = tuple(IDENTITY[j] for j in (0, 2, 4, 5, 3, 1))
assert INVERSES[21] == PI_PREDECESSOR
for branch in INVERSES:
    assert abs(det(tuple(tuple(Q(v) for v in row) for row in branch))) == 1
# Complementing a mask reverses the original gaps. Their density is symmetric.
for i, inverse in enumerate(INVERSES):
    assert inverse[::-1] == INVERSES[31 - i]


def pull(row, chart):
    return tuple(dot(row, column) for column in zip(*chart))


def compose(left, right):
    return tuple(pull(row, right) for row in left)


def section():
    within, between = within_interlevel()
    return Poly.root().restrict((*sorted(between), *sorted(within)))


def image_domain(poly, chart):
    rows_to_clip = []
    vertices = integer_vertices(poly.vertices)
    for row in IMAGE_CUTS:
        cut = integer_row(pull(row, chart))
        assert any(cut)  # Invertible affine branches do not collapse image walls.
        values = tuple(dot(cut, v) for v in vertices)
        if max(values) <= 0:
            return None
        if min(values) < 0:
            rows_to_clip.append((sum(v > 0 for v in values), cut))
    return poly.restrict(cut for _, cut in sorted(set(rows_to_clip)))


def determinant_integer(columns):
    # Bareiss elimination keeps the determinant calculation in exact integers.
    matrix = [list(row) for row in columns]
    sign, previous = 1, 1
    for i in range(5):
        pivot = next((j for j in range(i, 6) if matrix[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
            sign = -sign
        value = matrix[i][i]
        for j in range(i + 1, 6):
            for k in range(i + 1, 6):
                numerator = value * matrix[j][k] - matrix[j][i] * matrix[i][k]
                matrix[j][k], remainder = divmod(numerator, previous)
                assert not remainder
        previous = value
    return sign * matrix[-1][-1]


def python_sums(columns, scale, direct, rules):
    for size, _, nodes in rules:
        unit = scale * size
        acc = 0
        for node in nodes:
            g = tuple(dot(c, node) for c in columns)
            if direct:
                acc += prod(g)
            else:
                low = lows_from_gaps(g)
                for mask in MASKS[:16]:
                    high = [unit-v for v, chosen in zip(low, mask) if chosen]
                    rest = [v for v, chosen in zip(low, mask) if not chosen]
                    q = [unit, *high, *reversed(rest), 0]
                    acc += 2 * prod(a-b for a, b in zip(q, q[1:]))
        yield acc


def measure(poly, rules, direct=False, evaluator=None):
    """Integrate density of one inverse on the current output polytope."""
    total = Q(0)
    for simplex in poly.simplices():
        vertices = [poly.vertices[i] for i in simplex]
        scale = lcm(*(x.denominator for v in vertices for x in v))
        columns = tuple(zip(*(tuple(int(scale*x) for x in v) for v in vertices)))
        jacobian = Q(abs(determinant_integer(columns)), scale**6)
        assert jacobian > 0
        if evaluator is not None and scale <= 300000:
            assert evaluator.stdin is not None and evaluator.stdout is not None
            fields = (int(direct), scale, *(c for row in columns for c in row))
            evaluator.stdin.write(" ".join(map(str, fields)) + "\n")
            evaluator.stdin.flush()
            values = tuple(map(int, evaluator.stdout.readline().split()))
            assert len(values) == len(rules)
        else:
            values = python_sums(columns, scale, direct, rules)
        for (size, weight, _), value in zip(rules, values):
            total += factorial(11) * jacobian * weight * Q(value, (scale * size)**6)
    return total


def mapped(poly, chart):
    return Poly(tuple(tuple(dot(r, v) for r in chart) for v in poly.vertices), poly.facets)


def pullbacks(output, k, newly_entering=False, later_inverse_masks=None):
    # Keep geometrically possible partial inverse paths. For the final step,
    # the whole current domain has all 32 inverses, integrated as one density.
    frontier = [(output, IDENTITY)]
    for depth in range(1, k):
        accepted = []
        branches = INVERSES
        if newly_entering and depth == 1:
            branches = tuple(b for b in INVERSES if b != PI_PREDECESSOR)
        if depth >= 2 and later_inverse_masks is not None:
            branches = tuple(INVERSES[i] for i in later_inverse_masks)
        for poly, chart in frontier:
            poly = image_domain(poly, chart)
            if poly is not None:
                accepted.extend((poly, compose(branch, chart)) for branch in branches)
        frontier = accepted
    for poly, chart in frontier:
        poly = image_domain(poly, chart)
        if poly is not None:
            yield poly, chart


def integrate(output, k, rules, evaluator=None, newly_entering=False, later_inverse_masks=None):
    if k == 0:
        return measure(output, rules, direct=True, evaluator=evaluator), 1, 1
    assert not newly_entering or k >= 2
    domains = pullbacks(output, k, newly_entering, later_inverse_masks)
    history_count = 0
    # Distinct histories may map to the same final integration domain.
    copies = Counter()
    representatives = {}
    for poly, chart in domains:
        history_count += 1
        image = mapped(poly, chart)
        key = tuple(sorted(image.vertices))
        copies[key] += 1
        representatives.setdefault(key, image)
    value = sum(copies[key] * measure(poly, rules, evaluator=evaluator)
                for key, poly in representatives.items())
    return value, history_count * 32, len(representatives)


def ledger(k_max, evaluator=None, checkpoints=None):
    assert 0 <= k_max <= 4
    pi = section()
    rules = quadrature()
    previous = Q(0)
    entries = []
    for k in range(k_max + 1):
        total, branches, distinct = integrate(pi, k, rules, evaluator)
        assert total >= previous
        entry = {"step": k, "first_even": total - previous, "first_odd": total - previous,
                 "by_step_each_class": total, "outside_by_step_each_class": 1 - total,
                 "full_dimensional_inverse_paths": branches,
                 "distinct_final_integration_domains": distinct}
        print(f"{k}: {float(total):.12g}; {branches} paths; {distinct} domains", file=sys.stderr, flush=True)
        entries.append(entry)
        if checkpoints is not None:
            checkpoints.write_text(json.dumps(entries, default=str, indent=2) + "\n")
        previous = total
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--evaluator", type=Path)
    parser.add_argument("--checkpoints", type=Path)
    args = parser.parse_args()
    if args.evaluator:
        with subprocess.Popen([str(args.evaluator.resolve())], text=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE) as native:
            entries = ledger(args.depth, native, args.checkpoints)
            native.stdin.close()
            assert native.wait() == 0
    else:
        entries = ledger(args.depth, checkpoints=args.checkpoints)
    print(json.dumps(entries, default=str, indent=2))


if __name__ == "__main__":
    main()
