"""Replay the guarded order-8 families and compute their rational cube weights."""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from math import factorial
from pathlib import Path
import subprocess
import sys

ORDER = 8
REFERENCES = (
    {
        "name": "period1712",
        "source": (6353, 5202, 1589, 4459, 8536, 12967, 9032, 11185),
        "intervals": ((0, 8), (1414, 1446), (7214, 7478), (9454, 9670),
                      (10294, 11294), (17091, 17139), (21493, 24197), (29639, 31351)),
        "expected": (1712, 369, 194),
    },
    {
        "name": "period368",
        "source": (40763, 7002, 61702, 54124, 81085, 1790, 53786, 36842),
        "intervals": ((0, 8), (1602, 1618), (9422, 9574), (71560, 71688),
                      (93079, 93391), (128476, 130108), (142481, 144665), (171781, 172149)),
        "expected": (368, 66, 39),
    },
)


def step(x, forced=False):
    z = x[-1]
    sign = 1 if forced else (z > 0) - (z < 0)
    loss = 0 if forced else max(0, -1 - z)
    return (z - 3 * sign + 2 * loss, x[0] + 2 * loss,
            *(v - sign + loss for v in x[1:-1]))


def replay(reference):
    source = reference["source"]
    intervals = reference["intervals"]
    assert len(intervals) == ORDER
    assert all(0 <= start < end for start, end in intervals)
    assert all(a[1] <= b[0] for a, b in zip(intervals, intervals[1:]))
    wanted = {t for pair in intervals for t in pair}
    first_bad = {}
    snapshots = {}
    x = source
    for t in range(intervals[-1][1] + 1):
        if t in wanted:
            snapshots[t] = x
        if x[-1] <= 0:
            first_bad.setdefault((ORDER - 1 - t) % ORDER, t)
        x = step(x)
    assert len(first_bad) == ORDER - 1 and min(first_bad.values()) >= ORDER
    retired = tuple(first_bad)
    last = next(label for label in range(ORDER) if label not in first_bad)
    rays = []
    for index, (start, end) in enumerate(intervals):
        assert (end - start) % ORDER == 0
        target = (*retired, last)[index]
        if index < len(retired):
            assert end <= first_bad[target]
        difference = tuple(u - v for u, v in zip(snapshots[start], snapshots[end]))
        shift = start % ORDER
        ray = difference[shift:] + difference[:shift]
        assert min(ray) >= 0 and ray[target] > 0
        assert not any(ray[j] for j in retired[:index])
        assert all(not v or first_bad.get(j, end) >= end for j, v in enumerate(ray))
        rays.append(ray)
    assert all(j == last or v == 0 for j, v in enumerate(rays[-1]))
    slot = (last + intervals[-1][0]) % ORDER
    section = snapshots[intervals[-1][0]]
    return (*retired, last), rays, section, slot


def canonical_section(section, slot, period):
    state = list(section)
    state[slot] = 0
    state = tuple(state)
    visited = []
    for _ in range(period):
        if slot == ORDER - 1:
            visited.append(state)
        state = list(step(state, slot == ORDER - 1))
        slot = (slot + 1) % ORDER
        state[slot] = 0
        state = tuple(state)
    assert len(set(visited)) == len(visited) == period // ORDER
    return min(visited)


def phase_table(binary, period, descent, good, section):
    args = [binary, "component", str(ORDER)]
    request = " ".join(map(str, (ORDER - 1, 0, period, *section))) + "\n"
    rows = subprocess.run(args, input=request, text=True, capture_output=True,
                          check=True).stdout.splitlines()
    header = rows[0].split()
    assert (int(header[1]), int(header[3])) == (period, descent)
    assert rows[-1] == f"good {good} total {descent}"
    mark = dict(tuple(map(int, line.split()[:2])) for line in rows[1:-1])
    assert sorted(mark) == list(range(descent)) and sum(mark.values()) == good
    return [residue for residue, positive in sorted(mark.items()) if positive]


def eliminate(matrix, rhs=None):
    """Exact determinant, or exact solution when a right-hand side is supplied."""
    a = [list(map(Fraction, row)) for row in matrix]
    if rhs is not None:
        for row, value in zip(a, rhs):
            row.append(Fraction(value))
    n = len(a)
    determinant = Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return None if rhs is not None else Fraction(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            determinant = -determinant
        determinant *= a[j][j]
        a[j] = [value / a[j][j] for value in a[j]]
        for i in range(j + 1, n):
            multiplier = a[i][j]
            if multiplier:
                a[i] = [value - multiplier * leading for value, leading in zip(a[i], a[j])]
    if rhs is None:
        return determinant
    values = [Fraction(0)] * n
    for j in range(n - 1, -1, -1):
        values[j] = a[j][-1] - sum(a[j][i] * values[i] for i in range(j + 1, n))
    return tuple(values)


def dot(a, b):
    return sum(u * v for u, v in zip(a, b))


def coefficient_volume(rays):
    """Volume of {a>=0: Ua<=1} by exact vertices and pulling triangulation."""
    rows = list(dict.fromkeys(zip(*rays)))
    rows = [row for row in rows if not any(
        row != other and all(a <= b for a, b in zip(row, other)) for other in rows
    )]
    inequalities = [(row, 1) for row in rows]
    inequalities.extend(
        (tuple(-int(i == j) for i in range(ORDER)), 0) for j in range(ORDER)
    )
    vertices = set()
    for constraints in combinations(inequalities, ORDER):
        matrix, rhs = zip(*constraints)
        vertex = eliminate(matrix, rhs)
        if vertex is not None and all(dot(row, vertex) <= b for row, b in inequalities):
            vertices.add(vertex)
    vertices = sorted(vertices)
    walls = [frozenset(i for i, vertex in enumerate(vertices) if dot(row, vertex) == b)
             for row, b in inequalities]

    @lru_cache(maxsize=None)
    def triangulate(face):
        if len(face) == 1:
            return (tuple(face),)
        origin = min(face)
        intersections = {face & wall for wall in walls}
        proper = [part for part in intersections if part and part != face]
        return tuple((origin, *simplex)
                     for part in proper if origin not in part
                     and not any(part < other for other in proper)
                     for simplex in triangulate(part))

    simplices = triangulate(frozenset(range(len(vertices))))
    volume = Fraction(0)
    for simplex in simplices:
        assert len(simplex) == ORDER + 1
        anchor, *rest = (vertices[i] for i in simplex)
        matrix = [tuple(x - y for x, y in zip(vertex, anchor)) for vertex in rest]
        determinant = eliminate(matrix)
        assert determinant
        volume += abs(determinant) / factorial(ORDER)
    return volume, len(rows), len(vertices), len(simplices)


def main(binary):
    components = []
    covered = lower = loss_lower = Fraction(0)
    retirement_words = set()
    for ref in REFERENCES:
        retirement, rays, section, slot = replay(ref)
        period, descent, good = ref["expected"]
        assert descent == max(rays[-1])
        canonical = canonical_section(section, slot, period)
        residues = phase_table(binary, period, descent, good, canonical)
        volume, rows, vertices, simplices = coefficient_volume(rays)
        for rotation in range(ORDER):
            word = tuple((label + rotation) % ORDER for label in retirement)
            assert word not in retirement_words
            retirement_words.add(word)
        weight = ORDER * descent * volume
        marked = ORDER * good * volume
        unmarked = ORDER * (descent - good) * volume
        covered += weight
        lower += marked
        loss_lower += unmarked
        components.append({
            "name": ref["name"], "source": ref["source"],
            "retirement_order": retirement, "rays": rays,
            "determinant": str(abs(eliminate(rays))), "period": period, "descent": descent,
            "canonical_section": canonical, "good_residues": residues,
            "nondominated_cube_rows": rows, "vertices": vertices, "simplices": simplices,
            "unrotated_weight_per_residue": str(volume),
            "eight_rotations_cube_weight": str(weight),
            "eight_rotations_cube_weight_approx": f"{float(weight):.9g}",
            "cube_weight_of_zero_loss": str(marked),
            "cube_weight_of_positive_loss": str(unmarked),
        })
    assert covered == lower + loss_lower
    print(json.dumps({
        "components": components,
        "covered": str(covered), "covered_approx": f"{float(covered):.9g}",
        "uncovered": str(1 - covered), "zero_loss_lim_inf_lower": str(lower),
        "one_minus_zero_loss_lim_sup_lower": str(loss_lower),
        "zero_loss_lim_sup_upper": str(1 - loss_lower),
        "zero_loss_lower_approx": f"{float(lower):.9g}",
        "positive_loss_lower_approx": f"{float(loss_lower):.9g}",
    }, indent=2))


if __name__ == "__main__":
    main(str(Path(sys.argv[1]).resolve()))
