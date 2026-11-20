#!/usr/bin/env python3
"""Check the full period-177 certificate with exact integer and rational arithmetic."""

from fractions import Fraction
from itertools import product
from math import factorial
from pathlib import Path
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def plus(x, y, scale=1):
    return tuple(a + scale*b for a, b in zip(x, y))


def gaps(x):
    return plus((0, *x), (*x, 0), -1)


def digits(x, unit):
    return [(-unit, -unit)] + [
        (x[i] - unit*(i == 5), -x[5-i] - unit*(i != 5))
        for i in range(6)
    ]


def branch(x, mask, unit):
    assert 0 <= mask < 64
    pairs = digits(x, unit)
    slots = [-unit]
    comparisons = []
    for i, (a, b) in enumerate(pairs[1:]):
        swapped = (mask >> i) & 1
        high, low = (b, a) if swapped else (a, b)
        comparisons.append(high - low - unit*swapped)
        slots.extend((high, low))
    successor = tuple(slots[i] - slots[12-i] for i in range(6))
    return successor, tuple(comparisons)


def paired(x):
    mask = sum(1 << i for i, (a, b) in enumerate(digits(x, 1)[1:]) if a < b)
    return branch(x, mask, 1)[0]


def physical(x, base):
    """Recover the symmetric differences by exact whole-integer subtraction."""
    assert base % 7 == 0
    q = tuple(x[i] + (6-i)*(base//7) for i in range(6))
    assert base > q[0] > q[1] > q[2] > q[3] > q[4] > q[5] > 0
    word = sum(v * (base**(12-i) - base**i) for i, v in enumerate(q))
    out = []
    for _ in range(13):
        word, digit = divmod(word, base)
        out.append(digit)
    assert word == 0
    out.sort(reverse=True)
    return tuple(out[i] - out[12-i] - (6-i)*(base//7) for i in range(6))


def determinant(columns):
    matrix = [list(map(Fraction, row)) for row in zip(*columns)]
    det = Fraction(1)
    for i in range(6):
        pivot = next((j for j in range(i, 6) if matrix[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            matrix[pivot], matrix[i] = matrix[i], matrix[pivot]
            det = -det
        det *= matrix[i][i]
        for j in range(i+1, 6):
            ratio = matrix[j][i] / matrix[i][i]
            matrix[j] = [a-ratio*b for a, b in zip(matrix[j], matrix[i])]
    assert det.denominator == 1
    return det.numerator


class Bounds:
    def __init__(self):
        self.constant = 0
        self.slope = 0

    def include(self, center, rays):
        self.constant = max(self.constant, *(-7*g for g in gaps(center)))
        self.slope = max(self.slope, *(
            7*sum(max(-g, 0) for g in coeffs)
            for coeffs in zip(*(gaps(v) for v in rays))
        ))


def check_route(record):
    x = tuple(record["start"])
    rays = tuple(map(tuple, record["source_rays"]))
    count_rays = tuple(map(tuple, record["count_rays"]))
    assert len(rays) == len(count_rays) == 6
    det = determinant(rays)
    assert det != 0
    repeats = iter(zip(*count_rays))
    bounds = Bounds()
    margin = None
    forward_steps = 0

    def check_values(constant, slopes):
        nonlocal margin
        margin = min(constant) if margin is None else min(margin, *constant)
        assert min(constant) >= 0
        assert all(min(c) >= 0 for c in slopes)

    for stage in record["stages"]:
        if "copies" not in stage:
            forward_steps += len(stage["masks"])
            for mask in stage["masks"]:
                bounds.include(x, rays)
                x, guard = branch(x, mask, 1)
                transported = tuple(branch(v, mask, 0) for v in rays)
                check_values(guard, tuple(g for _, g in transported))
                rays = tuple(v for v, _ in transported)
            bounds.include(x, rays)
            continue

        repeats_for_ray = next(repeats)
        base_copies = stage["copies"]
        assert base_copies >= 1 and min(repeats_for_ray) >= 0
        forward_steps += base_copies*len(stage["masks"])
        growth = tuple(stage["growth"])
        # Low and high extremal copies have respectively one and all copies left.
        end = plus(x, growth, -base_copies)
        end_rays = tuple(plus(v, growth, -h) for v, h in zip(rays, repeats_for_ray))
        left, left_rays = plus(end, growth), end_rays
        right, right_rays = x, rays
        for mask in stage["masks"]:
            bounds.include(left, left_rays)
            bounds.include(right, right_rays)
            left, left_guard = branch(left, mask, 1)
            right, right_guard = branch(right, mask, 1)
            left_moved = tuple(branch(v, mask, 0) for v in left_rays)
            right_moved = tuple(branch(v, mask, 0) for v in right_rays)
            check_values(left_guard, tuple(g for _, g in left_moved))
            check_values(right_guard, tuple(g for _, g in right_moved))
            left_rays = tuple(v for v, _ in left_moved)
            right_rays = tuple(v for v, _ in right_moved)
        assert (left, left_rays) == (end, end_rays)
        assert (right, right_rays) == (plus(x, growth, -1), rays)
        x, rays = end, end_rays
        bounds.include(x, rays)

    assert next(repeats, None) is None
    assert x == tuple(record["tail"])
    assert all(v == (0,)*6 for v in rays)
    seen = {}
    while x not in seen:
        assert len(seen) < 1000
        seen[x] = len(seen)
        bounds.include(x, rays)
        x = paired(x)
    assert seen[x] == 0
    assert len(seen) == 177

    # R <= (B/8)-1 by these exact bounds and the chosen box.
    D = 16*bounds.slope
    B0 = max(72, 16*(bounds.constant+1))
    eta = Fraction(factorial(13), 8*9**6*D**6)
    control_base = max(B0, D) + (-max(B0, D)) % 7
    assert all(physical(x, control_base) == paired(x) for x in seen)
    print("route determinant in X/7", det)
    print("fixed route clock at zero parameters", forward_steps)
    print("least period", len(seen), "peak cycle radius", max(-7*min(gaps(y)) for y in seen))
    print("minimum signed comparison at box origins", margin)
    print("C, S, D, B0", bounds.constant, bounds.slope, D, B0)
    print("uniform original probability lower bound", eta)
    return control_base, D


def check_literal(record, base, D):
    tests = [(0,)*6, (1,)*6]
    tests.extend(tuple(int(i == j) for i in range(6)) for j in range(6))
    checked = 0
    for n in tests:
        state = tuple(record["start"])
        for coefficient, ray in zip(n, record["source_rays"]):
            state = plus(state, ray, coefficient)
        changes = iter(zip(*record["count_rays"]))
        length_bound = base//D
        assert max(n) <= length_bound
        for stage in record["stages"]:
            copies = 1
            if "copies" in stage:
                copies = stage["copies"] + sum(v*c for v, c in zip(n, next(changes)))
            for _ in range(copies):
                for mask in stage["masks"]:
                    next_state, guard = branch(state, mask, 1)
                    assert min(guard) >= 0
                    assert -7*min(gaps(state))*8 + 7 <= base
                    assert next_state == physical(state, base)
                    checked += 1
                    state = next_state
        assert state == tuple(record["tail"])
    print("literal subtraction", checked, "steps from", len(tests), "inputs at base", base)


if __name__ == "__main__":
    old = HERE.parent / "tl2_separation" / "verify.py"
    subprocess.run([sys.executable, "-I", "-B", str(old), "--route-only"], check=True)
    route = json.loads((HERE / "route.json").read_text())
    base, D = check_route(route)
    check_literal(route, base, D)
