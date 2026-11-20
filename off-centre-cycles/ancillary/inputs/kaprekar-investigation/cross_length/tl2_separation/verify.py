#!/usr/bin/env python3
"""Replay the two collar censuses, one original sample, and the affine cohort."""

import argparse
from collections import Counter
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ZERO = (0,) * 6


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def subtract(x, y):
    return tuple(a - b for a, b in zip(x, y))


def gaps(x):
    return subtract((0, *x), (*x, 0))


def slots(x, mask, unit):
    result = [-unit]
    for i in range(6):
        left = x[i] - unit * (i == 5)
        right = -x[5-i] - unit * (i < 5)
        result.extend((right, left) if mask & (1 << i) else (left, right))
    return result


def comparisons(x, mask, unit):
    values = slots(x, mask, unit)
    # Assign a numerical tie to the left member, using the seven-unit lattice.
    return tuple(values[2*i+1] - values[2*i+2] - unit * bool(mask & (1 << i))
                 for i in range(6))


def forced(x, mask, unit):
    values = slots(x, mask, unit)
    return subtract(values[:6], values[:6:-1])


def phi(x):
    left = slots(x, 0, 7)
    mask = sum(1 << i for i in range(6) if left[2*i+1] < left[2*i+2])
    return forced(x, mask, 7)


def actual(x, base):
    q = tuple((u+(6-i)*base)//7 for i, u in enumerate(x))
    assert all(7*v-(6-i)*base == u for i, (u, v) in enumerate(zip(x, q)))
    assert base > q[0] >= q[1] >= q[2] >= q[3] >= q[4] >= q[5] >= 1
    value = sum(v * (base**(12-i)-base**i) for i, v in enumerate(q))
    digits = []
    for _ in range(13):
        value, digit = divmod(value, base)
        digits.append(digit)
    assert value == 0
    digits.sort(reverse=True)
    return tuple(7*(digits[i]-digits[12-i])-(6-i)*base for i in range(6))


def rank(vectors):
    matrix = [list(map(Fraction, row)) for row in vectors]
    pivots = 0
    for column in range(6):
        pivot = next((j for j in range(pivots, len(matrix)) if matrix[j][column]), None)
        if pivot is None:
            continue
        matrix[pivot], matrix[pivots] = matrix[pivots], matrix[pivot]
        for j in range(pivots+1, len(matrix)):
            ratio = matrix[j][column] / matrix[pivots][column]
            matrix[j] = [a-ratio*b for a, b in zip(matrix[j], matrix[pivots])]
        pivots += 1
    return pivots


def check_route():
    data = json.loads((HERE / "class0-route.json").read_text())
    seed = tuple(data["x"])
    directions = tuple(tuple(record["growth"]) for record in data["words"])
    assert len(directions) == rank(directions) == 5
    assert all(value % 7 == 0 for row in (seed, *directions) for value in row)
    max_constant, max_linear = 0, 0
    for i, record in enumerate(data["words"]):
        original_slopes = directions[:i+1]
        slopes = original_slopes
        offset = seed
        for mask in record["masks"]:
            assert 0 <= mask < 64
            assert all(value >= 0 for v in slopes for value in comparisons(v, mask, 0))
            # The current word starts with at least one remaining copy.
            assert min(comparisons(add(offset, slopes[-1]), mask, 7)) >= 0
            max_constant = max(max_constant, *map(abs, gaps(offset)))
            max_linear = max(max_linear, *(sum(map(abs, coefficients))
                                           for coefficients in zip(*map(gaps, slopes))))
            offset = forced(offset, mask, 7)
            slopes = tuple(forced(v, mask, 0) for v in slopes)
        assert slopes == original_slopes
        assert offset == subtract(seed, original_slopes[-1])
    seen = {}
    state = seed
    while state not in seen:
        assert len(seen) < 500
        seen[state] = len(seen)
        state = phi(state)
    preperiod = seen[state]
    assert len(seen) - preperiod == data["period"] == 177
    assert next(j for x, j in seen.items() if -min(gaps(x)) <= 28) < preperiod
    assert max(-min(gaps(x)) for x in seen) == data["peak"]
    assert max(-min(gaps(x)) for x, j in seen.items() if j >= preperiod) == 28
    for state in seen:
        max_constant = max(max_constant, *map(abs, gaps(state)))
    denominator = 32 * max_linear
    cutoff = max(32 * max_constant, 16)
    base = cutoff + (-cutoff) % 7
    assert all(actual(x, base) == phi(x) for x in seen)
    # Explicit trajectory control through each affine word and fixed suffix.
    check_base = denominator * 21 + base
    state = seed
    for v in directions:
        for _ in range(2):
            state = add(state, v)
    for i in reversed(range(len(directions))):
        for _ in range(2):
            for mask in data["words"][i]["masks"]:
                successor = forced(state, mask, 7)
                assert successor == actual(state, check_base)
                state = successor
    assert state == seed
    constant = Fraction(factorial(13), 8 * 16**6 * denominator**5)
    print("rank", rank(directions), "word lengths", *(len(w["masks"]) for w in data["words"]))
    print("fixed suffix preperiod", preperiod, "least period", data["period"])
    print("route bounds C, S, D, B0", max_constant, max_linear, denominator, cutoff)
    print("original probability lower bound", constant, "/ B")


def check_premises(root):
    sources = (
        root / "kaprekar-investigation/cross_length/tl_labels/arrivals.cpp",
        root / "docs/drafts/research-direction/check-eighth-compass/kap.cpp",
    )
    expected = {
        27: Counter({4: 39552, 8: 392, 16: 108100, 33: 959524}),
        28: Counter({1: 9994, 4: 1734, 8: 1, 177: 1333175}),
        34: Counter({4: 134430, 8: 1454, 16: 383751, 33: 3318745}),
    }
    with tempfile.TemporaryDirectory(prefix="tl154-premise-") as directory:
        work = Path(directory)
        for source in sources:
            shutil.copyfile(source, work / source.name)
            subprocess.run(["clang++", "-O3", "-std=c++17", source.name, "-o", source.stem],
                           cwd=work, check=True)
        for radius, cutoff in ((27, 83), (28, 84), (34, 104)):
            out = subprocess.check_output(
                ["./arrivals", "census", str(radius), str(cutoff)], cwd=work, text=True)
            counts = Counter()
            for line in out.splitlines():
                r, guard, period, count, _ = map(int, line.split())
                assert (r, guard) == (radius, cutoff)
                assert period not in counts
                counts[period] = count
            assert counts == expected[radius]
            assert sum(counts.values()) == comb(radius+6, 6)
            print("collar", radius, sorted(counts.items()), "total", sum(counts.values()))
        out = subprocess.check_output(
            ["./kap", "pop", "13", "110005", "1", "1500", "151001"], cwd=work, text=True)
        sample = Counter(int(row.split()[2]) for row in out.splitlines())
        assert sample == {1: 11, 177: 1489}
        print("class-0 original sample at B=110005", sorted(sample.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=HERE.parents[2])
    parser.add_argument("--route-only", action="store_true")
    args = parser.parse_args()
    if not args.route_only:
        check_premises(args.root.resolve())
    check_route()
