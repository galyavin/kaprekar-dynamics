"""Verify the affine returns throughout the archived sorting cells."""

from collections import Counter
from fractions import Fraction as F
import gzip
import json
from pathlib import Path


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def reduced(v, basis):
    for i, row in sorted(basis.items()):
        d = v[i]
        if d:
            v = tuple(a - d * b for a, b in zip(v, row))
    return v


def add_equation(v, basis):
    v = reduced(v, basis)
    k = next((i for i, x in enumerate(v) if x), None)
    if k is not None:
        basis[k] = tuple(F(x, v[k]) for x in v)


def period(cell):
    # Homogeneous functions of the six initial gaps. Their sum represents 1.
    q = tuple(map(F, cell["q_witness"]))
    gap = sub((1, *q), (*q, 0))
    start = tuple(tuple(int(i >= j) for i in range(6)) for j in range(1, 6))
    one = (1,) * 6
    history = [start]
    basis = {}
    for face in cell["faces"]:
        source = history[-1]
        values = tuple(dot(row, gap) for row in source)
        groups = {}
        for x in values:
            counts = groups.setdefault(abs(2 * x - 1), [0, 0])
            counts[int(2 * x <= 1)] += 1
        zero = sum(groups.pop(0, []))
        assert face == [zero, [groups[i] for i in sorted(groups, reverse=True)], 1 in groups]
        rows = (one, *source, *(sub(one, row) for row in source))
        ordered = sorted(rows, key=lambda x: dot(x, gap), reverse=True)
        for a, b in zip(ordered, ordered[1:]):
            v = sub(a, b)
            if dot(v, gap) == 0:
                add_equation(v, basis)
        history.append(tuple(sub(ordered[i], ordered[-1 - i]) for i in range(5)))
    assert all(dot(sub(x, y), gap) == 0 for x, y in zip(history[-1], start))
    for x, y in zip(history[-1], start):
        add_equation(sub(x, y), basis)
    # Equal magnitude groups, including the endpoint, fix all ties in raw
    # sorted digits on this guard cell. Reduce each affine difference exactly.
    h = next(j for j, h in enumerate(history[1:], 1)
             if all(not any(reduced(sub(x, y), basis)) for x, y in zip(h, start)))
    return h, 5 - len(basis)


def main():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    archive = Path(__file__).resolve().parents[1] / "original" / "energy.json.gz"
    with gzip.open(archive, "rt", encoding="utf-8") as stream:
        cells = json.load(stream)
    assert len(cells) == 268, len(cells)
    periods = Counter()
    for cell in cells:
        assert cell["return_length"] == len(cell["faces"]), cell
        h, _ = period(cell)
        assert cell["return_length"] % h == 0 and h % cell["witness_period"] == 0, cell
        periods[h] += 1
    assert periods == {1: 16, 2: 14, 3: 80, 5: 153, 6: 2, 7: 2, 12: 1}, periods
    print(f"Exact affine return reduction verified all {len(cells)} archived cells.")
    print("Physical return counts", dict(sorted(periods.items())))


if __name__ == "__main__":
    main()
