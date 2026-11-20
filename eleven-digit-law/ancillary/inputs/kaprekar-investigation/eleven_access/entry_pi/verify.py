"""Exact finite checks for the analytic entry and stable-set arguments."""
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path

from regular_geometry import active_at_boundary, dot


def regular_rows():
    rows = json.loads(Path("catalogue.json").read_text())
    assert len(rows) == 191
    assert sum(row["paired"] for row in rows) == 30
    boundary_counts = Counter()
    for row in sorted(rows, key=lambda item: (-item["dimension"], item["word"])):
        word = row["word"]
        if not row["paired"] or not 0 < row["dimension"] < 5:
            continue
        if len(word) != 5 or any("-" in phase for phase in word):
            continue
        active = active_at_boundary(word)
        for i, (phase, comparison, form, direction) in enumerate(active):
            values = [dot(form, v) for _, _, _, v in active]
            if values[i] == 2 and min(values) >= 0:
                assert len(set(form)) > 1
                boundary_counts[row["dimension"]] += 1
                print(f"dimension {row['dimension']} {' '.join(word)}"
                      f" phase {phase + 1} d{comparison + 1}: {values}")
                break
        else:
            raise AssertionError(word)
    assert boundary_counts == {4: 2, 3: 6, 2: 7, 1: 3}


def gap_vector(q):
    return tuple(a - b for a, b in zip((0, *q), (*q, 0)))


def pair_step(gaps, positive):
    """Linearization obtained from the explicitly ordered digit offsets."""
    q = tuple(sum(gaps[i:]) for i in range(1, 6))
    digits = [0]
    for i, sign in enumerate((*positive, positive[1], positive[0])):
        pair = q[i], -q[4 - i]
        digits.extend(pair if sign else pair[::-1])
    return gap_vector(tuple(digits[i] - digits[10 - i] for i in range(5)))


def chart(alpha, eta, x, y, z):
    return (alpha, alpha + x, alpha + y, eta,
            -4 * alpha - eta - x - y - z, alpha + z)


def two_return():
    for first in (False, True):
        for second in (False, True):
            for j in range(5):
                alpha, eta, x, y, z = (int(i == j) for i in range(5))
                s = 0 if first else -z
                t = 0 if second else -y - s
                a = chart(alpha, eta, x, y, z)
                b = pair_step(a, (first, False, False))
                assert b[5] - b[0] == y + s
                result = pair_step(b, (second, True, True))
                expected = chart(alpha - s - t, eta - t,
                                 y + s + 2*t, z + 2*s + 2*t, x - z + t)
                assert result == expected
    print("Two return: all four exact branch matrices agree")


def constants():
    threshold = 1000
    theta = Q(199, 200)
    assert Q(101, 99)**2 < Q(25, 24)
    assert 1/(Q(99, 100)*(threshold - 1)) < Q(1, 100)
    assert Q(25, 24)**2 * Q(2, 3)**2 * 2/(1 - Q(4, threshold)) < theta**2
    assert Q(25, 24)**4 * Q(4, 3)**4 * Q(101, 99) * (1 + Q(2, threshold))**4/4 < theta**4
    assert (2 + 32*theta) / (1 - theta) / 10**5 <= Q(17, 250)
    assert Q(1, 10**7) / (1 - theta) / Q(1, 1000) <= Q(1, 50)
    assert 1 - Q(17, 250) - Q(1, 50) == Q(114, 125)
    scale = Q(1, 10**40)
    assert 8*10**10*scale < Q(1, 10**28)
    assert 896*Q(1, 10**12) < Q(1, 10**9)
    assert 20*Q(1, 10**9) < Q(1, 10**5)
    assert Q(1, 6) - 4*Q(11, 10**4) - Q(21, 10**4) - Q(1, 10**5) > Q(1, 7)
    assert Q(1, 10**4)**2 * scale**3 / (8*10**20) == Q(1, 8*10**148)
    print("Fractional-moment constants, neighborhood and volume bounds agree")


if __name__ == "__main__":
    regular_rows()
    two_return()
    constants()
