"""Independent P10 checks using the paper and released P4 affine data.

Run with python3 -B recompute.py. This module reads no certificates or scripts.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt


def folded_orbits(length):
    orbit_of = [0] * ((length + 1) // 2)
    orbit_sizes = {}
    for first in range(1, len(orbit_of)):
        if orbit_of[first]:
            continue
        slot = first
        size = 0
        while not orbit_of[slot]:
            orbit_of[slot] = first
            slot = min(2 * slot, length - 2 * slot)
            size += 1
        assert slot == first
        orbit_sizes[first] = size
    return orbit_of, orbit_sizes


def compare(gaps, width):
    return sum(gaps[-width:]) - sum(gaps[:width])


def prime_rows():
    prime = 1069
    assert all(prime % d for d in range(2, isqrt(prime) + 1))
    orbit_of, sizes = folded_orbits(prime)
    assert sizes == {1: 178, 3: 178, 7: 178}
    indicator_columns = [
        tuple(int(representative == label) for representative in orbit_of)
        for label in (3, 7)
    ]
    rows = {
        width: tuple(compare(column, width) for column in indicator_columns)
        for width in range(2, len(orbit_of) // 2 + 1)
    }
    expected = {2: (1, 0), 3: (2, 0), 4: (1, 1), 31: (-3, 0), 34: (-1, -1)}
    assert {width: rows[width] for width in expected} == expected
    assert len(rows) == 266 and (0, 0) not in rows.values()
    assert next(r for r in range(1, prime) if pow(2, r, prime) in (1, prime - 1)) == 178
    assert pow(2, 178, prime) == prime - 1
    assert pow(2, 356, prime**2) == 1 + 417 * prime

    # Enumerate the prime-square permutation. Anchor has the root weight.
    square_orbits, square_sizes = folded_orbits(prime**2)
    assert square_sizes[1] == 178 * prime
    square_orbits[0] = 1
    row534 = Counter(square_orbits[-534:])
    row534.subtract(square_orbits[:534])
    assert not any(row534.values())
    for width, (x, y) in expected.items():
        row = Counter(square_orbits[-width:])
        row.subtract(square_orbits[:width])
        assert row[3] == x and row[7] == y
        assert not any(value for label, value in row.items() if label not in (1, 3, 7))
    print(f"1069: orbit sizes {sizes}, {len(rows)} nonzero rows, witnesses {expected}")
    print(f"1069^2: root size {square_sizes[1]}, d_534 vanishes, five unit-offset rows unchanged")


def p4_scalar():
    # Published P4 sections 4.2 and 4.3, BB -> AB continuation.
    tangent = (4, 8, 7, 6)
    drift = (-12, -24, -12, -18)

    def a1_direction(q):
        u, v, w, _ = q
        return u + v - w

    decrement = -a1_direction(drift)
    value = a1_direction(tangent)
    image = {value * multiplier % decrement for multiplier in range(decrement)}
    kernel_step = next(multiplier for multiplier in range(1, decrement + 1)
                       if value * multiplier % decrement == 0)
    child = tuple(kernel_step * v + value * d for v, d in zip(tangent, drift))
    assert decrement == 24 and value == 5
    assert image == set(range(24)) and kernel_step == 24
    assert child == (36, 72, 108, 54)
    print(f"[G4] BB -> AB: guard image size {len(image)}, each exit weight {Fraction(1, len(image))}, "
          f"inherited exit direction {child}")


def from_gaps(gaps):
    return tuple(sum(gaps[i:]) for i in range(1, len(gaps)))


def to_gaps(q, base):
    return tuple(a - b for a, b in zip((base,) + q, q + (0,)))


def word_value(digits, base):
    value = 0
    for digit in digits:
        value = value * base + digit
    return value


def literal_step(gaps, base):
    q = from_gaps(gaps)
    k = len(q)
    digits = q + (0,) * (k + 1)
    difference = word_value(digits, base) - word_value(reversed(digits), base)
    unpacked = []
    for _ in digits:
        difference, digit = divmod(difference, base)
        unpacked.append(digit)
    assert difference == 0 and q[-1] > 0
    ordered = sorted(unpacked, reverse=True)
    out = tuple(a - b for a, b in zip(ordered[:k], reversed(ordered[-k:])))

    # Check physical level order against the actual digits.
    levels = tuple(
        (q[i] - (i == k - 1), base - q[k - 1 - i] - (i < k - 1))
        for i in range(k)
    )
    assert unpacked[::-1] == ([level[0] for level in levels] + [base - 1]
                              + [level[1] for level in levels])
    assert all(min(high) >= max(low) for high, low in zip(levels, levels[1:]))
    return to_gaps(out, base), tuple(unpacked[::-1])


def corner_source(base, label):
    # The first failed copy starts at c*=1+j iff c=28+j (mod 27).
    middle = 21 * base // 2000
    for low, u, v in product(range(3 * base // 1000, 3 * base // 1000 + 27),
                             range(middle, middle + 9), range(middle, middle + 9)):
        rhs = base - 4 * low - 2 * u - 2 * v - 3
        if (low - 28) % 27 == label and rhs % 9 == 0:
            assert base < 500 * low < 2 * base
            assert all(100 * a > base and 90 * a < base for a in (u, v))
            anchor = rhs // 9
            g = tuple(anchor + offset for offset in
                      (0, 0, 0, low + u + v + 4, 1, low + u, low + v, low, -2))
            assert sum(g) == base
            return low, g
    raise AssertionError("No representative of the residue on the open chart")


def nominal_update(gaps):
    # Equation (7), using only the primary balance and borrow.
    ordered = [gaps[0]] + [0] * (len(gaps) - 1)
    length = 2 * len(gaps) - 1
    for slot, value in enumerate(gaps[1:], 1):
        ordered[min(2 * slot, length - 2 * slot)] = value
    z = compare(gaps, 1)
    h = max(0, -1 - z)
    e = (z > 0) - (z < 0)
    for slot, adjustment in enumerate((e - h, -2 * e + h, e + h, -h)):
        ordered[slot] += adjustment
    return tuple(ordered)


def corner():
    base = 130243
    tally = Counter()
    simultaneous = []
    updates = 0
    for label in range(27):
        low, gaps = corner_source(base, label)
        full_copies = 1 + (low - 28) // 27
        assert (2 * low - 56) % 54 == 2 * label
        assert low - 27 * full_copies == 1 + label
        for copy in range(full_copies + 1):
            beginning = gaps
            for phase in range(244):
                failed = tuple(i for i in (2, 3, 4) if compare(gaps, i) < -1)
                successor, _ = literal_step(gaps, base)
                updates += 1
                if failed:
                    assert copy == full_copies
                    assert successor != nominal_update(gaps)
                    tally[failed] += 1
                    if len(failed) > 1:
                        simultaneous.append((1 + label, phase))
                    break
                assert successor == nominal_update(gaps)
                gaps = successor
            else:
                assert copy < full_copies
                shift = (12, 12, 12, -15, 12, -15, -15, -15, 12)
                assert gaps == tuple(value + delta for value, delta in zip(beginning, shift))
                continue
            break
    assert tally == {(2,): 19, (4,): 5, (2, 4): 3}
    assert simultaneous == [(8, 100), (26, 228), (27, 228)]
    weights = {failed: Fraction(tally[failed], 27) for failed in sorted(tally)}
    print(f"17-digit corner: B={base}, {updates} literal updates, "
          f"counts {dict(sorted(tally.items()))}, weights {weights}, simultaneous {simultaneous}")


def chamber_fixed():
    base = 1505
    # Raise the nonroot slot 3, with the four root slots kept equal.
    gaps = tuple(base * value // 35 for value in (6, 6, 6, 11, 6))
    assert sum(gaps) == base and 6 * min(gaps) >= base + 1
    assert compare(gaps, 1) == 0 and compare(gaps, 2) == 5 * base // 35
    successor, word = literal_step(gaps, base)
    assert successor == gaps
    subtraction = word_value(word, base)
    assert word_value(sorted(word, reverse=True), base) - word_value(sorted(word), base) == subtraction
    q = from_gaps(gaps)
    radius = max(abs(Fraction(value, base) - Fraction(5 - i, 5))
                 for i, value in enumerate(q, 1))
    assert radius == Fraction(3, 35)
    print(f"9-digit C fixed point: B={base}, gaps={gaps}, q={q}, "
          f"d_2={compare(gaps, 2)}, radius={radius}, word={word}, integer={subtraction}")


if __name__ == "__main__":
    prime_rows()
    p4_scalar()
    corner()
    chamber_fixed()
