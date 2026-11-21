"""Check the orbit-linear secondary rows directly from their two end sums."""

from math import isqrt


def secondary_rows(n):
    k = (n - 1) // 2
    # Label every positive slot by its forward orbit. The anchor shares
    # the orbit label of 1 because these profiles are constant on J.
    labels = bytearray(k + 1)
    orbits = []
    for first in range(1, k + 1):
        if labels[first]:
            continue
        index = len(orbits) + 1
        slot = first
        size = 0
        while not labels[slot]:
            labels[slot] = index
            size += 1
            slot = min(2 * slot, n - 2 * slot)
        assert slot == first
        orbits.append((first, size))
    labels[0] = labels[1]
    coefficients = [0] * (len(orbits) + 1)
    zeros = []
    small = {}
    final_row = (k + 1) // 2
    for i in range(1, final_row + 1):
        coefficients[labels[k + 1 - i]] += 1
        coefficients[labels[i - 1]] -= 1
        if i in (2, 3, 4, 31, 34):
            small[i] = tuple(coefficients[2:])
        if i >= 2 and not any(coefficients):
            zeros.append(i)
    print(f"n={n}, positive orbits (first, size)={orbits}")
    print(f"secondary rows={final_row - 1}, zero indices={zeros}")
    return small, zeros, orbits


def main():
    p = 1069
    assert all(p % d for d in range(2, isqrt(p) + 1))
    assert pow(2, 178, p) == p - 1
    assert pow(2, 356, p * p) == 1 + 417 * p
    expected = {2: (1, 0), 3: (2, 0), 4: (1, 1),
                31: (-3, 0), 34: (-1, -1)}
    small, zeros, orbits = secondary_rows(p)
    assert small == expected and not zeros
    assert orbits == [(1, 178), (3, 178), (7, 178)]
    small, zeros, orbits = secondary_rows(p * p)
    assert small == {i: coeff + (0, 0, 0) for i, coeff in expected.items()}
    assert [size for _, size in orbits] == [178 * p] * 3 + [178] * 3
    assert zeros == [534, 1070, 1603]


if __name__ == "__main__":
    main()
