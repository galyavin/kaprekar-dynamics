"""Exact folded-orbit, row, and literal-subtraction controls for Task 208."""

from collections import Counter
from functools import cache
from math import gcd, lcm


def factor(n):
    factors = Counter()
    p = 3
    while p * p <= n:
        while n % p == 0:
            factors[p] += 1
            n //= p
        p += 2
    if n > 1:
        factors[n] += 1
    return factors


@cache
def local_order(d):
    """Least positive exponent reaching either sign, and its sign."""
    exponent, residue = 1, 2 % d
    while residue not in (1, d - 1):
        exponent += 1
        residue = 2 * residue % d
    return exponent, 1 if residue == 1 else -1


def unit_orbit_count(d):
    phi = d
    for p in factor(d):
        phi -= phi // p
    order, _ = local_order(d)
    assert phi % (2 * order) == 0
    return phi // (2 * order)


def classified_stratum(factors):
    """Use the prime-power lifting and CRT criteria, away from modulus n."""
    if len(factors) == 1:
        (p, exponent), = factors.items()
        if exponent in (2, 3) and unit_orbit_count(p) == 1:
            if pow(2, p - 1, p * p) != 1:
                return exponent
    if len(factors) == 2 and all(exponent == 1 for exponent in factors.values()):
        p, q = factors
        r_p, sign_p = local_order(p)
        r_q, sign_q = local_order(q)
        if r_p == (p - 1) // 2 and r_q == (q - 1) // 2:
            exponent = lcm(r_p, r_q)
            if gcd(r_p, r_q) == 1 and sign_p**(exponent // r_p) != sign_q**(exponent // r_q):
                return 3
    return None


def orbit_partition(n):
    orbits = []
    seen = set()
    for j in range(1, (n + 1) // 2):
        if j in seen:
            continue
        orbit, slot = [], j
        while slot not in seen:
            seen.add(slot)
            orbit.append(slot)
            slot = min(2 * slot, n - 2 * slot)
        assert slot == j
        orbits.append(tuple(orbit))
    return orbits


def row_values(slots):
    values, running = [], 0
    for i in range(1, len(slots) // 2 + 1):
        running += slots[-i] - slots[i - 1]
        values.append(running)
    return tuple(values)


def prefix_values(slots):
    values, running = [], 0
    for j in range(len(slots) // 2):
        running += slots[2 * j + 1] - slots[j]
        values.append(running)
    return tuple(values)


def number(digits, base):
    value = 0
    for digit in digits:
        assert 0 <= digit < base
        value = base * value + digit
    return value


def fixed_point_control(n, marked, scale):
    m = len(marked)
    c = sum(marked)
    base = scale * (m * m + c)
    gaps = [scale * (m + entry) for entry in marked]
    assert (m + 1) * min(gaps) >= base + 1
    assert row_values(gaps)[0] == 0 and min(row_values(gaps)[1:]) >= 0
    differences, remaining = [], base
    for gap in gaps[:-1]:
        remaining -= gap
        differences.append(remaining)

    seed = differences + [0] * m
    value = number(seed, base) - number(reversed(seed), base)
    digits = []
    for _ in range(n):
        value, digit = divmod(value, base)
        digits.append(digit)
    assert value == 0
    digits.reverse()
    k = m - 1
    levels = list(zip(digits[:k], digits[k + 1:]))
    assert digits[k] == base - 1 and differences[-1] > 0
    assert all(min(a) >= max(b) for a, b in zip(levels, levels[1:]))
    sorted_digits = sorted(digits)
    assert [sorted_digits[-1 - i] - sorted_digits[i] for i in range(k)] == differences
    assert m * differences[0] - (m - 1) * base == scale * c > 0


def main():
    expected_strict = {9: (1,), 15: (1, 2, 1), 21: (1, 1, 1, 2)}
    observed_strict = {}
    low_counts, zero_counts = Counter(), Counter()
    composite_count = physical_count = indicator_count = row_checks = 0
    for n in range(9, 1002, 2):
        factors = factor(n)
        if sum(factors.values()) < 2:
            continue
        composite_count += 1
        m = (n + 1) // 2
        orbits = orbit_partition(n)
        layers = Counter(gcd(orbit[0], n) for orbit in orbits)
        for orbit in orbits:
            assert all(gcd(j, n) == gcd(orbit[0], n) for j in orbit)
        for divisor in range(3, n + 1, 2):
            if n % divisor == 0:
                assert layers[n // divisor] == unit_orbit_count(divisor)

        t = len(orbits)
        assert (t if t <= 3 else None) == classified_stratum(factors)
        for p in factors:
            marked = [int(j > 0 and j % p == 0) for j in range(m)]
            rows = row_values(marked)
            assert rows == prefix_values(marked)
            assert set(rows) <= {0, 1} and 0 < sum(marked) < m
            indicator_count += 1
            if t <= 3:
                for scale in (1, 7, 1009):
                    fixed_point_control(n, marked, scale)
                    physical_count += 1

        if t > 3:
            continue
        low_counts[t] += 1
        assert layers[1] == 1
        if n in expected_strict:
            nonunits = [int(j > 0 and gcd(j, n) != 1) for j in range(m)]
            rows = row_values(nonunits)[1:]
            assert rows == expected_strict[n]
            observed_strict[n] = rows
            print(f"strict n={n}: orbits={orbits}, secondary rows={rows}")
            continue

        forced = 2 if n % 3 else 4
        for orbit in orbits[1:]:
            support = set(orbit)
            offset = [int(j in support) for j in range(m)]
            rows = row_values(offset)
            assert rows == prefix_values(offset)
            assert rows[forced - 1] == 0
            row_checks += 1
        zero_counts[forced] += 1

    assert observed_strict == expected_strict
    p = 1093
    assert factor(p) == {p: 1}
    assert pow(2, p - 1, p * p) == 1
    assert unit_orbit_count(p * p) == p * unit_orbit_count(p)
    print(f"index lift at p={p}: {unit_orbit_count(p)} -> {unit_orbit_count(p*p)}")
    print(f"odd composites through 1001: {composite_count}; t=2: {low_counts[2]}; t=3: {low_counts[3]}")
    print(f"forced d_2: {zero_counts[2]}; forced d_4: {zero_counts[4]}; basis checks: {row_checks}")
    print(f"prime-divisor indicators: {indicator_count}; literal physical fixed-point checks: {physical_count}")

    progression_primes = []
    for p in range(7, 334, 8):
        if p <= 7 or factor(p) != {p: 1}:
            continue
        assert (p - 1) // 2 % 2 == 1
        if local_order(p) == ((p - 1) // 2, 1):
            assert len(orbit_partition(3 * p)) == 3
            assert pow(2, (p - 1) // 2, p) == 1 and pow(2, (p - 1) // 2, 3) == 2
            progression_primes.append(p)
    print(f"progression controls p up to 333 with p>7: {progression_primes}")


if __name__ == "__main__":
    main()
