"""Exact checks of the P10 prime rows, quadratic endpoint and finite sieve."""

from collections import Counter
from fractions import Fraction
from itertools import accumulate
from math import gcd, isqrt, prod


def primes_below(limit):
    prime = bytearray(b"\x01") * limit
    prime[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit - 1) + 1):
        if prime[divisor]:
            start = divisor * divisor
            prime[start::divisor] = b"\x00" * len(prime[start::divisor])
    return [n for n in range(2, limit) if prime[n]]


def slot_orbits(p):
    """Name a folded permutation orbit by its least positive slot."""
    k = (p - 1) // 2
    labels = [1] + [0] * k
    for first in range(1, k + 1):
        if labels[first] != 0:
            continue
        slot = first
        while labels[slot] == 0:
            labels[slot] = first
            slot = min(2 * slot, p - 2 * slot)
        assert slot == first
    return labels


def quotient_cosets(p):
    """Construct multiplicative cosets on all nonzero residues."""
    powers = [1]
    value = 2
    while value != 1:
        powers.append(value)
        value = 2 * value % p
    subgroup = set(powers) | {p - value for value in powers}
    labels = [0] * p
    for rep in range(1, p):
        if labels[rep] != 0:
            continue
        for multiplier in subgroup:
            labels[rep * multiplier % p] = rep
    return labels, (p - 1) // len(subgroup), len(powers)


def compare_walk_rows(p, slots, cosets):
    k = (p - 1) // 2
    assert slots[1:] == cosets[1:k + 1]
    literal = Counter()
    walk = Counter()
    for i in range(1, (k + 1) // 2 + 1):
        # Grow the two disjoint blocks in P10's definition.
        literal[slots[k + 1 - i]] += 1
        literal[slots[i - 1]] -= 1
        # Grow the alternating walk in its integer order of residues.
        if i > 1:
            walk[cosets[2 * i - 2]] -= 1
        walk[cosets[2 * i - 1]] += 1
        # The full slot row sums to zero. Its identity coordinate is
        # one less than the walk's, for the anchored first term.
        literal[1] += 1
        assert literal == walk, (p, i)
        literal[1] -= 1
        assert sum(literal.values()) == 0
    return (k + 1) // 2 - 1


def reduced_class_number(p):
    """Count Gauss-reduced primitive forms of discriminant -4p."""
    count = 0
    bound = isqrt(4 * p // 3)
    for z in range(bound // 2 + 1):
        b = 2 * z
        for a in range(max(1, b), bound + 1):
            c, remainder = divmod(p + z * z, a)
            if remainder or a > c or gcd(gcd(a, b), c) != 1:
                continue
            # Count both signs of b unless the form is on a boundary.
            count += 1 if b == 0 or b == a or a == c else 2
    return count


def compare_quadratic_rows(p, slots):
    k = (p - 1) // 2
    endpoint = k // 2
    assert p % 8 == 1
    chi = [-1] * p
    chi[0] = 0
    for n in range(1, k + 1):
        chi[n * n % p] = 1
    assert all((label == 1) == (chi[j] == 1)
               for j, label in enumerate(slots[1:], 1))
    prefix = list(accumulate(chi))
    w = [k // 2 if label == 1 else -(k // 2 + 1) for label in slots]
    assert sum(w) == 0
    w_prefix = list(accumulate(w, initial=0))
    assert prefix[k] == 0
    class_number = reduced_class_number(p)
    assert prefix[endpoint] * 2 == class_number
    # This weighted character sum is independent of reduced-form counting.
    weighted = sum(
        n * chi[n % p] * (1 if n % 4 == 1 else -1)
        for n in range(1, 4 * p, 2)
    )
    h, remainder = divmod(-weighted, 4 * p)
    assert remainder == 0 and h == class_number
    raw, legendre = [], []
    for i in range(1, endpoint + 1):
        row = w_prefix[-1] - w_prefix[k + 1 - i] - w_prefix[i]
        value = prefix[2 * i - 1] - 2 * prefix[i - 1]
        assert value % 2 == 1
        assert 2 * row == (k + 1) * (value - 1), (p, i)
        raw.append(row)
        legendre.append(value)
    assert raw[0] == 0 and legendre[1] == chi[3]
    assert legendre[-1] == 1 - h and raw[-1] < 0
    weak = max(raw[1:]) <= 0 or min(raw[1:]) >= 0
    strict = max(raw[1:]) < 0 or min(raw[1:]) > 0
    assert weak == (max(legendre[1:]) <= 1)
    assert strict == (max(legendre[1:]) <= -1)
    assert not strict or chi[3] == -1
    return weak, strict, endpoint - 1


def prime_rows():
    prime_counts = Counter()
    walk_rows = quadratic_rows = 0
    index_two = weak_count = 0
    strict_primes = []
    for p in primes_below(20_000):
        if p < 7:
            continue
        cosets, index, order = quotient_cosets(p)
        assert (index == 2) == (p % 8 == 1 and order == (p - 1) // 2)
        if p < 6_000 or index == 2:
            slots = slot_orbits(p)
            assert len(set(slots)) == index
        if p < 6_000:
            prime_counts[index] += 1
            walk_rows += compare_walk_rows(p, slots, cosets)
        if index == 2:
            index_two += 1
            weak, strict, rows = compare_quadratic_rows(p, slots)
            quadratic_rows += rows
            weak_count += weak
            if strict:
                strict_primes.append(p)
    print(f"walk: 7 <= p < 6000, {sum(prime_counts.values())} primes, "
          f"{prime_counts[1]} full-orbit, {walk_rows} secondary rows")
    print(f"index-two conversion: all primes 7 <= p < 20000")
    print(f"quadratic: p < 20000, {index_two} primes, "
          f"{quadratic_rows} secondary rows, {index_two} class numbers, "
          f"{weak_count} weak, strict {strict_primes}")


def finite_sieve(last_input):
    odd_primes = primes_below(last_input + 1)[1:]
    parity_masks = []
    for n in range(1, last_input + 1):
        mask = 0
        for position, prime in enumerate(odd_primes):
            while n % prime == 0:
                n //= prime
                mask ^= 1 << position
        parity_masks.append(mask)

    configurations = 1 << len(odd_primes)
    parity = [0] * configurations
    for mask in range(1, configurations):
        parity[mask] = parity[mask // 2] ^ (mask % 2)
    nonpositive = nonnegative = all_zero = 0
    for signs in range(configurations):
        values = [1 - 2 * parity[signs & mask] for mask in parity_masks]
        prefix = list(accumulate(values, initial=0))
        alternating = list(accumulate((1 - 2 * (j % 2)) * value
                                      for j, value in enumerate(values)))
        defects = [prefix[2 * i - 1] - 2 * prefix[i - 1] - 1
                   for i in range(2, (last_input + 1) // 2 + 1)]
        assert defects == [value - 1 for value in alternating[2::2]]
        nonpositive += max(defects) <= 0
        nonnegative += min(defects) >= 0
        all_zero += max(defects) == min(defects) == 0
    one_sided = nonpositive + nonnegative - all_zero
    print(f"sieve: input {last_input}, nonpositive {nonpositive}, "
          f"nonnegative {nonnegative}, all zero {all_zero}, "
          f"one-sided {one_sided}/{configurations}")

    if last_input == 15:
        modulus = 8 * prod(odd_primes)
        frequencies = Counter()
        for residue in range(1, modulus, 8):
            if gcd(residue, modulus) != 1:
                continue
            pattern = 0
            for position, prime in enumerate(odd_primes):
                if pow(residue, (prime - 1) // 2, prime) == prime - 1:
                    pattern |= 1 << position
            frequencies[pattern] += 1
        assert len(frequencies) == configurations
        assert len(set(frequencies.values())) == 1
        print(f"CRT: {configurations} sign patterns, each with "
              f"{frequencies[0]} reduced classes modulo {modulus}")


def artin_lower_bound():
    product = prod(1 - Fraction(1, prime * (prime - 1))
                   for prime in primes_below(20))
    lower = product * Fraction(18, 19)
    assert lower > Fraction(1, 3)
    print(f"Artin constant lower bound: {lower} > 1/3")


if __name__ == "__main__":
    prime_rows()
    for cutoff in (15, 31, 47):
        finite_sieve(cutoff)
    artin_lower_bound()
