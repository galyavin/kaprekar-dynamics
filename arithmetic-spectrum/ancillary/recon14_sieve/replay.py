"""Independent exact replay of the index-three and index-four census and prefix checks."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import factorial
from pathlib import Path
import json
import random


SOURCE = Path(__file__).resolve().parent.parent / "arith"
SHORT_INPUTS = (7, 11, 19, 23)
LONG_INPUTS = (101, 301, 1001)
BANDS = ((7, 1599), (1601, 4999), (5001, 9999))


def smallest_factors(limit):
    factors = [0] * (limit + 1)
    for prime in range(2, limit + 1):
        if factors[prime] == 0:
            for multiple in range(prime, limit + 1, prime):
                if factors[multiple] == 0:
                    factors[multiple] = prime
    return factors


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def interval(inequalities, lower, upper):
    """Intersect bounds with constraints of the form slope * u >= constant."""
    for slope, constant in inequalities:
        if slope > 0:
            lower = max(lower, Fraction(constant, slope))
        elif slope < 0:
            upper = min(upper, Fraction(constant, slope))
        elif constant > 0:
            return None
        if lower > upper:
            return None
    return lower, upper


def feasible_square(inequalities):
    """Maintain a feasible point, repairing it on each newly violated line."""
    accepted = [(1, 0, -1), (-1, 0, -1), (0, 1, -1), (0, -1, -1)]
    x, y = 0, 0
    for a, b, constant in inequalities:
        if a * x + b * y < constant:
            if a == b == 0:
                return False
            # Every feasible repair can be joined to the previous witness by a
            # segment in the old polygon, so a repair exists on this line.
            if a:
                origin = (Fraction(constant, a), 0)
                direction = (-Fraction(b, a), 1)
            else:
                origin = (0, Fraction(constant, b))
                direction = (1, 0)
            bounds = interval(
                ((dot((c, d), direction), rhs - dot((c, d), origin))
                 for c, d, rhs in accepted),
                -1, 1,
            )
            if bounds is None:
                return False
            low, high = bounds
            parameter = max(low, min(high, 0))
            x = origin[0] + parameter * direction[0]
            y = origin[1] + parameter * direction[1]
        accepted.append((a, b, constant))
    return True


def half_space(rows, dim, strict=False):
    # A nonzero normal can be scaled to one of the cube faces. For the
    # strict test maximize the margin epsilon in [0, 1] on the full cube.
    # At a vertex, Cramer's rule bounds its denominator by (d+1)! M**d,
    # since the margin column has entries of modulus at most one.
    scale = 1
    if strict:
        if not all(any(row) for row in rows):
            return False
        magnitude = max(1, max(abs(entry) for row in rows for entry in row))
        scale = factorial(dim + 1) * magnitude ** dim
    margin = int(strict)
    for axis, sign in product(range(dim), (-1, 1)):
        constraints = [
            tuple(scale * entry for i, entry in enumerate(row) if i != axis)
            + (margin - scale * sign * row[axis],)
            for row in rows
        ]
        if dim == 2:
            if interval(constraints, -1, 1) is not None:
                return True
        elif feasible_square(constraints):
            return True
    return False


def walk(labels, t):
    current = [0] * t
    rows = []
    for n in range(1, len(labels)):
        current[labels[n]] += 1 if n % 2 else -1
        if n >= 3 and n % 2:
            rows.append(tuple(current[1:]))
    return rows


def orbit_labels(prime):
    labels = [0] + [-1] * ((prime - 1) // 2)
    t = 0
    for initial in range(1, len(labels)):
        if labels[initial] >= 0:
            continue
        slot = initial
        while labels[slot] < 0:
            labels[slot] = t
            doubled = 2 * slot
            slot = min(doubled, prime - doubled)
        assert slot == initial
        t += 1
    return t, labels


def literal_rows(labels, t):
    running = [0] * t
    rows = []
    for i in range(1, len(labels) // 2 + 1):
        running[labels[-i]] += 1
        running[labels[i - 1]] -= 1
        if i >= 2:
            rows.append(tuple(running[1:]))
    return rows


def records(path):
    return {record["n"]: record for record in map(json.loads, path.read_text().splitlines())}


def prime_replay(t, primes, frozen, supplemental):
    counts = Counter()
    initial = Counter()
    bands = Counter()
    totals = Counter()
    prefix = Counter()
    row_total = 0
    observed = set()
    for prime in primes:
        index, orbit = orbit_labels(prime)
        if index % t:
            continue
        assert index == frozen[prime]["t"]
        if index == t:
            observed.add(prime)
            literal = literal_rows(orbit, t)
            assert literal == walk(orbit, t)
            row_total += len(literal)
            weak = half_space(literal, t - 1)
            strict = weak and half_space(literal, t - 1, strict=True)
            classification = "strict" if strict else "scarred" if weak else "zero"
            expected = supplemental.get(prime, frozen[prime])["cls"]
            assert classification == expected, (prime, classification, expected)
            initial[frozen[prime]["cls"]] += 1
            counts[classification] += 1
            band = next((low, high) for low, high in BANDS if low <= prime <= high)
            bands[band, "total"] += 1
            bands[band, "nonzero"] += weak

        max_input = min(LONG_INPUTS[-1], 2 * ((prime + 1) // 4) - 1)
        # This also generates labels on the larger divisible-index superset.
        residues = [pow(n, (prime - 1) // t, prime) for n in range(max_input + 1)]
        residue_label = {1: 0}
        labels = [0] + [residue_label.setdefault(r, len(residue_label)) for r in residues[1:]]
        quotient = walk(labels, t)
        if index == t:
            relabel = {orbit[n]: labels[n] for n in range(1, max_input + 1)}
            assert len(relabel) == len(set(relabel.values()))
            assert all(relabel[orbit[n]] == labels[n] for n in range(1, max_input + 1))
        for cutoff in (*SHORT_INPUTS, *LONG_INPUTS):
            if max_input < cutoff:
                continue
            survives = half_space(quotient[:(cutoff - 1) // 2], t - 1)
            totals[cutoff, "divisible"] += 1
            prefix[cutoff, "divisible"] += survives
            if index == t:
                totals[cutoff, "exact"] += 1
                prefix[cutoff, "exact"] += survives

    assert observed == {p for p in primes if frozen[p]["t"] == t}
    print(f"t={t}: {len(observed)} primes, {row_total} literal/quotient rows matched")
    classes = ("strict", "scarred", "zero", "unknown")
    print("  Classes (strict, scarred, zero, unknown)")
    print(f"  initial {tuple(initial[c] for c in classes)}")
    print(f"  final   {tuple(counts[c] for c in classes)}")
    for low, high in BANDS:
        print(f"  p={low:4d}..{high}: {bands[(low, high), 'nonzero']}/{bands[(low, high), 'total']} nonzero")
    return totals, prefix


def model_replay(t, cutoffs):
    limit = max(cutoffs)
    factors = smallest_factors(limit)
    odd_primes = [p for p in range(3, limit + 1) if factors[p] == p]
    if cutoffs == SHORT_INPUTS:
        assignments = product(range(t), repeat=len(odd_primes))
        total = t ** len(odd_primes)
    else:
        # Replay the recorded input draw stream with an independent solver.
        total = 1024
        generator = random.Random(209034 + t)
        assignments = ([generator.randrange(t) for _ in odd_primes] for _ in range(total))
    successes = Counter()
    for assignment in assignments:
        prime_labels = {2: 0, **dict(zip(odd_primes, assignment))}
        labels = [0] * (limit + 1)
        for n in range(2, limit + 1):
            factor = factors[n]
            labels[n] = (prime_labels[factor] + labels[n // factor]) % t
        rows = walk(labels, t)
        for cutoff in cutoffs:
            if not half_space(rows[:(cutoff - 1) // 2], t - 1):
                break
            successes[cutoff] += 1
    return {cutoff: (successes[cutoff], total) for cutoff in cutoffs}


def main():
    frozen = records(SOURCE / "certificates.jsonl")
    supplemental = records(SOURCE / "planar-certificates.jsonl")
    factors = smallest_factors(9999)
    primes = [n for n in range(7, 10_000, 2) if factors[n] == n]
    for t in (3, 4):
        totals, prefix = prime_replay(t, primes, frozen, supplemental)
        model = model_replay(t, SHORT_INPUTS) | model_replay(t, LONG_INPUTS)
        print("  Input | model | exact-index primes | divisible-index primes")
        for cutoff in (*SHORT_INPUTS, *LONG_INPUTS):
            model_pass, model_total = model[cutoff]
            exact_pass, exact_total = prefix[cutoff, "exact"], totals[cutoff, "exact"]
            divides_pass, divides_total = prefix[cutoff, "divisible"], totals[cutoff, "divisible"]
            print(f"  {cutoff:5d} | {model_pass}/{model_total} | {exact_pass}/{exact_total} | {divides_pass}/{divides_total}")


if __name__ == "__main__":
    main()
