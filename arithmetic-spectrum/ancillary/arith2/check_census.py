"""Check the prime walk and compare exact finite prefixes with Task 207's census."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path
import random


SOURCE = Path(__file__).resolve().parent.parent / "arith"
PREFIXES = (7, 11, 19, 23)
SAMPLE_PREFIXES = (101, 301, 1001)
SAMPLE_SIZE = 1024
SAMPLE_SEED = 209034


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def cross(left, right):
    a, b, c = left
    d, e, f = right
    return b * f - c * e, c * d - a * f, a * e - b * d


def survives(rows, dim):
    points = set(filter(any, rows))
    if not points:
        return True
    if dim == 2:
        normals = ((b, -a) for a, b in points)
    else:
        normals = (
            vector for left, right in combinations(points, 2)
            if any(vector := cross(left, right))
        )
    have_normal = False
    for normal in normals:
        have_normal = True
        projections = [dot(point, normal) for point in points]
        if min(projections) >= 0 or max(projections) <= 0:
            return True
    return not have_normal


def walk(labels, dim):
    point = [0] * (dim + 1)
    for n, label in enumerate(labels, 1):
        point[label] += 1 if n % 2 else -1
        if n > 1 and n % 2:
            yield tuple(point[1:])


def check_walk_and_certificate(record):
    p, t = record["n"], record["t"]
    k = (p - 1) // 2
    subgroup = {1, p - 1}
    value = 2
    while value not in subgroup:
        subgroup.update((value, p - value))
        value = value * 2 % p

    coset_labels = [0] * p
    count = 0
    covered = set()
    for n in range(1, k + 1):
        if n in covered:
            continue
        coset = {n * h % p for h in subgroup}
        for a in coset:
            coset_labels[a] = count
        covered.update(coset)
        count += 1
    assert count == t and len(covered) == p - 1 and len(subgroup) == 2 * record["r"]

    # Compute the P10 slot labels by folded successor, independently of cosets.
    slot_labels = [0] + [-1] * k
    count = 0
    for n in range(1, k + 1):
        if slot_labels[n] >= 0:
            continue
        j = n
        while slot_labels[j] == -1:
            slot_labels[j] = count
            j = min(2 * j, p - 2 * j)
        assert j == n
        count += 1
    assert count == t

    horizon = (k + 1) // 2
    rows = list(walk(coset_labels[1:2 * horizon], t - 1))
    literal = [0] * t
    for i in range(1, horizon + 1):
        literal[slot_labels[k + 1 - i]] += 1
        literal[slot_labels[i - 1]] -= 1
        if i >= 2:
            assert rows[i - 2] == tuple(literal[1:]), (p, i)

    cert = record["cert"]
    cls = record["cls"]
    if cls in ("strict", "scarred"):
        assert any(cert["x"]) and len(cert["x"]) == t - 1
        assert min(dot(row, cert["x"]) for row in rows) >= (cls == "strict")
    if "tie" in cert:
        assert cls == "scarred" and not any(rows[cert["tie"] - 2])
    if "relation" in cert:
        witness_rows = []
        totals = [0] * (t - 1)
        for i, weight in cert["relation"]:
            assert weight > 0 and 2 <= i <= horizon
            row = rows[i - 2]
            witness_rows.append(row)
            for j, coordinate in enumerate(row):
                totals[j] += weight * coordinate
        assert not any(totals)
        if cls == "zero":
            if t == 3:
                assert any(a * d != b * c for (a, b), (c, d) in combinations(witness_rows, 2))
            else:
                assert any(dot(cross(a, b), c) for a, b, c in combinations(witness_rows, 3))
    elif cls == "zero":
        raise AssertionError(f"Missing spanning relation at {p}")
    end = 2 * horizon - 1
    assert prefix_decisions(rows, t - 1, (end,))[end] == (cls != "zero"), p
    return len(rows)


def smallest_prime_factors(largest):
    spf = [0] * (largest + 1)
    for n in range(2, largest + 1):
        if spf[n] == 0:
            for multiple in range(n, largest + 1, n):
                if spf[multiple] == 0:
                    spf[multiple] = n
    return spf


def model_walk(t, spf, primes, pattern):
    labels = [0] * len(spf)
    for q, label in zip(primes, pattern):
        labels[q] = label
    for n in range(3, len(labels)):
        if spf[n] != n:
            labels[n] = (labels[spf[n]] + labels[n // spf[n]]) % t
    return walk(labels[1:], t - 1)


def exact_model(t):
    largest = max(PREFIXES)
    spf = smallest_prime_factors(largest)
    primes = [n for n in range(3, largest + 1) if spf[n] == n]
    counts = Counter()
    for pattern in product(range(t), repeat=len(primes)):
        rows = list(model_walk(t, spf, primes, pattern))
        for n in PREFIXES:
            if survives(rows[:(n - 1) // 2], t - 1):
                counts[n] += 1
            else:
                break
    return {n: Fraction(counts[n], t ** len(primes)) for n in PREFIXES}


def clip_face(vertices, inequality):
    """Intersect a face polygon with one inequality using homogeneous coordinates."""
    values = [dot(vertex, inequality) for vertex in vertices]
    if min(values) >= 0:
        return vertices
    if max(values) < 0:
        return []
    result = []
    for j, vertex in enumerate(vertices):
        previous, before, current = vertices[j - 1], values[j - 1], values[j]
        if before * current < 0:
            intersection = tuple(abs(current) * a + abs(before) * b for a, b in zip(previous, vertex))
            divisor = gcd(*intersection)
            result.append(tuple(coordinate // divisor for coordinate in intersection))
        if current >= 0:
            result.append(vertex)
    return list(dict.fromkeys(result))


def prefix_decisions(rows, dim, cutoffs):
    # A nonzero feasible normal can be scaled onto one face of [-1,1]^dim.
    # In dimension two the unused polygon coordinate has coefficient zero.
    square = [(-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
    faces = {(axis, sign): square for axis in range(dim) for sign in (-1, 1)}
    seen = set()
    decisions = {}
    for n, row in enumerate(rows, 1):
        if row not in seen and any(row) and faces:
            seen.add(row)
            for axis, sign in list(faces):
                free = [row[i] for i in range(dim) if i != axis] + [0]
                inequality = (*free[:2], sign * row[axis])
                polygon = clip_face(faces[axis, sign], inequality)
                if polygon:
                    faces[axis, sign] = polygon
                else:
                    del faces[axis, sign]
        if (cutoff := 2 * n + 1) in cutoffs:
            decisions[cutoff] = bool(faces)
    return decisions


def sampled_model(t):
    largest = max(SAMPLE_PREFIXES)
    spf = smallest_prime_factors(largest)
    primes = [n for n in range(3, largest + 1) if spf[n] == n]
    generator = random.Random(SAMPLE_SEED + t)
    counts = Counter()
    for draw in range(SAMPLE_SIZE):
        pattern = [generator.randrange(t) for _ in primes]
        rows = list(model_walk(t, spf, primes, pattern))
        decisions = prefix_decisions(rows, t - 1, (*SAMPLE_PREFIXES, *PREFIXES))
        checks = (*PREFIXES, 101) if draw < 12 else PREFIXES
        for n in checks:
            assert decisions[n] == survives(rows[:(n - 1) // 2], t - 1)
        for n in SAMPLE_PREFIXES:
            counts[n] += decisions[n]
    return counts


def arithmetic_prefix(p, t, n):
    labels = {1: 0}
    for a in range(1, n + 1):
        value = pow(a, (p - 1) // t, p)
        if value not in labels:
            labels[value] = len(labels)
        yield labels[value]


def format_classes(records):
    counts = Counter(row["cls"] for row in records)
    return tuple(counts[name] for name in ("strict", "scarred", "zero", "unknown"))


def main():
    data = [json.loads(line) for line in (SOURCE / "certificates.jsonl").read_text().splitlines()]
    supplement = {
        row["n"]: row
        for row in map(json.loads, (SOURCE / "planar-certificates.jsonl").read_text().splitlines())
    }
    raw = (SOURCE / "census.txt").read_text().splitlines()
    primes = [row for row in data if row["band"] == "census" and row["bin"].startswith("P")]
    expected_primes = [
        p for p in range(7, 10000, 2)
        if all(p % q for q in range(2, isqrt(p) + 1))
    ]
    assert [row["n"] for row in primes] == expected_primes
    for t in (3, 4):
        initial = [row for row in primes if row["t"] == t]
        registered = next(line for line in raw if line.startswith(f"census exact P t={t} "))
        strict, scarred, full, zero, unknown, *_ = map(int, registered.split()[4:])
        assert full == 0 and (strict, scarred, zero, unknown) == format_classes(initial)
        verified = [supplement.get(row["n"], row) for row in initial]
        assert all(row["cls"] != "unknown" for row in verified)
        row_count = sum(map(check_walk_and_certificate, verified))
        print(f"t={t} {len(verified)} primes, {row_count} secondary rows replayed")
        print(f"  initial (strict, scarred, zero, unknown): {format_classes(initial)}")
        print(f"  final   (strict, scarred, zero, unknown): {format_classes(verified)}")
        for low, high in ((7, 1599), (1601, 4999), (5001, 9999)):
            band = [row for row in verified if low <= row["n"] <= high]
            strict, scarred, _, _ = format_classes(band)
            count = strict + scarred
            print(f"  {low}–{high}: {count}/{len(band)} nonzero ({count / len(band):.2%})")
        model = exact_model(t)
        sampled = sampled_model(t)
        print("  N | model (exact through 23, sampled thereafter) | exact-index prime prefix | divisible-index prime prefix")
        for n in (*PREFIXES, *SAMPLE_PREFIXES):
            eligible = [row for row in primes if row["t"] % t == 0 and row["n"] >= 2 * n + 1]
            counts = Counter()
            for row in eligible:
                p = row["n"]
                rows = walk(arithmetic_prefix(p, t, n), t - 1)
                survives_prefix = prefix_decisions(rows, t - 1, (n,))[n]
                counts["divisible_size"] += 1
                counts["divisible_alive"] += survives_prefix
                if row["t"] == t:
                    counts["exact_size"] += 1
                    counts["exact_alive"] += survives_prefix
                    if row["cls"] in ("strict", "scarred") or p in supplement:
                        assert survives_prefix, (p, n)
            prime_exact = f'{counts["exact_alive"]}/{counts["exact_size"]}'
            prime_divisible = f'{counts["divisible_alive"]}/{counts["divisible_size"]}'
            share = model[n] if n in model else Fraction(sampled[n], SAMPLE_SIZE)
            fraction = str(share) if n in model else f"{sampled[n]}/{SAMPLE_SIZE}"
            print(f"  {n:4} | {fraction} ({float(share):.2%}) | {prime_exact} | {prime_divisible}")


if __name__ == "__main__":
    main()
