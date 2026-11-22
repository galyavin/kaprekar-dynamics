"""Replay Task 247 certificates with Python's standard library."""

from fractions import Fraction
from functools import reduce
from itertools import chain
from math import gcd, isqrt, lcm
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent
LIMIT = 10_000
SHORT_LIMIT = 101
EXAMPLES = {25: 1, 27: 2, 35: 0, 39: 5, 97: 0, 119: 23, 277: 91}


class Length:
    def __init__(self, n):
        self.k = (n - 1) // 2
        self.m = self.k + 1
        self.permutation = [min(2 * j, n - 2 * j) for j in range(self.m)]
        self.orbits = []
        self.labels = [-1] * self.m
        for start in range(1, self.m):
            if self.labels[start] >= 0:
                continue
            orbit = []
            label = len(self.orbits)
            j = start
            while self.labels[j] < 0:
                self.labels[j] = label
                orbit.append(j)
                j = self.permutation[j]
            assert j == start
            self.orbits.append(orbit)
        self.r = len(self.orbits[0])
        self.comp = tuple(sorted(chain.from_iterable(self.orbits[1:])))
        self.q = len(self.comp)
        self.s = len(self.orbits) - 1
        self.indices = range(2, self.m // 2 + 1)

    def invariant_rows(self):
        """Return literal-support flags and rows in complement-orbit offsets."""
        suffix_count = 0
        row = [0] * self.s
        result = {}
        for i in range(1, self.m // 2 + 1):
            for j, sign in ((i - 1, -1), (self.m - i, 1)):
                label = self.labels[j]
                if label > 0:
                    suffix_count += 1
                    row[label - 1] += sign
            if i >= 2:
                result[i] = (suffix_count > 0, tuple(row))
        return result

    def rotated_rows(self, indices):
        indices = tuple(indices)
        pos = self.comp
        for t in range(self.r):
            for i in indices:
                row = tuple(int(j >= self.m - i) - int(j < i) for j in pos)
                yield (t, i), row
            pos = tuple(self.permutation[j] for j in pos)
        assert pos == self.comp


def dot(row, point):
    return sum(a * b for a, b in zip(row, point))


def primitive(values):
    values = tuple(values)
    denominator = lcm(*(v.denominator for v in values))
    return primitive_integer(int(v * denominator) for v in values)


def primitive_integer(values):
    values = tuple(values)
    divisor = reduce(gcd, values, 0) or 1
    return tuple(v // divisor for v in values)


def echelon(rows):
    """Fraction-free echelon basis over Q with primitive integer rows."""
    basis = {}
    for row in rows:
        for index in sorted(basis):
            if row[index]:
                pivot = basis[index]
                row = primitive_integer(pivot[index] * a - row[index] * b
                                        for a, b in zip(row, pivot))
        index = next((j for j, val in enumerate(row) if val), None)
        if index is not None:
            basis[index] = row
    return basis


def kernel(equations, width):
    """Return a basis as vectors in the original coordinates."""
    basis = echelon(equations)
    columns = []
    for free in range(width):
        if free in basis:
            continue
        vector = [Fraction(0)] * width
        vector[free] = Fraction(1)
        for index, row in sorted(basis.items(), reverse=True):
            vector[index] = -dot(row, vector) / row[index]
        columns.append(tuple(vector))
    return columns


def relation_value(rows, weights):
    return tuple(sum(w * a for w, a in zip(weights, column))
                 for column in zip(*rows))


def cone_certificate(length, certificate):
    """Prove the implicit rows, their rank, and the reduced row forms."""
    relint, relation = certificate["interior"], certificate["relation"]
    assert len(relint) == length.s
    rows = length.invariant_rows()
    values = {i: dot(row, relint) for i, (_, row) in rows.items()}
    assert all(value >= 0 for value in values.values())
    implicit = tuple(i for i, value in values.items() if value == 0)
    assert len(relation) == len(implicit)
    assert all(value > 0 for value in relation)
    assert not any(relation_value([rows[i][1] for i in implicit], relation))
    equations = [row for _, row in length.rotated_rows(implicit)]
    span = kernel(equations, length.q)
    return implicit, span


def distinct_normals(length, implicit, span):
    rows = {}
    for _, row in length.rotated_rows(i for i in length.indices if i not in implicit):
        normal = primitive(dot(row, vector) for vector in span)
        assert any(normal)
        rows[normal] = None
    return tuple(rows)


def phase_zero_rows(length, implicit, span):
    return {
        i: primitive(dot(row, vector) for vector in span)
        for (t, i), row in length.rotated_rows(i for i in length.indices if i not in implicit)
        if t == 0
    }


def verify_faces(length, implicit, span, certificate):
    rows = distinct_normals(length, implicit, span)
    phase_zero = phase_zero_rows(length, implicit, span)
    results = {"facet": [], "nonfacet": [], "apex": []}
    assert set(map(int, certificate)) == set(phase_zero)
    for i, row in phase_zero.items():
        witness = certificate[str(i)]
        other = [v for v in rows if v != row]
        point = witness["point"]
        assert len(point) == len(span) and dot(row, point) == 0
        slacks = [dot(v, point) for v in other]
        assert all(slack >= 0 for slack in slacks)
        active = [v for v, slack in zip(other, slacks) if slack == 0]
        weights = witness["weights"]
        assert len(weights) == len(active) and all(w > 0 for w in weights)
        if active:
            assert witness["scale"] > 0
            assert relation_value(active, weights) == tuple(witness["scale"] * a for a in row)
        rank = len(echelon([row, *active]))
        if rank == 1:
            results["facet"].append(i)
        elif rank == len(span):
            results["apex"].append(i)
        else:
            results["nonfacet"].append(i)
    return results


def is_prime(n):
    return all(n % d for d in range(2, isqrt(n) + 1))


def verify_prime_walk(n, length, rows):
    """Compare the literal support and entries with the retained P11 walk."""
    walk = [0] * length.s
    touched_complement = False
    for b in range(1, length.k + 1):
        label = length.labels[b]
        if label > 0:
            walk[label - 1] += 1 if b % 2 else -1
            touched_complement = True
        if b % 2 and (b + 1) // 2 in rows:
            i = (b + 1) // 2
            assert rows[i] == (touched_complement, tuple(walk))
    if length.s == 1:
        prefix = 0
        for b in range(1, length.k + 1):
            chi = 1 if pow(b, length.k, n) == 1 else -1
            prefix += chi if b % 2 else -chi
            if b % 2 and (b + 1) // 2 in rows:
                assert rows[(b + 1) // 2][1] == ((1 - prefix) // 2,)
        assert rows[max(rows)][1][0] > 0


def verify_census(certificates):
    certs = {record[0]: record[1:] for record in certificates}
    assert len(certs) == len(certificates)
    candidates = set()
    classified = []
    prime_count = 0
    index_two = []
    for n in range(7, LIMIT, 2):
        length = Length(n)
        if not length.q or length.labels[3] != 0:
            continue
        candidates.add(n)
        rows = length.invariant_rows()
        assert not rows[2][0]
        prime = is_prime(n)
        if prime:
            verify_prime_walk(n, length, rows)
        case, data = certs[n]
        if case == "+":
            assert len(data) == length.s
            assert all(dot(row, data) > 0 for keep, row in rows.values() if keep)
            classified.append(n)
            prime_count += prime
            if prime and length.s == 1:
                index_two.append(n)
        else:
            assert case == "-" and data
            assert all(rows[i][0] and w > 0 for i, w in data)
            assert not any(relation_value([rows[i][1] for i, _ in data], [w for _, w in data]))
    assert set(certs) == candidates
    assert (len(classified), prime_count) == (274, 94)
    print(f"Census: {len(candidates)} candidates certified, {len(classified)} rank zero "
          f"({prime_count} prime, {len(classified) - prime_count} composite)")
    print(f"Additional index-two primes below {LIMIT}: {len(index_two)} "
          f"({', '.join(map(str, index_two))})")
    for n in (17, 41, 97, 3343, 9931):
        length = Length(n)
        rows = length.invariant_rows()
        verify_prime_walk(n, length, rows)
        print(f"Prime {n}: index {length.s + 1}, deleted secondary rows "
              f"{[i for i, (keep, _) in rows.items() if not keep]}, "
              f"strict={n in (17, 41)}")
        if n in (17, 41):
            assert length.s == 1 and all(keep and row[0] > 0 for keep, row in rows.values())
        else:
            assert n in classified


def verify_short(certificates):
    records = {int(n): record for n, record in certificates.items()}
    assert set(records) == set(range(7, SHORT_LIMIT + 1, 2)) | set(EXAMPLES)
    print("n    r   k-r   c   d  | vanishes | facets | apex-only | other nonfacets")
    for n in sorted(records):
        length = Length(n)
        implicit, span = cone_certificate(length, records[n])
        c = length.q - len(span)
        if n in EXAMPLES:
            assert c == EXAMPLES[n]
        if n <= SHORT_LIMIT and c:
            results = verify_faces(length, implicit, span, records[n]["faces"])
            print(f"{n:<4}{length.r:>3}{length.q:>6}{c:>4}{len(span):>4}  | "
                  + " | ".join(str(list(value)) for value in
                               (implicit, results["facet"], results["apex"], results["nonfacet"])))
        elif n in EXAMPLES:
            print(f"{n:<4}{length.r:>3}{length.q:>6}{c:>4}{len(span):>4}  | {list(implicit)}")
        if n == 277:
            # Even all orbit-invariant Gordan relations miss one pin here.
            rows = length.invariant_rows()
            relations = kernel(zip(*(rows[i][1] for i in implicit)), len(implicit))
            rotated = dict(length.rotated_rows(implicit))
            combinations = (primitive(relation_value([rotated[t, i] for i in implicit], lam))
                            for t in range(length.r) for lam in relations)
            assert len(echelon(combinations)) == 90
    print("At 277 the rotated orbit-invariant Gordan sums have rank 90 over Q.")


if __name__ == "__main__":
    verify_short(json.loads((HERE / "short_certificates.json").read_text()))
    verify_census([json.loads(line) for line in (HERE / "census_certificates.jsonl").read_text().splitlines()])
