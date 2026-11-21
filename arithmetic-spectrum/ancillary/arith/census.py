"""Generate and replay exact P10 chamber certificates for the frozen Task 207 panel."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from math import gcd, isqrt, lcm, prod
import json
from pathlib import Path
import time

MODULUS = 65521  # prime, used only to certify full rank over Q
BANDS = ((10_001, 20_001), (90_001, 100_001))
BINS = ("P1", "P2", "P3-5", "P6+", "C2-3", "C4-12", "C13+/2", "C13+/3+")
CLASSES = ("strict", "scarred", "full-orbit", "zero", "unknown")


@cache
def factors(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            result.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        result.append((n, 1))
    return tuple(result)


@cache
def layer(d):
    phi = prod((p - 1) * p ** (e - 1) for p, e in factors(d))
    order = phi
    for p, _ in factors(phi):
        while order % p == 0 and pow(2, order // p, d) == 1:
            order //= p
    if order % 2 == 0 and pow(2, order // 2, d) == d - 1:
        order //= 2
    return order, phi // (2 * order)


def orbit_count(n):
    divisors = [1]
    for p, e in factors(n):
        divisors = [d * p ** a for d in divisors for a in range(e + 1)]
    return sum(layer(d)[1] for d in divisors if d > 1)


def sample_bin(n, t):
    omega = sum(e for _, e in factors(n))
    if omega == 1:
        if t <= 2:
            return f"P{t}"
        return "P3-5" if t <= 5 else "P6+"
    if t <= 3:
        return "C2-3"
    if t <= 12:
        return "C4-12"
    return "C13+/2" if omega == 2 else "C13+/3+"


def panel():
    for n in range(7, 10_000, 2):
        t = orbit_count(n)
        yield n, t, "census", sample_bin(n, t)
    for start, stop in BANDS:
        counts = Counter()
        for n in range(start, stop, 2):
            t = orbit_count(n)
            kind = sample_bin(n, t)
            if counts[kind] == 4:
                continue
            counts[kind] += 1
            yield n, t, str(start), kind
            if sum(counts.values()) == 4 * len(BINS):
                break


def orbit_labels(n):
    """Number positive orbits in order of their least slot; 0 denotes J."""
    k = (n - 1) // 2
    labels = [0] + [-1] * k
    t = 0
    for j in range(1, k + 1):
        if labels[j] != -1:
            continue
        slot = j
        while labels[slot] == -1:
            labels[slot] = t
            slot = min(2 * slot, n - 2 * slot)
        if slot != j:
            raise ValueError(f"folded successor is not a permutation at {n}")
        t += 1
    return labels, t


def profile_rows(labels, offsets):
    """Exact literal comparisons with the value on the root set shifted to 0."""
    profile = [0, *offsets]
    k = len(labels) - 1
    delta = 0
    for i in range(1, (k + 1) // 2 + 1):
        delta += profile[labels[k + 1 - i]] - profile[labels[i - 1]]
        if i >= 2:
            yield delta


def selected_rows(labels, indices, dim):
    wanted = set(indices)
    k = len(labels) - 1
    current = [0] * (dim + 1)
    answer = {}
    for i in range(1, max(wanted) + 1):
        current[labels[k + 1 - i]] += 1
        current[labels[i - 1]] -= 1
        if i in wanted:
            answer[i] = current[1:]
    return answer


def rank_indices(rows):
    """Return the original indices of an independent set modulo MODULUS."""
    matrix = [[int(value) % MODULUS for value in row] for row in rows]
    if not matrix:
        return []
    indices = list(range(len(matrix)))
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        indices[rank], indices[pivot] = indices[pivot], indices[rank]
        pivot_row = matrix[rank]
        inverse = pow(pivot_row[column], -1, MODULUS)
        for row in matrix[rank + 1:]:
            multiple = row[column] * inverse % MODULUS
            if multiple:
                for j in range(column, len(row)):
                    row[j] = (row[j] - multiple * pivot_row[j]) % MODULUS
        rank += 1
        if rank == len(matrix):
            break
    return indices[:rank]

def validate(row):
    n, t, cert = row["n"], row["t"], row["cert"]
    cls = row["cls"]
    labels, exact_t = orbit_labels(n)
    assert t == exact_t == orbit_count(n), n
    if cls == "unknown":
        assert isinstance(cert.get("reason"), str) and cert["reason"], n
        return
    if cls == "full-orbit":
        assert t == 1, n
        return
    assert t > 1 and cls in CLASSES, n
    if cls in ("strict", "scarred"):
        x = cert["x"]
        assert len(x) == t - 1 and all(type(a) is int for a in x) and any(x), n
        assert min(profile_rows(labels, x)) >= (1 if cls == "strict" else 0), n
    if cls == "strict":
        return

    if "tie" in cert:
        i = cert["tie"]
        assert cls == "scarred" and 2 <= i <= len(labels) // 2, n
        assert not any(selected_rows(labels, [i], t - 1)[i]), (n, i)
        return

    terms = cert["relation"]
    assert terms and len({i for i, _ in terms}) == len(terms), n
    assert all(2 <= i <= len(labels) // 2 and type(c) is int and c > 0 for i, c in terms), n
    exact = selected_rows(labels, [i for i, _ in terms], t - 1)
    assert all(sum(c * exact[i][j] for i, c in terms) == 0 for j in range(t - 1)), n
    if cls == "zero":
        assert len(rank_indices([exact[i] for i, _ in terms])) == t - 1, n


def index_two_identity(row):
    """Compare literal rows with Legendre sums and the reduced-form class number."""
    p = row["n"]
    k = (p - 1) // 2
    assert p % 8 == 1 and row["r"] == k // 2, p
    squares = {a * a % p for a in range(1, k + 1)}
    partial = [0]
    for a in range(1, k + 1):
        partial.append(partial[-1] + (1 if a in squares else -1))
    labels, t = orbit_labels(p)
    assert t == 2 and all((j in squares) == (labels[j] == 0) for j in range(1, k + 1)), p
    defects = [partial[2 * i - 1] - 2 * partial[i - 1] - 1 for i in range(2, k // 2 + 1)]
    assert [2 * d for d in profile_rows(labels, [-(k + 1)])] == [
        (k + 1) * d for d in defects
    ], p
    discriminant = -4 * p
    h = 0
    # Enumerate reduced primitive positive definite forms [a,b,c].
    for a in range(1, isqrt((-discriminant) // 3) + 1):
        for b in range(-a + a % 2, a + 1, 2):
            c, rem = divmod(b * b - discriminant, 4 * a)
            if rem or c < a or (b < 0 and (-b == a or a == c)):
                continue
            if gcd(a, b, c) == 1:
                h += 1
    assert defects[-1] == -h < 0, p
    largest = max(defects)
    expected = "strict" if largest < 0 else "scarred" if largest == 0 else "zero"
    assert row["cls"] in (expected, "unknown"), p
    return len(defects)


def low_orbit_check(row):
    n = row["n"]
    expected = "strict" if n in (9, 15, 21) else "scarred"
    assert row["cls"] in (expected, "unknown"), n
    if expected == "scarred":
        i = 2 if n % 3 else 4
        labels, t = orbit_labels(n)
        assert not any(selected_rows(labels, [i], t - 1)[i]), (n, i)


def walk_matrix(labels, dim):
    """Find candidate witnesses using alternating folded-prefix counts."""
    import numpy as np

    stop = len(labels) // 2
    matrix = np.empty((stop - 1, dim), dtype=np.int32)
    counts = np.zeros(dim + 1, dtype=np.int32)
    for i in range(2, stop + 1):
        counts[labels[2 * i - 2]] -= 1
        counts[labels[2 * i - 1]] += 1
        matrix[i - 2] = counts[1:]
    return matrix


def solve_lp(deadline, objective, **kwargs):
    from scipy.optimize import linprog

    ttl = min(2.0, deadline - time.monotonic())
    if ttl <= 0:
        return None
    result = linprog(objective, method="highs", options={"time_limit": ttl}, **kwargs)
    return result.x if result.success else None


def integer_profile(candidate, labels, strict):
    if candidate is None:
        return None
    for bits in (8, 16, 24, 32, 40):
        x = [round(float(a) * (1 << bits)) for a in candidate]
        common = gcd(*x)
        if common:
            x = [a // common for a in x]
            if min(profile_rows(labels, x)) >= int(strict):
                return x
    return None


def solve_exact(columns, target):
    """Solve a rational linear system with independent supplied columns."""
    count = len(columns)
    matrix = [
        [Fraction(int(col[i])) for col in columns] + [Fraction(int(rhs))]
        for i, rhs in enumerate(target)
    ]
    for j in range(count):
        pivot = next((i for i in range(j, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            return None
        matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
        row = matrix[j]
        matrix[j] = row = [a / row[j] for a in row]
        for i, other in enumerate(matrix):
            if i != j and other[j]:
                f = other[j]
                matrix[i] = [a - f * b for a, b in zip(other, row)]
    if any(row[-1] for row in matrix[count:]):
        return None
    return [row[-1] for row in matrix[:count]]


def primal_relation(rows, reference, deadline):
    """Find nonnegative weights on primitive rows cancelling reference rows."""
    import numpy as np

    count, dim = rows.shape
    if reference:
        system = rows.T
        target = -rows[reference].sum(axis=0)
    else:
        system = np.vstack((rows.T, np.ones(count, dtype=np.int64)))
        target = [0] * dim + [1]
    candidate = solve_lp(
        deadline, np.ones(count), A_eq=system, b_eq=target, bounds=(0, None)
    )
    if candidate is None:
        return None
    support = [int(i) for i in np.flatnonzero(candidate > 1e-10)]
    weights = solve_exact([system[:, i] for i in support], target)
    if weights is None or any(a < 0 for a in weights):
        return None
    relation = {i: a for i, a in zip(support, weights) if a}
    for i in reference:
        relation[i] = relation.get(i, Fraction()) + 1
    return relation


def relation_certificate(relation, representatives, divisors):
    items = [
        (int(representatives[i]) + 2, a / int(divisors[representatives[i]]))
        for i, a in sorted(relation.items())
    ]
    multiple = lcm(*(a.denominator for _, a in items))
    return [[i, int(a * multiple)] for i, a in items]


def classify(n, t, deadline):
    if t == 1:
        return "full-orbit", {}
    if time.monotonic() >= deadline:
        return "unknown", {"reason": "campaign deadline"}
    labels, exact_t = orbit_labels(n)
    assert t == exact_t, n
    dim = t - 1
    if min(profile_rows(labels, [1] * dim)) > 0:
        return "strict", {"x": [1] * dim}

    import numpy as np

    weak = None
    if len(factors(n)) > 1 or factors(n)[0][1] > 1:
        representatives = [labels.index(i) for i in range(1, t)]
        weak = [int(j % factors(n)[0][0] == 0) for j in representatives]
    else:
        for x in ([1] * dim, [-1] * dim):
            if min(profile_rows(labels, x)) >= 0:
                weak = x
                break
    matrix = walk_matrix(labels, dim)
    divisors = np.gcd.reduce(matrix, axis=1)
    ties = np.flatnonzero(divisors == 0)
    obstruction = {"tie": int(ties[0]) + 2} if len(ties) else None
    if obstruction and weak is not None:
        return "scarred", {"x": weak, **obstruction}

    active = np.flatnonzero(divisors)
    if not len(active):
        return "scarred", {"x": [1] * dim, **obstruction}
    primitive, positions = np.unique(
        matrix[active] // divisors[active, None], axis=0, return_index=True
    )
    representatives = active[positions]
    count = len(primitive)

    if obstruction is None:
        by_row = {tuple(int(a) for a in row): i for i, row in enumerate(primitive)}
        for row, i in by_row.items():
            j = by_row.get(tuple(-a for a in row))
            if j is not None:
                obstruction = {"relation": relation_certificate(
                    {i: Fraction(1), j: Fraction(1)}, representatives, divisors
                )}
                break

    if obstruction is None:
        candidate = solve_lp(
            deadline, np.zeros(dim), A_ub=-primitive, b_ub=-np.ones(count),
            bounds=(None, None)
        )
        x = integer_profile(candidate, labels, strict=True)
        if x is not None:
            return "strict", {"x": x}
        relation = primal_relation(primitive, [], deadline)
        if relation:
            obstruction = {"relation": relation_certificate(relation, representatives, divisors)}
    if obstruction and weak is not None:
        return "scarred", {"x": weak, **obstruction}

    if weak is None:
        candidate = solve_lp(
            deadline, -primitive.sum(axis=0), A_ub=-primitive,
            b_ub=np.zeros(count), bounds=(-1, 1)
        )
        weak = integer_profile(candidate, labels, strict=False)
        if weak is not None and obstruction is not None:
            return "scarred", {"x": weak, **obstruction}

    reference = rank_indices(primitive)
    if len(reference) == dim:
        relation = primal_relation(primitive, reference, deadline)
        if relation:
            return "zero", {"relation": relation_certificate(relation, representatives, divisors)}
    return "unknown", {"reason": "no exact classification within candidate limits"}


def generate(path, seconds):
    deadline = time.monotonic() + seconds
    with path.open("x") as f:
        for n, t, band, kind in panel():
            cls, cert = classify(n, t, deadline)
            row = {"n": n, "t": t, "r": layer(n)[0], "band": band, "bin": kind,
                   "cls": cls, "cert": cert}
            try:
                validate(row)
            except AssertionError:
                row["cls"] = "unknown"
                row["cert"] = {"reason": "candidate failed exact validation"}
            f.write(json.dumps(row, separators=(",", ":")) + "\n")


def replay(path, overlap_path):
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    identities = [(r["n"], r["t"], r["band"], r["bin"]) for r in rows]
    assert identities == list(panel()), "panel changed or census incomplete"
    for row in rows:
        assert row["r"] == layer(row["n"])[0], row["n"]
        validate(row)
    by_n = {r["n"]: r for r in rows}
    overlap = json.loads(overlap_path.read_text())
    for earlier in overlap:
        row = by_n[earlier["n"]]
        assert row["r"] == earlier["r"] and row["t"] == earlier["s"] + 1, row["n"]
        assert row["cls"] in (earlier["cls"], "unknown"), row["n"]
    print("verified", len(rows), "panel rows;", len(overlap), "overlap rows")
    index_two = [r for r in rows if r["bin"] == "P2"]
    print("index two:", len(index_two), "class numbers;",
          sum(index_two_identity(r) for r in index_two), "secondary rows")
    low_orbit = [r for r in rows if r["bin"] == "C2-3"]
    for row in low_orbit:
        low_orbit_check(row)
    print("low-orbit composites:", len(low_orbit), "checked;",
          "strict lengths", [r["n"] for r in low_orbit if r["cls"] == "strict"])

    def print_counts(label, selected):
        counts = Counter(r["cls"] for r in selected)
        ties = sum(r["cls"] == "scarred" and "tie" in r["cert"] for r in selected)
        print(label, *(counts[c] for c in CLASSES), ties, counts["scarred"] - ties)

    for label, select in (
        ("census primes", lambda r: r["band"] == "census" and r["bin"][0] == "P"),
        ("census composites", lambda r: r["band"] == "census" and r["bin"][0] == "C"),
        ("census 7–1599", lambda r: r["band"] == "census" and r["n"] <= 1599),
        ("census 1601–4999", lambda r: r["band"] == "census" and 1599 < r["n"] <= 4999),
        ("census 5001–9999", lambda r: r["band"] == "census" and r["n"] >= 5001),
    ):
        print_counts(label, [r for r in rows if select(r)])
    for band in ("census", *(str(start) for start, _ in BANDS)):
        for kind in BINS:
            print_counts(f"{band} {kind}", [r for r in rows if r["band"] == band and r["bin"] == kind])
    for n in (2, 3, 4, 5):
        for category in ("P", "C"):
            print_counts(f"census exact {category} t={n}", [
                r for r in rows if r["band"] == "census" and r["bin"][0] == category and r["t"] == n
            ])
    print("unknown", [(r["n"], r["cert"]["reason"]) for r in rows if r["cls"] == "unknown"])


def planar_scarred_witness(labels):
    matrix = walk_matrix(labels, 2)
    ties = [i for i, row in enumerate(matrix, 2) if not any(row)]
    if not ties:
        return None
    for row in matrix:
        a, b = map(int, row)
        common = gcd(a, b)
        if not common:
            continue
        for x in ([-b // common, a // common], [b // common, -a // common]):
            if min(profile_rows(labels, x)) >= 0:
                return {"x": x, "tie": ties[0]}
    return None


def recover_planar(path, output):
    """Certify unresolved planar cases without replacing the initial census."""
    with output.open("x") as f:
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if row["cls"] != "unknown" or row["bin"][0] != "P" or row["t"] != 3:
                continue
            labels, _ = orbit_labels(row["n"])
            cert = planar_scarred_witness(labels)
            if cert is None:
                continue
            row["cls"], row["cert"] = "scarred", cert
            validate(row)
            f.write(json.dumps(row, separators=(",", ":")) + "\n")
            print(row["n"], row["cls"], cert)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    counter = subparsers.add_parser("generate")
    counter.add_argument("path", type=Path)
    counter.add_argument("--seconds", type=float, default=2_400)
    checker = subparsers.add_parser("replay")
    checker.add_argument("path", type=Path)
    checker.add_argument("overlap", type=Path)
    rescue = subparsers.add_parser("recover-planar")
    rescue.add_argument("path", type=Path)
    rescue.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.mode == "generate":
        generate(args.path, args.seconds)
    elif args.mode == "replay":
        print("columns: strict scarred full-orbit zero unknown | scarred-tie scarred-relation")
        replay(args.path, args.overlap)
    else:
        recover_planar(args.path, args.output)


if __name__ == "__main__":
    main()
