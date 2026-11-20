"""Reconstruct leading maps, phase normals and covector certificates."""
import argparse
from collections import Counter, deque
from fractions import Fraction as Q
from gzip import open as gzip_open
from itertools import combinations, permutations, product
import json
from math import gcd, lcm
from pathlib import Path


N = 5
ZERO = (0,) * N
I = tuple(tuple(int(i == j) for j in range(N)) for i in range(N))
CENTER = tuple(Q(5-i, 6) for i in range(N))


def dot(v, w):
    return sum(a*b for a, b in zip(v, w))


def difference(v, w):
    return tuple(a-b for a, b in zip(v, w))


def multiply(a, b):
    bt = tuple(zip(*b))
    return tuple(tuple(dot(row, col) for col in bt) for row in a)


def row_times(v, a):
    return tuple(dot(v, col) for col in zip(*a))


def echelon(rows, n=N):
    rows = [list(map(Q, row)) for row in rows if any(row)]
    pivots = []
    for j in range(n):
        k = len(pivots)
        pivot = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        rows[k] = [v/rows[k][j] for v in rows[k]]
        for i, row in enumerate(rows):
            if i != k and row[j]:
                rows[i] = [v-row[j]*w for v, w in zip(row, rows[k])]
        pivots.append(j)
    return tuple(tuple(row) for row in rows[:len(pivots)]), tuple(pivots)


def normal_basis(rows):
    reduced, _ = echelon(rows)
    return reduced


def null_basis(rows):
    reduced, pivots = echelon(rows)
    results = []
    for j in range(N):
        if j not in pivots:
            v = list(I[j])
            for row, p in zip(reduced, pivots):
                v[p] = -row[j]
            results.append(tuple(v))
    return tuple(results)


def primitive(v):
    divisor = gcd(*v)
    if next((a for a in v if a), 0) < 0:
        divisor = -divisor
    assert divisor != 0
    return tuple(x//divisor for x in v)


def integer_row(v):
    scale = lcm(*(a.denominator for a in v))
    return primitive(tuple(int(x*scale) for x in v))


def leading(q):
    sorted_digits = sorted((1, *q, *(1-x for x in q)), reverse=True)
    return tuple(sorted_digits[j]-sorted_digits[-1-j] for j in range(N))


def pairing(q):
    levels = tuple((q[j], 1-q[4-j]) for j in range(N))
    return (q[-1] > 0 and max(levels[0]) < 1
            and all(min(a) > max(b) for a, b in zip(levels, levels[1:])))


def signed_groups(q):
    by_value = {}
    for j, value in enumerate(q):
        by_value.setdefault(abs(value-Q(1, 2)), []).append(j)
    return tuple(tuple(by_value[k]) for k in sorted(by_value, reverse=True))


def jacobian(selected):
    # Full leading sort: positive complement on top, negative on bottom.
    upper = [tuple(s*v for v in I[j]) for j, s in selected]
    return tuple(tuple(a+b for a,b in zip(v,w)) for v,w in zip([ZERO]+upper,upper))


def signature(q):
    groups = signed_groups(q)
    sign = tuple((x > Q(1, 2))-(x < Q(1, 2)) for x in q)
    counts = [[sum(sign[j] == 1 for j in group), sum(sign[j] == -1 for j in group)]
              for group in groups if sign[group[0]]]
    return [sign.count(0), counts, any(x in (0, 1) for x in q)]


def phase_data(q, expected):
    groups = signed_groups(q)
    sign = tuple((x > Q(1, 2))-(x < Q(1, 2)) for x in q)
    assert signature(q) == expected
    ties = []
    for group in groups:
        j = group[0]
        if sign[j]:
            for i in group[1:]:
                ties.append(tuple(sign[i]*a-sign[j]*b for a,b in zip(I[i], I[j])))
        else:
            ties.extend(I[i] for i in group)
    choices = []
    for orders in product(*(permutations(group) for group in groups)):
        order = sum(orders, ())
        for signs in product(*((sign[j],) if sign[j] else (1,-1) for j in order)):
            choices.append(tuple(zip(order, signs)))
    return tuple(jacobian(choice) for choice in choices), tuple(ties)


def family(row):
    witness = q = tuple(map(Q, row["witness"]))
    phases, equations = [], []
    chart = I
    points = []
    for face in row["faces"]:
        points.append(q)
        choices, ties = phase_data(q, face)
        phases.append(choices)
        equations.extend(row_times(c, chart) for c in ties)
        chart = multiply(choices[0], chart)
        q = leading(q)
    assert q == witness
    assert row["paired"] == all(map(pairing, points))
    period = row["period"]
    assert len(points) % period == 0 and len(set(points)) == period
    assert all(q == points[i % period] for i,q in enumerate(points))
    equations.extend(difference(a,b) for a,b in zip(chart,I))
    basis = null_basis(equations)
    assert len(basis) == row["dimension"]
    return phases, basis, equations, tuple(points)


def find_orbit(phases, seed, coefficient_cap=4, size_cap=64):
    path = {(0, seed)}
    queue = deque(path)
    while queue:
        step, v = queue.popleft()
        prev = (step-1) % len(phases)
        for matrix in phases[prev]:
            w = row_times(v,matrix)
            assert gcd(*w) == 1
            state = prev, primitive(w)
            if state not in path:
                if max(map(abs,w)) > coefficient_cap or len(path) == size_cap:
                    return None
                path.add(state)
                queue.append(state)
    return path


def verify_orbit(phases, tangents, seed, coefficient_cap=4, size_cap=64):
    assert all(dot(seed,t) == 0 for t in tangents)
    orbit = find_orbit(phases, seed, coefficient_cap, size_cap)
    if orbit is not None:
        assert all(dot(v,t) == 0 for p,v in orbit if p == 0 for t in tangents)
    return orbit


def seeds(tangents, equations):
    direct = {primitive(v) for v in equations if any(v)}
    for v in product((0,1,-1), repeat=N):
        if next((a for a in v if a), 0) == 1 and all(dot(v,t) == 0 for t in tangents):
            direct.add(v)
    return sorted(direct, key=lambda v:(sum(map(abs,v)),v))


def inverse(a):
    rows = [tuple(x)+y for x,y in zip(a,I)]
    reduced, pivots = echelon(rows)
    assert len(pivots) == N
    result = tuple(tuple(int(v) for v in row[N:]) for row in reduced)
    assert multiply(a,result) == I
    return result


def power(a,exponent):
    result = I
    while exponent:
        if exponent % 2:
            result = multiply(result,a)
        exponent //= 2
        if exponent:
            a = multiply(a,a)
    return result


def generators(phases):
    """The base return and single adjacent-choice deviations in its frame."""
    prefix = I
    variations = set()
    for options in phases:
        branch = options[0]
        inv = inverse(multiply(branch,prefix)) if len(options) > 1 else None
        for option in options[1:]:
            variations.add(multiply(inv, multiply(option,prefix)))
        prefix = multiply(branch,prefix)
    return prefix, tuple(sorted(variations))


def orbit_subspace(phases, tangents):
    """Linearly constrain all phase-zero normal rows with finite orbit."""
    direct, alternatives = generators(phases)
    full = tuple(sorted(set((direct,*alternatives))-{I}))
    invariant = normal_basis(tangents)

    # Root-of-unity eigenvalues of an integral 5x5 word have algebraic
    # degree at most five. A periodic row therefore has period dividing 120.
    seen = {I}
    for word in (*full, *(multiply(a,b) for a,b in combinations(full,2))):
        if word in seen:
            continue
        seen.add(word)
        displaced = tuple(difference(a,b) for a,b in zip(power(word,120),I))
        invariant = normal_basis((*invariant,*zip(*displaced)))
        if len(invariant) == N:
            return ()
    while True:
        expanded = normal_basis((*invariant,*(tuple(dot(row,v) for row in a)
                                               for a in full for v in invariant)))
        if len(expanded) == len(invariant):
            return null_basis(expanded)
        invariant = expanded


def contracted_pairing(points):
    contracted = tuple(tuple(c+(v-c)/8 for c,v in zip(CENTER,q)) for q in points)
    for i,q in enumerate(contracted):
        assert pairing(q)
        endpoints = (1,*q,0)
        assert min(a-b for a,b in zip(endpoints,endpoints[1:])) >= Q(7,48)
        assert leading(q) == contracted[(i+1) % len(points)]


def archived_cells(records, energy):
    with gzip_open(energy,"rt") as source:
        high = [row for row in json.load(source)
                if row["location"] == "interior" and row["meets_high_energy"]]
    assert len(high) == 5
    overlap = [i for i,row in enumerate(records) if not row["paired"]
               and any(row["faces"] == h["faces"] for h in high)]
    assert overlap == [176,177,179]
    return overlap


def denominator_25(records):
    known = tuple(Q(x,25) for x in (22,19,16,10,6))
    cycle = (known,leading(known),leading(leading(known)))
    word = [signature(q) for q in cycle]
    assert leading(cycle[-1]) == cycle[0]
    matches = {}
    for i,row in enumerate(records):
        faces = row["faces"]
        for j in range(row["period"]):
            shifted = [faces[(j+k) % len(faces)] for k in range(len(faces))]
            if len(faces) % 3 == 0 and shifted == word * (len(faces)//3):
                matches[i] = j, shifted
                break
    assert set(matches) == {147}
    j, shifted = matches[147]
    q = tuple(map(Q, records[147]["witness"]))
    for _ in range(j):
        q = leading(q)
    _, tangents, _, points = family(dict(records[147], witness=q, faces=shifted))
    normals = ((0,1,0,0,1),(1,0,-1,0,-1),(2,0,0,0,1))
    assert len(tangents) == 2 and all(dot(n,t) == 0 for n in normals for t in tangents)
    assert tuple(dot(n,points[0]) for n in normals) == (1,0,2)
    assert tuple(dot(n,known) for n in normals) == (1,0,2)
    return 147


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    parser.add_argument("--catalogue", type=Path, default=here.parent/"entry_pi/catalogue.json")
    parser.add_argument("--energy", type=Path, default=here.parent/"global_work/evidence/energy.json.gz")
    parser.add_argument("--astra", type=Path, required=True)
    parser.add_argument("--diagnose", action="store_true")
    args = parser.parse_args()
    import ast
    import re

    reference = {}
    format_ = r"^(\d+) dimension (\d+) period (\d+) row (\([^)]+\)) orbit size up to sign (\d+)$"
    for line in args.astra.read_text().splitlines():
        match = re.match(format_,line)
        if match:
            i, d, p, v, size = match.groups()
            reference[int(i)] = int(d), int(p), ast.literal_eval(v), int(size)
    assert len(reference) == 54
    records = json.loads(args.catalogue.read_text())
    passed, remaining, paired = {}, {}, []
    for i,row in enumerate(records):
        phases, tangents, equations, points = family(row)
        if row["paired"]:
            contracted_pairing(points)
            paired.append(i)
            continue
        if i in reference:
            d,p,v,size = reference[i]
            assert (row["dimension"], row["period"]) == (d,p)
            orbit = verify_orbit(phases,tangents,v)
            assert orbit is not None and len(orbit) == size
        for seed in seeds(tangents,equations):
            orbit = verify_orbit(phases,tangents,seed)
            if orbit is not None:
                passed[i] = (row["dimension"], seed, len(orbit))
                break
        else:
            remaining[i] = (phases,tangents,equations)

    assert set(reference) == set(passed)
    print("Catalogue cells:", len(records))
    print("Wholly strictly paired (all contracted itineraries verified):", len(paired))
    print("Astra certificates independently verified:", len(reference))
    print("Independent strict-budget search:", len(passed),
          dict(sorted(Counter(d for d,_,_ in passed.values()).items())))
    print("Extra rows in strict-budget search:", sorted(set(passed)-set(reference)))
    print("Uncertified rows by dimension and period:",
          dict(sorted(Counter((records[i]["dimension"],records[i]["period"])
                              for i in remaining).items())))
    if args.diagnose:
        found = {}
        unresolved = []
        zero_dim = set()
        for i,(phases,tangents,_) in remaining.items():
            candidate_space = orbit_subspace(phases,tangents)
            if not candidate_space:
                zero_dim.add(i)
                continue
            # If basis rows close, all rows in the resulting space have finite orbit.
            for vector in candidate_space:
                seed = integer_row(vector)
                orbit = verify_orbit(phases,tangents,seed,coefficient_cap=64,size_cap=4096)
                if orbit is not None:
                    found[i] = (seed,len(orbit),len(candidate_space))
                    break
            else:
                unresolved.append((i,len(candidate_space)))
        print("Additional closed orbits:", found)
        print("Other possible finite-orbit rows:", unresolved)
        print("Zero upper dimension for phase-zero normal finite orbits:", len(zero_dim))
        assert len(found) == 17 and not unresolved and len(zero_dim) == 90
        assert all(dim == 1 for _,_,dim in found.values())
        triangle = archived_cells(records,args.energy)
        patch = denominator_25(records)
        assert set(triangle) | {patch} <= zero_dim
        print("Already-null Lemma M triangle cells:", triangle)
        print("Denominator-25 whole-cell extension:", patch)
        print("Residual unpaired rows:", len(set(remaining)-set(found)-set(triangle)-{patch}))


if __name__ == "__main__":
    main()
