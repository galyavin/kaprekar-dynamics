"""Finite covector checks with jointly feasible phase-locked words."""

from collections import defaultdict, deque
from itertools import combinations
from math import gcd

from guarded import I, N, difference, dot, integer_row, multiply, normal_basis, null_basis
from guarded import pull, row_times
from catalogue_checks import power, primitive, seeds


def add_fixed_constraints(columns, W):
    difference_from_identity = tuple(difference(a, b) for a, b in zip(power(W, 120), I))
    return normal_basis((*columns, *zip(*difference_from_identity)))


def monoid_space(matrices, tangent):
    """Necessary candidate space under freely concatenating these matrices."""
    matrices = sorted(set(matrices) - {I})
    columns = normal_basis(tangent)
    for W in matrices:
        columns = add_fixed_constraints(columns, W)
        if len(columns) == N:
            return ()
    for A, B in combinations(matrices, 2):
        columns = add_fixed_constraints(columns, multiply(A, B))
        if len(columns) == N:
            return ()
    while True:
        expanded = normal_basis((*columns, *(tuple(dot(row, v) for row in A)
                                             for A in matrices for v in columns)))
        if len(expanded) == len(columns):
            return null_basis(columns)
        columns = expanded


def guarded_orbit(cones, phases, tangents, phase, seed, row_cap=128, state_cap=2048):
    """Close exact backwards states, subsuming narrower cones with the same row."""
    state = phase, seed, ()
    covered = defaultdict(list, {(phase, seed): [()]})
    queue = deque([state])
    added = 1
    while queue:
        target, row, guard = queue.popleft()
        if guard not in covered[target, row]:
            continue
        prev = (target - 1) % len(phases)
        for A, local in phases[prev]:
            pulled_guard = (*local, *pull(guard, A))
            if not cones.feasible(pulled_guard):
                continue
            image = row_times(row, A)
            assert gcd(*image) == 1
            w = primitive(image)
            assert all(dot(w, t) == 0 for t in tangents[prev])
            if (prev, w) not in covered and len(covered) == row_cap:
                return None, ("rows", len(covered), added)
            prior = covered[prev, w]
            if any(cones.contains(old, pulled_guard) for old in prior):
                continue
            reduced = cones.facets(pulled_guard)
            prior[:] = [old for old in prior if not cones.contains(reduced, old)]
            prior.append(reduced)
            queue.append((prev, w, reduced))
            added += 1
            if added >= state_cap:
                return None, ("cones", len(covered), added)
    cover = [(p, row, guard) for (p, row), guards in sorted(covered.items())
             for guard in sorted(guards)]
    return cover, ("closed", len(covered), added)


def candidate_rows(basis, tangent, guards):
    """Candidate search only. The basis bounds every rational finite orbit."""
    if not basis:
        return ()
    basis = tuple(integer_row(v) for v in basis)
    normal = normal_basis(basis)
    equations = null_basis(normal)
    direct = {*basis, *seeds(equations, guards)}
    direct = {v for v in direct if all(dot(v, t) == 0 for t in (*tangent, *equations))}
    return tuple(sorted(direct, key=lambda v: (sum(map(abs, v)), v)))


def verify_cover(cones, phases, tangents, phase, seed, cover):
    indexed = defaultdict(list)
    for p, row, guard in cover:
        guard = tuple(map(tuple, guard))
        row = tuple(row)
        assert row == primitive(row) and gcd(*row) == 1
        assert all(dot(row, t) == 0 for t in tangents[p])
        assert cones.feasible(guard)
        indexed[p, row].append(guard)
    assert () in indexed[phase, tuple(seed)]
    for (p, row), guards in list(indexed.items()):
        previous = (p - 1) % len(phases)
        for guard in guards:
            for A, local in phases[previous]:
                domain = (*local, *pull(guard, A))
                if cones.feasible(domain):
                    image = row_times(row, A)
                    assert gcd(*image) == 1
                    assert any(cones.contains(outer, domain)
                               for outer in indexed.get((previous, primitive(image)), ()))
    return sorted(row for p, row in indexed if p == phase)
