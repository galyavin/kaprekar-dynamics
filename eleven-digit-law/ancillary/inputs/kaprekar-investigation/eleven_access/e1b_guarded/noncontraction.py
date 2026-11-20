"""Exact finite-time transverse noncontraction, supplementary to covectors."""

from guarded import I, N, difference, dot, multiply, negate, normal_basis
from guarded import pull, row_times, signed_primitive


def verify_norm(rows, tangents):
    assert all(dot(r, t) == 0 for r in rows for t in tangents)
    assert len(normal_basis(rows)) + len(tangents) == N


def survivors(cones, phases, start, norm, rounds):
    """Test all strict words with every later norm smaller than the initial one."""
    signed_rows = sorted(set((*norm, *(negate(r) for r in norm))))
    counts = [0] * (rounds + 1)
    for maximum in signed_rows:
        inequalities = tuple(signed_primitive(difference(maximum, r))
                             for r in signed_rows if r != maximum)
        if not cones.feasible(inequalities):
            continue
        stack = [(0, I, inequalities)]
        while stack:
            time, A, guards = stack.pop()
            if time % len(phases) == 0:
                counts[time // len(phases)] += 1
            if time == rounds * len(phases):
                continue
            phase = (start + time) % len(phases)
            for B, local in phases[phase]:
                conditions = (*guards, *pull(local, A))
                if not cones.feasible(conditions):
                    continue
                W = multiply(B, A)
                if (time + 1) % len(phases) == 0:
                    smaller = tuple(difference(maximum, row_times(r, W)) for r in signed_rows)
                    if not all(map(any, smaller)):
                        continue
                    conditions = tuple(sorted(set((*conditions, *map(signed_primitive, smaller)))))
                    if not cones.feasible(conditions):
                        continue
                stack.append((time + 1, W, conditions))
    return counts

