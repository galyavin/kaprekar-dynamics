"""Anchor and gap-majorization certificates for isolated periodic orbits."""

from guarded import I, N, difference, dot, multiply, negate, pull, signed_primitive


GAPS = (negate(I[0]), *(difference(a, b) for a, b in zip(I, I[1:])), I[-1])


def gap_values(q):
    endpoints = (1, *q, 0)
    return tuple(a - b for a, b in zip(endpoints, endpoints[1:]))


def majorization_guards(q):
    """Strict inequalities at group boundaries, after omitting their null hyperplanes."""
    g = gap_values(q)
    guards = [negate(I[0])]
    total = (0,) * N
    for value in sorted(set(g))[:-1]:
        total = tuple(a + sum(GAPS[j][i] for j, v in enumerate(g) if v == value)
                      for i, a in enumerate(total))
        guards.append(signed_primitive(negate(total)))
    return tuple(guards)


def gap_order_profiles(tangents, points):
    """All gap-order profiles of a one-dimensional relative-interior sorting cell."""
    from fractions import Fraction as Q

    triples = []
    breakpoints = set()
    for q, (v,) in zip(points, tangents):
        gaps = gap_values(q)
        slopes = tuple(dot(row, v) for row in GAPS)
        assert min(gaps) > 0
        triples.extend(zip(gaps, slopes))

        # All magnitudes are signed on strict nonzero levels; tangent ties persist.
        magnitudes = []
        for a, b in zip(q, v):
            a -= Q(1, 2)
            if a:
                signed = (abs(a), b if a > 0 else -b)
                triples.append(signed)
                magnitudes.append(signed)
            else:
                assert b == 0
                magnitudes.append((Q(0), Q(0)))
        magnitudes.sort(reverse=True)
        for (a, b), (c, d) in zip(magnitudes, magnitudes[1:]):
            if a != c:
                triples.append((a - c, b - d))
            else:
                assert b == d

        for j, a in enumerate(gaps):
            for k in range(j):
                b = slopes[j] - slopes[k]
                if b:
                    breakpoints.add((gaps[k] - a) / b)

    lo = max(-a / b for a, b in triples if b > 0)
    hi = min(-a / b for a, b in triples if b < 0)
    assert lo < 0 < hi
    critical = [lo, *sorted(t for t in breakpoints if lo < t < hi), hi]
    representatives = [*critical[1:-1], *((a+b)/2 for a, b in zip(critical, critical[1:]))]
    return tuple((t, tuple(tuple(a+t*b for a, b in zip(q, v))
                           for q, (v,) in zip(points, tangents)))
                 for t in sorted(representatives))


def constrained_counts(cones, phases, start, guards_by_phase, length, cap=8000):
    """Exact finite words with joint strict sorting and mixed limit constraints."""
    strict, closed = guards_by_phase[start]
    counts = [0] * (length + 1)
    stack = [(0, I, strict, closed)]
    while stack:
        step, matrix, strict, closed = stack.pop()
        if not cones.check(strict, closed=closed):
            continue
        counts[step] += 1
        if counts[step] == cap:
            return counts, False
        if step == length:
            continue
        phase = (start + step) % len(phases)
        next_strict, next_closed = guards_by_phase[(phase+1) % len(phases)]
        for A, local in phases[phase]:
            next_matrix = multiply(A, matrix)
            strict_guards = (*strict, *pull(local, matrix), *pull(next_strict, next_matrix))
            closed_guards = (*closed, *pull(next_closed, next_matrix))
            stack.append((step+1, next_matrix, tuple(sorted(set(strict_guards))),
                          tuple(sorted(set(closed_guards)))))
    return counts, True
