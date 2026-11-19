"""Independent borrowed-digit orders and exact per-phase work maxima."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product

from exhaust import N, dot, gaps_from_q, lin, q_from_gaps


def physical(q, offsets):
    return all(a != b or u >= v for a, b, u, v in
               zip(q, q[1:], offsets, offsets[1:]))


def borrow_forms(face, offsets):
    # Offsets are doubled. Sort scores by (leading value, doubled offset).
    # The five signed pairs are the coefficients of the original subtraction.
    q = q_from_gaps(face.exemplar)
    borrowed = sorted([(1, -2)]
                      + [(p, d - 2 * (i == 4)) for i, (p, d) in enumerate(zip(q, offsets))]
                      + [(1 - p, -d - 2 * (i != 0)) for i, (p, d) in enumerate(zip(q, offsets))],
                      reverse=True)
    leading_gap_shift = gaps_from_q(tuple(borrowed[i][1] - borrowed[-1 - i][1]
                                          for i in range(5)), total=0)
    source_gap_shift = gaps_from_q(offsets, total=0)
    output_as_source = [0] * N
    for i, label in enumerate(face.permutation):
        output_as_source[label] = leading_gap_shift[i]
    # Gradient factor two cancels our use of doubled residuals.
    return lin(output_as_source, source_gap_shift)


def residual_groups(q):
    groups = {}
    for i, p in enumerate(q):
        groups.setdefault(abs(2 * p - 1), []).append(i)
    return tuple((radius, tuple(indices)) for radius, indices in sorted(groups.items()))


def compressed_central(size):
    for parity in (0, 1):
        limit = 2 * size - parity
        for ascending in combinations_with_replacement(range(-limit, limit + 1, 2), size):
            magnitudes = sorted(set(map(abs, ascending)))
            permitted_first = (1,) if parity else (0, 2)
            if magnitudes[0] not in permitted_first:
                continue
            if any(b - a > 2 for a, b in zip(magnitudes, magnitudes[1:])):
                continue
            yield ascending[::-1]


@lru_cache(maxsize=None)
def work_forms(face):
    q = q_from_gaps(face.exemplar)
    assert face.exemplar[0] > 0 and q[-1] > 0
    baseline = borrow_forms(face, (0,) * 5)
    separate_groups = []
    translation_checks = []
    for radius, indices in residual_groups(q):
        size = len(indices)
        distinct = {}
        if radius == 0:
            patterns = compressed_central(size)
        else:
            sides = tuple(1 if q[i] > Fraction(1, 2) else -1 for i in indices)
            patterns = (tuple(2 * x * side for x, side in zip(levels, sides))
                        for levels in product(range(size), repeat=size)
                        if min(levels) == 0 and max(levels) + 1 == len(set(levels)))
            translation = [0] * 5
            for i, side in zip(indices, sides):
                translation[i] = side
            translation_checks.append(lin(borrow_forms(face, translation), baseline))
        for pattern in patterns:
            offsets = [0] * 5
            for i, level in zip(indices, pattern):
                offsets[i] = level
            if not physical(q, offsets):
                continue
            delta = lin(borrow_forms(face, offsets), baseline)
            distinct.setdefault(delta, tuple(offsets))
        assert distinct
        separate_groups.append(tuple(sorted(distinct.items())))
    return baseline, tuple(separate_groups), tuple(translation_checks)


def maximum_phase_work(face, gap):
    baseline, groups, translations = work_forms(face)
    assert all(not dot(t, gap) for t in translations)
    w = dot(baseline, gap)
    offsets = [0] * 5
    for forms in groups:
        delta, pattern = max(forms, key=lambda item: dot(item[0], gap))
        w += dot(delta, gap)
        offsets = [x + y for x, y in zip(offsets, pattern)]
    assert physical(q_from_gaps(face.exemplar), offsets)
    assert dot(borrow_forms(face, offsets), gap) == w
    return w, tuple(offsets)


def analyze_work(word, faces, vertices):
    checks = []
    for start in vertices:
        phase_gaps = start
        phases = []
        for phase in word:
            face = faces[phase]
            value, offsets = maximum_phase_work(face, phase_gaps)
            phases.append((value, offsets))
            phase_gaps = tuple(phase_gaps[i] for i in face.permutation)
        assert phase_gaps == start
        bound = sum(value for value, _ in phases)
        scale = len(word) * max(abs(x - Fraction(1, N)) for x in start)
        if scale == 0:
            assert bound == 0
        checks.append({"point": start, "phases": phases, "sum": bound, "scaled_norm": scale,
                       "constant": -bound / scale if scale else None})
    candidates = [v["constant"] for v in checks if v["constant"] is not None]
    return {"constant": min(candidates) if candidates else None, "tests": checks}
