#!/usr/bin/env python3
"""Check the length-35 certificate from the paired gap update and prepare controls."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from math import factorial, prod
from pathlib import Path


N = 35
K = 17
M = K + 1
PERM = tuple(min(2 * j, N - 2 * j) for j in range(M))
ROOT = []
slot = 1
while slot not in ROOT:
    ROOT.append(slot)
    slot = PERM[slot]
assert slot == 1
R = len(ROOT)
J = set(ROOT) | {0}
COMPLEMENT = sorted(set(range(M)) - J)
DIM = range(R)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b, strict=True))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def scale(scalar, x):
    return tuple(scalar * v for v in x)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def sup(x):
    return max(map(abs, x))


def comparison(gaps, row):
    return sum(gaps[-row:]) - sum(gaps[:row])


def folds(support):
    unseen = set(support)
    cycles = []
    while unseen:
        position = min(unseen)
        cycle = []
        while position in unseen:
            cycle.append(position)
            unseen.remove(position)
            position = PERM[position]
        assert position == cycle[0]
        cycles.append(cycle)
    return cycles


# Relative coordinates are in physical fold order throughout the replay.
# At a phase, the label attached to index j is (phase + R - 1 - j) mod R.
# Derive physical reads from prefix/suffix rows before assigning these labels.
READ_POSITIONS = []
for row in (1, 2):
    signed = {j: 1 for j in range(M - row, M)}
    signed.update({j: -1 for j in range(row)})
    assert set(signed) <= J and sum(signed.values()) == 0
    READ_POSITIONS.append(tuple((ROOT.index(j), c) for j, c in signed.items() if j))


def form_at(phase, row_index):
    result = [0] * R
    for j, c in READ_POSITIONS[row_index]:
        result[(phase + R - 1 - j) % R] += c
    return tuple(result)


FORMS = [form_at(phase, row) for row in range(2) for phase in DIM]


def branch_threshold(index):
    return 1 if index < R else -1


def make_kick(changes):
    physical = dict(changes)
    return tuple(physical.get(j, 0) - physical.get(0, 0) for j in ROOT)


# Primitive borrowing correction and the two differences of adjacent h rows.
# An h_i has effects across the end of pair i and the start of pair i+1.
EFFECT_SIGN = make_kick([(0, 1), (1, -2), (2, 1)])
EFFECT_H1, EFFECT_H2 = [
    make_kick(zip(range(2 * i - 2, 2 * i + 2), (-1, 1, 1, -1)))
    for i in (1, 2)
]


def step(x, y, z):
    sign = (y > 0) - (y < 0)
    h1 = max(0, -1 - y)
    h2 = max(0, -1 - z)
    rotated = (x[-1],) + x[:-1]
    relative = tuple(v + sign * s + h1 * a + h2 * b for v, s, a, b in
                     zip(rotated, EFFECT_SIGN, EFFECT_H1, EFFECT_H2, strict=True))
    return relative, h2


def observations(x):
    return [sum(x[j] * c for j, c in read) for read in READ_POSITIONS]


def trace(x, steps, tracked):
    assert steps % R == 0
    minima = [None] * (2 * R)
    evaluated = [0, 0, 0, 0, 0]
    largest = sup(x)
    exported = 0
    for t in range(steps):
        phase = t % R
        y, z = observations(x)
        for index, value in ((phase, y), (phase + R, z)):
            if minima[index] is None or value < minima[index]:
                minima[index] = value
        if phase in tracked:
            evaluated[(y > 0) - (y < 0) + 1] += 1
        else:
            y = 1
        if phase + R in tracked:
            evaluated[3 + (z >= -1)] += 1
        else:
            z = -1
        x, h2 = step(x, y, z)
        largest = max(largest, sup(x))
        exported += h2
    return x, minima, largest, exported, evaluated


class Span:
    """Exact elimination of a small set of rational row vectors."""

    def __init__(self):
        self.basis = {}

    def reduce(self, vector):
        current = tuple(map(Fraction, vector))
        for pivot, row in sorted(self.basis.items()):
            current = subtract(current, scale(current[pivot], row))
        return current

    def insert(self, vector):
        row = self.reduce(vector)
        pivot = next(i for i, a in enumerate(row) if a)
        self.basis[pivot] = scale(1 / row[pivot], row)

    def indices(self):
        return {i for i, form in enumerate(FORMS) if not any(self.reduce(form))}


def terminal(x):
    seen = {}
    path = []
    losses = []
    while x not in seen:
        seen[x] = len(path)
        path.append(x)
        x, loss = step(x, *observations(x))
        losses.append(loss)
    entry = seen[x]
    assert not any(losses[entry:])
    return entry, path[entry:], sum(losses[:entry]), max(map(sup, path))


def certify_route(certificate):
    span = Span()
    tracked = set()
    x = (0,) * R
    stages = certificate['stages']
    assert len(stages) == R
    audited = []
    determinant = 1
    for stage in stages:
        before = x
        assert stage['prefix'] >= 0 and stage['period'] > 0
        a, low_prefix, size_prefix, h_prefix, prefix_counts = trace(before, stage['prefix'], tracked)
        b, low_word, size_word, h_word, word_counts = trace(a, stage['period'], tracked)
        ray = subtract(a, b)[::-1]
        shift = tuple(stage['shift'])
        assert ray == tuple(stage['ray'])
        assert all(dot(FORMS[i], ray) == dot(FORMS[i], shift) == 0 for i in tracked)
        assert low_prefix == stage['prefix_minimum'] and low_word == stage['word_minimum']
        assert (h_prefix, h_word) == (stage['prefix_export'], stage['word_export'])
        size = max(size_prefix, size_word)
        assert size == stage['largest_nominal']
        x = add(b, shift[::-1])
        span.insert(FORMS[stage['added']])
        next_tracked = span.indices()
        fresh = sorted(next_tracked - tracked)
        assert fresh == stage['fresh']
        slopes = {f: dot(FORMS[f], ray) for f in fresh}
        assert min(slopes.values()) > 0
        assert all(low_word[f] + dot(FORMS[f], shift) >= branch_threshold(f) for f in fresh)
        determinant *= slopes[stage['added']]
        audited.append(dict(ray=ray, shift=shift, fresh=fresh, prefix=low_prefix, word=low_word,
                            size=size, hp=h_prefix, hw=h_word,
                            counts=add(prefix_counts, word_counts)))
        tracked = next_tracked
    assert tracked == set(range(2 * R))
    assert determinant == certificate['determinant_in_read_coordinates']
    print(f'read-coordinate determinant {determinant}')

    # Solve floors from sufficient affine inequalities. In earlier copies, a
    # negative coefficient can have magnitude below (floor+1)*L, while a positive
    # one can be dropped. The form being installed has a strictly positive slope.
    for j, audit in enumerate(audited):
        needed = [0]
        for f in audit['fresh']:
            read = FORMS[f]
            slope = dot(read, audit['ray'])
            for i in range(j + 1):
                previous = audited[i]
                segments = ['prefix', 'word'] if i < j else ['prefix']
                for segment in segments:
                    minimum = previous[segment][f]
                    if minimum is None:
                        continue
                    translated = minimum + sum(dot(read, v['shift']) for v in audited[i:j + 1])
                    worst_earlier = sum((v['floor'] + 1) * min(0, dot(read, v['ray']))
                                        for v in audited[i:j])
                    residual = max(0, branch_threshold(f) - translated) - worst_earlier
                    needed.append(-(-residual // slope))
        audit['floor'] = max(needed)
        assert audit['floor'] == stages[j]['coefficient_floor']

    # Independently check every forced guard at the exact worst endpoints of the
    # offsets throughout each segment. Each worst slack has form A*L + C.
    # Nonnegative A and A+C establish it over every integer L >= 1.
    for i, audit in enumerate(audited):
        minimum_slack = None
        guard_cases = 0
        for j in range(i, R):
            for f in audited[j]['fresh']:
                read = FORMS[f]
                for segment in ('prefix', 'word'):
                    low = audit[segment][f]
                    if low is None:
                        continue
                    slope, constant = 0, low - branch_threshold(f)
                    for k in range(i, j + 1):
                        direction = audited[k]
                        linear = dot(read, direction['ray'])
                        constant += dot(read, direction['shift'])
                        if linear < 0:
                            slope += (direction['floor'] + 1) * linear
                            constant -= linear
                        elif k > i or segment == 'prefix':
                            slope += direction['floor'] * linear
                    assert slope >= 0 and slope + constant >= 0, (i, j, segment, f)
                    if minimum_slack is None or slope + constant < minimum_slack:
                        minimum_slack = slope + constant
                    guard_cases += 1
        checked = sum(value is not None for segment in ('prefix', 'word')
                      for index, value in enumerate(audit[segment])
                      if index in {f for a in audited[i:] for f in a['fresh']})
        assert checked == guard_cases
        primary_exact, secondary_exact = sum(audit['counts'][:3]), sum(audit['counts'][3:])
        steps = stages[i]['prefix'] + stages[i]['period']
        print(f"stage {i+1:2}  new {audit['fresh']!s:24}  K {audit['floor']:2}"
              f"  guards {guard_cases:2} min slack {minimum_slack:3}"
              f"  exact reads {primary_exact:6}/{secondary_exact:<6}"
              f"  forced reads {steps-primary_exact:6}/{steps-secondary_exact:<6}")

    assert x == tuple(certificate['terminal_source'])
    entry, cycle, export, terminal_size = terminal(x)
    assert (entry, len(cycle), export) == (certificate['terminal_entry_from_release'],
                                          certificate['terminal_period'],
                                          certificate['terminal_export_to_entry'])
    assert cycle == list(map(tuple, certificate['terminal_cycle']))
    print(f'release {x}; first entry {entry}; least period {len(cycle)}; transient export {export}')

    max_relative = max(terminal_size, *(a['size'] for a in audited))
    max_relative += sum(sup(a['shift']) + (a['floor'] + 1) * sup(a['ray']) for a in audited)
    spread = 2 * max_relative
    total_loss = export + sum(a['hp'] + (a['floor'] + 2) * a['hw'] for a in audited)
    elapsed = entry + len(cycle) + sum(s['prefix'] + (a['floor'] + 2) * s['period']
                                       for a, s in zip(audited, stages, strict=True))
    assert (spread, total_loss, elapsed) == (certificate['spread_coefficient'],
                                              certificate['export_coefficient'],
                                              certificate['time_through_a_period_coefficient'])
    last = audited[-1]
    assert last['floor'] >= 1 and last['hw'] == 461994 and stages[-1]['period'] == 596448
    assert 7 in last['fresh'] and dot(FORMS[7], last['ray']) == 77454
    print(f'all-scale upper bounds: spread {spread} L; export {total_loss} L; time {elapsed} L')
    return audited


def transport(certificate):
    root_only = [i for i in range(1, M // 2 + 1)
                 if set(range(i)) | set(range(M - i, M)) <= J]
    assert root_only == [1, 2]
    assert folds(COMPLEMENT) == [[5, 10, 15], [7, 14]]
    print(f'root orbit {ROOT}, complementary orbits {folds(COMPLEMENT)}, root-only rows {root_only}')
    profile = tuple(Fraction(3, 55) if j in J else Fraction(16, 275) for j in range(M))
    assert sum(profile) == 1 and all(profile[j] == profile[PERM[j]] for j in range(M))
    rows = [comparison(profile, i) for i in range(1, M // 2 + 1)]
    delta = min(profile) - Fraction(1, M + 1)
    radii = [abs(sum(profile[i:]) - Fraction(M - i, M)) for i in range(1, M)]
    radius = max(radii)
    assert radius == Fraction(1, 198) and radii.index(radius) + 1 == 5
    assert (profile, rows, delta, radius) == (
        tuple(map(Fraction, certificate['lift_profile'])),
        list(map(Fraction, certificate['lift_rows'])),
        Fraction(certificate['lift_delta_pair']), Fraction(certificate['lift_radius']))
    gamma = min(rows[2:])
    budget = max(certificate['spread_coefficient'], certificate['export_coefficient'])
    assert budget == certificate['common_coefficient']
    error = M + 2 * budget + 1
    d = 4 * error / min(delta, gamma / M, radius / M)
    assert d == certificate['lift_D'] and d.denominator == 1
    d = int(d)
    threshold = 2 * d
    assert threshold == certificate['certified_base_threshold']

    # At every B >= threshold use EL/B <= E/D. Check positive pair margin
    # (covers the borrow offset), each individual remaining row, and the cut.
    deviation = Fraction(error, d)
    min_gap = min(profile) - deviation
    assert (M + 1) * min_gap * threshold >= threshold + 1
    assert min_gap * threshold >= 2
    assert all(row - 2 * i * deviation > 0 for i, row in enumerate(rows, 1) if i >= 3)
    cut = Fraction(certificate['lift_cut'])
    assert cut == radius / 2 and radius - (M - 5) * deviation > cut
    assert min_gap - Fraction(1, threshold) >= Fraction(1, 38)
    assert factorial(N) == certificate['original_probability_lower_bound_numerator']
    assert threshold ** K * 38 ** M == certificate['original_probability_lower_bound_denominator']
    print(f'd rows {[str(a) for a in rows]}, pairing slack {delta}, radius {radius}, D {d}')
    print(f'word bound {factorial(N) / (threshold ** K * 38 ** M):.12e}')
    return profile, d


def controls(audited, profile, d):
    # Bases and offsets chosen here, independently of the published controls.
    print('# base L (18 physical gaps, slots 0 through 17)')
    for base, sample in ((2 * d + 1601, 0), (3 * d + 71393, 1)):
        length = base // d
        x = (0,) * R
        for j, a in enumerate(audited):
            variation = (sample + 2 * j + j // 2) % length
            repeat = a['floor'] * length + variation
            x = add(x, add(a['shift'], scale(repeat, a['ray']))[::-1])
        g = [None] * M
        centre = base * profile[0] - Fraction(sum(x), R + 1)
        g[0] = int(centre // 1) + (sample + 1) % length
        for j, relative in zip(ROOT, x, strict=True):
            g[j] = g[0] + relative
        for index, slot in enumerate(COMPLEMENT[:-1]):
            g[slot] = int((base * profile[slot]) // 1) + (index + sample) % length
        g[COMPLEMENT[-1]] = base - sum(value for value in g if value is not None)
        assert min(g) >= 2
        print(base, length, *g)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--controls', type=Path)
    args = parser.parse_args()
    contents = args.certificate.read_bytes()
    print(f'certificate SHA-256 {sha256(contents).hexdigest()}')
    certificate = json.loads(contents)
    audited = certify_route(certificate)
    profile, d = transport(certificate)
    if args.controls:
        from contextlib import redirect_stdout

        with args.controls.open('w') as output, redirect_stdout(output):
            controls(audited, profile, d)


if __name__ == '__main__':
    main()
