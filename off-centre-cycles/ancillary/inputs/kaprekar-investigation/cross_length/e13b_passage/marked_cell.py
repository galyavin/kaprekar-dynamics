#!/usr/bin/env python3
"""A full-dimensional class-0 route into one already certified outer key."""

from fractions import Fraction as Q
from collections import defaultdict
import csv
import gzip
from math import factorial, prod
from pathlib import Path
import sys

from gate_cell import gaps, guards, slots, update
from premise import add, advance, centered, differences


B0 = 140000
SEED = (90594, 64043, 49868, 17220, 10073, 13923)
OUTER = (-1988, 1484, 1491, -1974, 315, 322, 350)
START, END = 13, 14108
WALLS = ((13, 6, 1370, 0, 15), (8233, 12, 8, 2, 14), (8354, 30, 42, 7, 15),
         (9616, 18, 56, 11, 18), (10724, 36, 89, 12, 15))
LAST_START, LAST_LENGTH = 13984, 18
INCREMENTS = ((42, 35, 28, 21, 14, 7), (84, 70, 7, -7, 28, 14),
              (210, 28, 91, -91, -77, -112), (126, 105, 84, -84, 42, 21),
              (147, 84, 28, -28, -84, 63), (84, 21, 7, -7, -21, 63))
RATES = (Q(9584, 7), Q(169, 21), Q(902, 21), Q(1192, 21), Q(826, 9), Q(2063, 63))
WIDTH_DIVISOR = 12


def replay():
    q = differences(B0, gaps(SEED))
    path = {}
    masks = {}
    for t in range(START, END + 1):
        g = centered(B0, q)
        x = tuple(-sum(g[:i]) for i in range(1, 7))
        x_slots = slots(x, 0, 7)
        mask = sum((x_slots[2*i+1] < x_slots[2*i+2]) << i for i in range(6))
        assert min(guards(B0, x, mask, 7)) >= 0
        q = advance(B0, q)
        assert gaps(update(x, mask, 7)) == centered(B0, q)
        if t < END:
            assert min(g) < -2000
        path[t], masks[t] = x, mask
    assert gaps(path[END]) == OUTER
    return path, masks


def follow(x, word):
    for m in word:
        x = update(x, m, 0)
    return x


def excess_slope(slope, perturbations):
    return slope - sum(abs(c) * r for c, r in zip(perturbations, RATES)) / WIDTH_DIVISOR


def slack_positive(offset, slope, perturbations):
    return offset >= 0 and excess_slope(slope, perturbations) >= 0


def exact_forms(form, offset, alpha, *columns):
    at_seed = form(B0, offset, 7)
    at_slope = form(B0, alpha, 0)
    directions = tuple(form(0, c, 0) for c in columns)
    return zip(at_seed, at_slope, *directions)


def valid_word(path, masks, at, word, directions, outside=True):
    for j, m in enumerate(word):
        assert m == masks[at + j]
        x = path[at + j]
        rows = exact_forms(lambda b, a, d: guards(b, a, m, d), x, *directions)
        assert all(slack_positive(a, b, c) for a, b, *c in rows)
        if outside:
            rows = exact_forms(lambda b, a, d: tuple(-g for g in gaps(a)), x, *directions)
            assert any(slack_positive(a - 2001, b, c) for a, b, *c in rows)
        directions = tuple(update(v, m, 0) for v in directions)
    return directions


def verify_route():
    path, masks = replay()
    word_at = lambda start, length: tuple(masks[t] for t in range(start, start + length))
    columns = (SEED, *INCREMENTS)
    bridge = START
    for i, (t, length, copies, phase, row) in enumerate(WALLS):
        columns = valid_word(path, masks, bridge, word_at(bridge, t-bridge), columns)
        word = word_at(t, length)
        drift = add(path[t + length], path[t], -1)
        mask = word[phase]
        h = guards(B0, path[t+phase], mask, 7)[row]
        decrement = -guards(0, follow(drift, word[:phase]), mask, 0)[row]
        assert decrement > 0 and copies == 1 + h // decrement
        assert follow(drift, word) == drift and tuple(follow(c, word) for c in columns) == columns
        rates = tuple(Q(guards(0, follow(c, word[:phase]), mask, 0)[row], decrement) for c in columns)
        assert rates == (RATES[i], *(int(j == i) for j in range(6)))
        columns_at_end = tuple(add(c, drift, r) for c, r in zip(columns, rates))
        if i == 0:
            source_phase = columns
            for m in word:
                assert m == 0
                normalized = guards(B0, source_phase[0], m, 0)
                base_free = [guards(0, c, m, 0) for c in source_phase[1:]]
                assert all(excess_slope(n, perturbations) > 0 for n, *perturbations in zip(normalized, *base_free))
                source_phase = tuple(update(c, m, 0) for c in source_phase)
            x_at_exit = path[t + copies * length][0]
            gate_gap = excess_slope(64 * columns_at_end[0][0] - B0, [64*c[0] for c in columns_at_end[1:]])
            assert gate_gap > 0 and 64 * x_at_exit - B0 - 63 * B0 + 63 * gate_gap > 0
        # Every phase stays outside the outer section during the whole repeat.
        # Use the same witnessing gap at the first and last copies.
        first, last = columns, columns_at_end
        for p, m in enumerate(word):
            assert valid_word(path, masks, t+p, (m,), first, False)
            assert valid_word(path, masks, t+(copies-1)*length+p, (m,), last, False)
            offsets = (path[t+p], path[t+(copies-1)*length+p])
            gap_forms = [list(exact_forms(lambda b, a, s: tuple(-g for g in gaps(a)), x, *cs))
                         for x, cs in zip(offsets, (first, last))]
            assert any(all(slack_positive(h0 - 2001, hn, perturbations)
                           for h0, hn, *perturbations in rows) for rows in zip(*gap_forms))
            first = tuple(update(v, m, 0) for v in first)
            last = tuple(update(v, m, 0) for v in last)
        columns = columns_at_end
        bridge = t + copies * length

    columns = valid_word(path, masks, bridge, word_at(bridge, LAST_START-bridge), columns)
    word = word_at(LAST_START, LAST_LENGTH)
    copies, phase = divmod(END - LAST_START, LAST_LENGTH)
    drift = add(path[LAST_START+LAST_LENGTH], path[LAST_START], -1)
    rising = gaps(follow(drift, word[:phase]))[0]
    before = -2000 - gaps(path[LAST_START+phase])[0]
    assert rising > 0 and copies == (before + rising - 1) // rising
    rates = tuple(-Q(gaps(follow(c, word[:phase]))[0], rising) for c in columns)
    assert rates == (RATES[-1], 0, 0, 0, 0, 0, 1)
    assert follow(drift, word) == drift and tuple(follow(c, word) for c in columns) == columns
    zero = tuple(add(c, drift, rate) for c, rate in zip(columns, rates))
    assert all(not any(c) for c in zero)
    valid_word(path, masks, LAST_START, word, columns, False)
    valid_word(path, masks, LAST_START + (copies-1)*LAST_LENGTH, word, zero, False)
    valid_word(path, masks, LAST_START + copies*LAST_LENGTH, word[:phase], zero, False)
    for j in range(LAST_LENGTH):
        phase_drift = gaps(follow(drift, word[:j]))
        start_gaps = gaps(path[LAST_START+j])
        last_copy = copies - int(j >= phase)
        assert any(g + last_copy * d < -2000 and d >= 0 for g, d in zip(start_gaps, phase_drift))

    assert all(slack_positive(a, b, c) for a, b, *c in exact_forms(
        lambda base, x, shift: tuple(25*(base+g) - 7*base for g in gaps(x)), SEED, SEED, *INCREMENTS))
    # Strictly below B/8 in g0 for large enough B, hence disjoint from Task 158.
    source_bound = B0 + gaps(SEED)[0]
    perturbed_source_bound = source_bound + sum(abs(gaps(c)[0]) * r for c, r in zip(INCREMENTS, RATES)) / WIDTH_DIVISOR
    assert 8 * perturbed_source_bound < 7 * B0
    assert B0 + 8*gaps(SEED)[0] + 63*B0 + 63*(8*perturbed_source_bound - 7*B0) < 0
    assert all(v % 7 == 0 for c in INCREMENTS for v in c)
    assert all(63 * rate == int(63 * rate) for rate in RATES)
    print(f'Class-0 source X={SEED}; six slopes={RATES}; first outer key G={OUTER}.')


def verify_mark(archive):
    with gzip.open(archive / 'encounters.tsv.gz', 'rt') as stream:
        matches = [r for r in csv.DictReader(stream, delimiter='\t')
                   if (r['base'], r['draw_id'], r['section'], r['event']) == ('140000', '502', 'R2000', 'entrance')]
    assert len(matches) == 1
    assert matches[0]['mark'] == '177' and tuple(map(int, matches[0]['tangent'].split(','))) == OUTER
    q = differences(B0, OUTER)
    seen, peak = {}, 0
    while q not in seen:
        seen[q] = len(seen)
        g = centered(B0, q)
        peak = max(peak, -min(g))
        q = advance(B0, q)
        assert len(seen) <= 32768 and peak <= 4096
    assert len(seen) - seen[q] == 177
    print(f'Archived outer key: preperiod={seen[q]}, period=177, suffix peak R={peak}.')


def weight():
    """Exact polynomial integral over the six original count coordinates."""
    widths = tuple(rate / WIDTH_DIVISOR for rate in RATES)
    indices = tuple(range(6))
    poly = {(0,) * 6: Q(1)}
    for gap, *directions in zip(gaps(SEED), *(gaps(v) for v in INCREMENTS)):
        nxt = defaultdict(Q)
        for exponent, c in poly.items():
            nxt[exponent] += c * (B0 + gap) / (7*B0)
            for i, d in enumerate(directions):
                power = tuple(k + (j == i) for j, k in zip(indices, exponent))
                nxt[power] += c * d / (7*B0)
        poly = nxt
    integral = sum(c * prod(w**k / (k+1) for w, k in zip(widths, power))
                   for power, c in poly.items() if not any(k % 2 for k in power))
    volume = prod(2*w / B0 for w in widths)
    return factorial(13) * volume * integral


if __name__ == '__main__':
    verify_route()
    verify_mark(Path(sys.argv[1]))
    mass = weight()
    print(f'Exact class-0 limiting original mass of this marked cell: {mass} (~{float(mass):.6g})')
