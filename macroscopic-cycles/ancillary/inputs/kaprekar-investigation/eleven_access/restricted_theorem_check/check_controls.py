"""Independently derive a strict inherited path and the marked box profiles."""
import argparse
import json
from fractions import Fraction as Q
from itertools import product
from math import factorial, prod
from pathlib import Path

import cover_exact
from cover_exact import Reader
from exact_model import (CENTRE, IDENTITY, Sheet, ZERO, combine, dot, exit_chart,
                         leading, pivot, pull, sort_step, stage, total)

ACCESS = Path(__file__).resolve().parents[1]
IDENTITY_CHART = tuple(f[:6] for f in IDENTITY)
IDENTITY_BASIS = tuple(f[:5] for f in IDENTITY)


def class_six():
    gaps = tuple(map(Q, ('12379624856777/184832467793424', '82001661370921/539094697730820',
                         '2173374886826093/12938272745539680', '9688666835593/40180971259440',
                         '389387410061627/2156378790923280', '494841812881277/2587654549107936')))
    assert sum(gaps) == 1
    q = tuple(sum(gaps[i:]) for i in range(1, 6))
    chart = IDENTITY_CHART
    phase = q
    for _ in range(5):
        assert all(phase[i] + phase[4-i] > 1 for i in range(5))
        assert all(phase[i+1] + phase[4-i] < 1 for i in range(4))
        a, b, c, d, e = phase
        phase = a, a-e, b-e, b-d, c-d
    assert phase == q
    sheet = Sheet(0, 6, ZERO, ZERO, IDENTITY_BASIS)
    active = []
    weight = Q(1)
    rank_list = [len(sheet.basis)]
    chosen = (5, 1, 5, 9, 16, 5, 6, None)
    for rem in chosen:
        guards, read = stage(chart, q, sheet)
        assert read
        shorts, length, failures, chart, sheet, drift, falling = read
        race = sorted(falling, key=lambda f: leading(f, q))
        assert len(race) == 1 or leading(race[0], q) < leading(race[1], q)
        winning = race[0]
        guard, decrement = falling[winning]
        active += list(guards)
        active += [combine(1, other, -1, winning) for other in race[1:]]
        chart = exit_chart(chart, drift, winning)
        children = list(sheet.exits(guard, decrement, drift))
        rank = {len(child.basis) for _, _, child in children}
        assert len(rank) == 1
        rank_list.extend(rank)
        if rem is not None:
            (index, sheet), = [(i, child) for r, i, child in children
                               if r == rem and (6 - child.base) % child.stride == 0]
            weight /= index
        print('rank control', len(rank_list) - 1, shorts, length, failures,
              tuple(6 * d for d in drift), decrement * 6, rem, rank_list[-1])
    assert rank_list == [5, 4, 3, 2, 2, 2, 1, 1, 0]
    assert weight == Q(1, 155520)
    assert all(chart[i] == (0, 0, 0, 0, 0, CENTRE[i]) for i in range(5))
    # Bound a cube in the original differences by every strict source guard.
    delta = min([min(gaps) / 4, *(leading(f, q) / (2 * sum(map(abs, f[:5])))
                                 for f in active if any(f[:5]))])
    lower = weight * factorial(11) * (2 * delta) ** 5 * prod(g / 2 for g in gaps)
    assert lower > 0
    print('rank control strict source half-width', delta, 'mass lower bound', lower)


BOX = tuple(tuple(Q(t, 10**5) for t in ends) for ends in
            ((94040, 94130), (73770, 73850), (57940, 58015), (20, 80), (200, 300)))
BOX_VERTICES = tuple(product(*BOX))
BOX_POINT = tuple((a+b)/2 for a, b in BOX)
ANCHOR = IDENTITY_CHART[0]
CONSTANT = (0, 0, 0, 0, 0, 1)
A_COUNT = (Q(1, 22), 0, Q(1, 22), Q(-1, 4), Q(-5, 22), Q(-3, 44))
P_FORM = combine(1, (0, 0, 1, Q(-1, 2), 0, 0), -2, A_COUNT)
R_FORM = combine(1, (0, 1, 0, Q(-3, 2), -1, 0), -4, A_COUNT)
Q_R_CHART = (combine(-1, P_FORM, Q(3, 2), CONSTANT),
             R_FORM, P_FORM, combine(-1, P_FORM, 1, CONSTANT),
             combine(1, R_FORM, Q(-1, 2), CONSTANT))


ORDER_WORDS = dict(S='5061728394a', U='5061782394a', V='5061273894a',
                   W='5067182934a', C='5067128934a')


def guarded_switch(chart, sheet, name):
    for symbol in name:
        q = tuple(leading(f, BOX_POINT) for f in chart)
        def sign(f):
            value = leading(f, q) or total(f, sheet.point, sheet.base)
            return (value > 0) - (value < 0)
        word, guards, _, order = sort_step(IDENTITY, sign)
        assert ''.join(format(j, 'x') for j in order) == ORDER_WORDS[symbol]
        assert all(min(leading(pull(g, chart), x) for x in BOX_VERTICES) >= 0 for g in guards)
        chart = tuple(pull(f, chart) for f in word)
        sheet = sheet.push(word)
    return chart, sheet


def check_return(start_chart, start_sheet, end_chart, end_sheet, drift, name, failures=()):
    """Certify the pinned long word at first and last repetitions."""
    current = IDENTITY
    initial = tuple(leading(f, BOX_POINT) for f in start_chart)
    end_point = combine(1, end_sheet.point, -1, drift)

    def sign(f):
        value = leading(f, initial) or total(f, start_sheet.point, start_sheet.base)
        return (value > 0) - (value < 0)

    for phase, symbol in enumerate(name):
        current, guards, paired, order = sort_step(current, sign)
        assert paired == (phase not in failures)
        assert ''.join(format(j, 'x') for j in order) == ORDER_WORDS[symbol]
        for guard in guards:
            for chart, sheet, point in ((start_chart, start_sheet, start_sheet.point),
                                        (end_chart, end_sheet, end_point)):
                form = pull(guard, chart)
                assert all(leading(form, x) >= 0 for x in BOX_VERTICES)
                if not any(form):
                    assert total(guard, point, sheet.base) >= 0
                    assert not dot(guard, (*sheet.slope, sheet.stride))
                    assert all(not dot(guard, c) for c in sheet.basis)
    assert tuple(pull(f, start_chart) for f in current) == start_chart
    for row, d, own in zip(current, drift, IDENTITY):
        displacement = combine(1, row, -1, own)
        assert total(displacement, start_sheet.point, start_sheet.base) == d
        assert not dot(displacement, (*start_sheet.slope, start_sheet.stride))
        assert all(not dot(displacement, c) for c in (*start_sheet.basis, drift))


def early_box_prefix():
    """Check all repeated failures on cohort 0, overshoot one and first pin 17."""
    q_r = Sheet(5, 6, (23, 0, 0, 5, -1), (9, 0, 0, 6, -3),
                ((-1, 0, 1, -1, 0), (0, 1, 0, 0, 1)))
    late = (-20, -4, -2, 2, -4)
    prepare = (-10, -2, 0, 2, -2)
    regular = (-5, -4, -3, -2, -1)
    align = (-5, 0, 1, 2, -1)
    n = (0, 0, 0, Q(1, 6), 0, 0)
    m = (0, 0, 0, Q(5, 12), Q(1, 2), 0)
    assert all(leading(f, p) > 0 for f in (n, m, A_COUNT) for p in BOX_VERTICES)

    def back(chart, sheet, count, drift):
        reverse = tuple(-d for d in drift)
        return exit_chart(chart, reverse, count), Sheet(sheet.base, sheet.stride, sheet.point,
                                                        sheet.slope, (*sheet.basis, reverse))

    f_chart, q_f = back(Q_R_CHART, q_r, A_COUNT, late)
    a_chart, q_a = back(f_chart, q_f, m, prepare)
    e_chart, q_e = a_chart, Sheet(q_a.base, q_a.stride, combine(1, q_a.point, -1, align),
                                  q_a.slope, q_a.basis)
    original, q_0 = back(e_chart, q_e, n, regular)
    a, b, c, d, e = IDENTITY_CHART
    assert original == (a, b, c, combine(1, CONSTANT, -1, combine(1, c, 1, e)),
                        combine(1, b, Q(-1, 2), combine(1, CONSTANT, 1, d)))
    check_return(original, q_0, e_chart, q_e, regular, 'S'*5)
    arrived, lifted = guarded_switch(e_chart, q_e, 'SUSSS')
    assert arrived == a_chart and lifted == q_a.reduced()
    check_return(a_chart, q_a, f_chart, q_f, prepare, 'SUSSS'*2)
    check_return(f_chart, q_f, Q_R_CHART, q_r, late, 'VUSSSSUSSS'*2, failures=(0, 10))
    basis = q_0.reduced().basis
    assert prod(col[i] for i, col in enumerate(basis)) == 132
    print('marked early box prefix', 'two failures per late word', 'origin index', 132,
          'lowest late count per base', min(leading(A_COUNT, p) for p in BOX_VERTICES))


def box_checkpoint(base, r2):
    """Reconstruct the primary source coset and its first paired checkpoint."""
    r = 2
    q = (3 * ((base - 1) // 2) + 17, r, 0, base, r - (base - 3) // 2)
    sheet = Sheet(base, 72, q, (108, 0, 0, 72, -36),
                  ((-1, 0, 1, -1, 0), (0, 6, 0, 0, 6)))
    routes = json.loads((ACCESS / 'box_mod12/routes.json').read_text())
    (labels,) = [r for r in routes if (r['cohort'], r['z'], r['r1']) == (0, 1, 17)]
    chart, sheet = guarded_switch(Q_R_CHART, sheet, labels['switch'])
    j = (2, 2, 0, -2, 4, -3, 0)
    assert dot(j, sheet.slope + (sheet.stride,)) % 72 == 0
    value, vector, rest = pivot(sheet.basis, j)
    assert value == 2
    diff = (r2 - total(j, sheet.point, base)) % 72
    assert diff % 2 == 0
    initial = combine(1, sheet.point, diff // 2, vector)
    assert total(j, initial, base) % 72 == r2 % 72
    sheet = Sheet(base, sheet.stride, initial, sheet.slope, (*rest, tuple(36 * x for x in vector)))
    drift = (-20, -4, 0, 4, -4)
    time = tuple(v / 72 for v in map(Q, pull(j, chart)))
    count = (total(j, initial, base) - r2) // 72
    end = combine(1, initial, count, drift)
    shift = combine(1, sheet.slope, dot(j, (*sheet.slope, sheet.stride)) // 72, drift)
    cols = tuple(combine(1, c, dot(j, c) // 72, drift) for c in sheet.basis)
    out = exit_chart(chart, drift, time)
    end_sheet = Sheet(base, sheet.stride, end, shift, cols).reduced()
    check_return(chart, sheet, out, end_sheet, drift, labels['paired'])
    return out, end_sheet


def marked_profiles(reader, base, r2):
    chart, sheet = box_checkpoint(base, r2)
    initial = combine(1, ANCHOR, -1, chart[0])
    pending = [(chart, sheet, initial, Q(1))]
    profiles = {}
    max_count = 0
    while pending:
        chart, sheet, elapsed, weight = pending.pop()
        max_count += 1
        assert max_count < 5000
        kind, read, guards = reader(chart, BOX_POINT, sheet)
        assert kind == "stage", (kind, base, r2, chart, sheet)
        assert all(min(leading(f, x) for x in BOX_VERTICES) >= 0 for f in guards)
        shorts, length, failures, chart, sheet, drift, falling = read
        assert not failures
        time = min(falling, key=lambda f: leading(f, BOX_POINT))
        assert all(min(leading(combine(1, other, -1, time), x) for x in BOX_VERTICES) >= 0
                   for other in falling)
        guard, decrement = falling[time]
        elapsed = combine(1, elapsed, length, time)
        out = exit_chart(chart, drift, time)
        if not any(f[:5] != ZERO for f in out):
            assert tuple(f[5] for f in out) == CENTRE
            profiles[elapsed] = profiles.get(elapsed, Q(0)) + weight
            continue
        kids = [(index, child) for _, index, child in sheet.exits(guard, decrement, drift)
                if (base - child.base) % child.stride == 0]
        assert kids and sum(Q(1, index) for index, _ in kids) == 1
        pending.extend((out, child, elapsed, weight/index) for index, child in kids)
    assert sum(profiles.values()) == 1
    candidates = [tuple(map(Q, x)) for x in (('29/54','-38/27','-20/27'),
                  ('139/144','-133/72','-35/36'), ('863/918','-836/459','-440/459'),
                  ('-11/324','-133/162','-35/81'))]
    result = {}
    for row, mass in profiles.items():
        (coeff,) = [c for c in candidates if combine(1, ANCHOR, 1, combine(c[1], P_FORM, c[2], R_FORM))
                    == combine(1, row, -c[0], (0, 0, 0, 0, 0, 1))]
        result[coeff[0]] = mass
    return result


def box(reader):
    single = {Q(29, 54): Q(1)}
    short = {Q(-11, 324): Q(1)}
    long = {Q(139, 144): Q(2, 3), Q(863, 918): Q(1, 3)}
    expected = {5: [single, long], 11: [long, short],
                17: [long, single], 23: [short, long]}
    pools = {}
    for base in range(5, 72, 6):
        profiles = [marked_profiles(reader, base, r2) for r2 in (603, 639)]
        assert profiles == expected[base % 24]
        print('box', base, profiles)
        pool = tuple(sorted((k, (profiles[0].get(k, 0) + profiles[1].get(k, 0)) / 2)
                            for k in profiles[0].keys() | profiles[1].keys()))
        mod = base % 12
        if mod in pools:
            assert pools[mod] == pool
        pools[mod] = pool
    print('pooled box', pools)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    args = parser.parse_args()
    class_six()
    early_box_prefix()
    cover_exact.SHORT_CAP = 500
    reader = Reader(args.executable)
    box(reader)
    reader.close()
