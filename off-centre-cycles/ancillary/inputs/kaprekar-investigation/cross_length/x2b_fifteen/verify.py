"""Exact guard arithmetic and ordered-digit replay for Task 157."""
import sys
import csv
import gzip
from collections import Counter
from fractions import Fraction
from itertools import product
from math import floor, gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'x2_seventeen'))
import verify as root_certificate  # noqa: E402

ALPHA = Fraction(47, 400)
BETA = ALPHA + Fraction(1, 50)
SCALE = 100_000
ROOT = (0, 1, 2, 4, 7)
COMPLEMENT = (3, 5, 6)


def premise():
    results = tuple(map(root_certificate.check_route, root_certificate.ROUTES))
    assert results == ((4, 0, 177), (16, 4, 179), (52, 12, 155), (48, 18, 124), (52, 13, 72))
    root_certificate.tail()
    base = 1_111_013
    repeats = (3, 4, 6, 8)
    offsets = (9, 6, 2, 10)
    root_certificate.check_physical(base, (3, 4, 5, 7), offsets)
    error = Fraction(191, SCALE)
    assert Fraction(14, 135) - error > Fraction(101, 1000)
    assert Fraction(1, 60) - 8*error > 0
    assert Fraction(1, 9) - Fraction(14, 135) - error > Fraction(1, 200)
    print('Task 155 exact premise replay:', base, repeats, offsets, 'four repeat counts, guards, and margins')


def differences(g):
    return tuple(sum(g[i:]) for i in range(1, len(g)))


def gaps(word, base):
    digits = sorted(word)
    q = tuple(hi - lo for lo, hi in zip(digits[:len(digits)//2], reversed(digits)))
    return (base - q[0], *(a - b for a, b in zip(q, q[1:])), q[-1])


def literal(word, base):
    digits = sorted(word)
    value = 0
    for lo, hi in zip(digits, reversed(digits)):
        value = value*base + hi - lo
    image = []
    for _ in digits:
        value, digit = divmod(value, base)
        image.append(digit)
    assert value == 0
    return tuple(reversed(image))


def fold(i, m):
    return min(2*i, 2*m - 1 - 2*i)


def secondary(g):
    return tuple(sum(g[-i:]) - sum(g[:i]) for i in range(2, len(g)//2 + 1))


def paired_slack(g):
    q = differences(g)
    base = sum(g)
    k = len(q)
    pairs = [(q[i] - int(i == k-1), base - q[k-1-i] - int(i < k-1)) for i in range(k)]
    return min(min(left) - max(right) for left, right in zip(pairs, pairs[1:]))


def profile_gaps(base, rel, offsets):
    g = [0] * 8
    anchor = floor(ALPHA*base - Fraction(sum(rel), 5)) + offsets[0]
    for i, d in zip(ROOT, rel):
        g[i] = anchor + d
    for i, u in zip(COMPLEMENT[:-1], offsets[1:]):
        g[i] = floor(BETA*base) + u
    g[COMPLEMENT[-1]] = base - sum(g)
    return tuple(g)


def lifted_step(g):
    z = g[-1] - g[0]
    e = (z > 0) - (z < 0)
    h = max(0, -1 - z)
    image = [0] * len(g)
    for i, d in enumerate(g):
        image[fold(i, len(g))] = d
    image[0] += e - h
    image[1] += -2*e + h
    image[2] += e + h
    image[3] -= h
    return tuple(image), h


def physical15(base, parameters, offsets):
    L = base//SCALE
    rel = (0, *root_certificate.initial_root(parameters))
    g = profile_gaps(base, rel, offsets)
    word = root_certificate.word_with_gaps(g)
    assert len(set(word)) == 15 and 0 <= word[0] < word[-1] < base
    r0, r1, p2, p3 = parameters
    repetitions = (r0, r1, 1, p2 + 1, p3 + 1)
    words = ''.join(route.word * r for route, r in zip(root_certificate.ROUTES, repetitions))
    entrance = len(words) + 35
    seen = {}
    loss = 0
    for phase in range(entrance + 8):
        assert g not in seen and g == gaps(word, base)
        seen[g] = phase
        assert 9*min(g) > base + 1 and paired_slack(g) > 0 and min(secondary(g)) > 0
        assert 200*(base - 8*g[0]) > 8*base
        x = tuple(g[i] - g[0] for i in ROOT[1:])
        assert root_certificate.spread(x) <= 179*L
        assert all(abs(g[i] - ALPHA*base) < 191*L for i in ROOT)
        assert all(abs(g[i] - BETA*base) <= 60*L for i in COMPLEMENT)
        if phase < len(words):
            assert root_certificate.branch(x[-1]) == words[phase]
        elif phase == len(words):
            assert x == (0, 0, 18, 2)
        g, h = lifted_step(g)
        word = literal(word, base)
        if phase >= entrance:
            assert not h
        loss += h
    assert g == gaps(word, base) and seen[g] == entrance
    assert min(tuple(d[i] - d[0] for i in ROOT[1:]) for d, p in seen.items() if p >= entrance) == (-1, 2, -1, 0)
    assert loss == 4*r1 + 18*(p2 + 1) + 13*(p3 + 1) + 22 <= 53*L
    cycle = word
    for phase in range(1, 9):
        word = literal(word, base)
        assert (word == cycle) == (phase == 8)


def fifteen():
    profile = tuple(ALPHA if i in ROOT else BETA for i in range(8))
    assert sum(profile) == 1 and secondary(profile) == (Fraction(1, 50), Fraction(1, 25), Fraction(1, 50))
    assert tuple(fold(i, 8) for i in ROOT) == (0, 2, 4, 7, 1)
    assert tuple(fold(i, 8) for i in COMPLEMENT) == (6, 5, 3)
    error = Fraction(191, SCALE)
    assert 9*(ALPHA - error)*SCALE > SCALE + 1
    assert 8*error < Fraction(1, 50)
    assert Fraction(1, 8) - ALPHA - error == Fraction(559, 100_000) > Fraction(1, 200)
    assert (ALPHA - error - Fraction(1, 11))*SCALE > 1
    physical15(SCALE, (0,)*4, (0,)*3)
    for b in range(2*SCALE, 2*SCALE + 16):
        for p in product((0, 1), repeat=4):
            physical15(b, p, (b % 2, 1, 0))
    physical15(10_000_003, (72, 86, 99, 97), (92, 45, 0))
    print('15-digit exact profile:', ALPHA, BETA, 'secondary margins', secondary(profile))
    print('15-digit literal-subtraction lifts: 258 across all bases mod 16, including the exact lower endpoint')


def translations():
    root = (0, 1, 2, 4, 8)
    assert gcd(40, 49) == 1
    outside = tuple(2000 if i in root else 1900 for i in range(9))
    assert secondary(outside) == (-100, -200, -200) and paired_slack(outside) > 0
    prior = gaps(literal((0,)*9 + differences(outside), sum(outside)), sum(outside))
    for a, base_shift in ((4, 40), (5, 49)):
        t = tuple(a + (i not in root) for i in range(9))
        assert sum(t) == base_shift and all(t[i] == t[fold(i, 9)] for i in range(9))
        assert secondary(t) == (1, 2, 2) and 10*min(t) >= base_shift
        shifted = tuple(d + u for d, u in zip(outside, t))
        image = gaps(literal((0,)*9 + differences(shifted), sum(shifted)), sum(shifted))
        mismatch = tuple(y - x - u for x, y, u in zip(prior, image, t))
        assert mismatch == (0, 0, 1, -1, 1, -1, 0, 0, 0)
        print('17-digit base shift:', base_shift, 'secondary increment', secondary(t), 'outside-chamber mismatch', mismatch)


def outside_predecessor(g):
    # Inverse of the paired branch that puts b before a in pairs 4 and 5.
    a, b, c, d, e, f, h, i, j = g
    exchange = a - b + c - d + e - f - 5
    return (a - 1, c - 1, e, i - exchange, j, h + exchange, f, d, b + 2)


def branch_image(g, borrow):
    base = sum(g)
    q = differences(g)
    ordered = [base - int(borrow)]
    for i in range(8):
        a = q[i] - int(borrow and i == 7)
        b = base - q[7-i] - int(borrow and i < 7)
        ordered.extend((b, a) if i in (3, 4) else (a, b))
    image = tuple(ordered[i] - ordered[16-i] for i in range(8))
    return (base - image[0], *(a - b for a, b in zip(image, image[1:])), image[-1]), ordered


def outside_feeder():
    for original in [(0,)*9] + [tuple(int(i == j) for i in range(9)) for j in range(9)]:
        predecessor = outside_predecessor(original)
        assert sum(original) == sum(predecessor)
        assert branch_image(predecessor, True)[0] == original
        assert predecessor[-1] - predecessor[0] == original[1] - original[0] + 3

    alpha, beta = root_certificate.ALPHA, root_certificate.BETA
    difference = beta - alpha
    profile = (alpha, alpha, alpha, alpha + 3*difference, alpha, alpha - difference, beta, beta, alpha)
    predicted, levels = branch_image(profile, False)
    assert predicted == tuple(alpha if i in root_certificate.ROOT else beta for i in range(9))
    assert min(levels[2*i] - levels[2*i+1] for i in range(8)) == Fraction(19, 270)
    assert tuple(levels[2*i+1] - levels[2*i+2] for i in range(8)) == (
        0, difference, 2*difference, 2*difference, 2*difference, 2*difference, difference, 0)
    assert secondary(profile) == (difference, 2*difference, -2*difference)
    limit = 10_000_000
    error = Fraction(1342, limit)
    assert (min(profile) - error - Fraction(1, 12))*limit > 1
    assert (Fraction(19, 270) - 18*error)*limit > 1
    assert (difference - 18*error)*limit > 1
    assert (2*difference - 8*error)*limit > 1
    b = limit + 17
    target = root_certificate.initial_gaps(b, (0,)*4, (0,)*4)
    predecessor = outside_predecessor(target)
    predicted, levels = branch_image(predecessor, True)
    assert predicted == target and secondary(predecessor)[2] < -1
    assert all(x > y for x, y in zip(levels, levels[1:]))
    assert gaps(literal(root_certificate.word_with_gaps(predecessor), b), b) == target
    print('17-digit outside feeder: affine inverse, exact margins and literal one-step replay at base', b)


def archived_small_radius():
    sets = []
    for scale, expected in ((2, (139, 133)), (3, (147, 135))):
        identities = Counter()
        with gzip.open(root_certificate.ARCHIVE / f'raw/n17-s{scale}.tsv.gz', 'rt') as source:
            for row in csv.DictReader(source, delimiter='\t'):
                base = int(row['base'])
                radius = int(row['radius'])
                if int(row['period']) != 8 or 1000*radius > 9*base:
                    continue
                center = tuple(map(int, row['centered_cycle'].split(',')))
                q = tuple((x + (8-i)*base)//9 for i, x in enumerate(center))
                first = q
                assert root_certificate.centered(q, base) == center
                roots, phases = [], []
                for _ in range(8):
                    g = root_certificate.gaps(q, base)
                    if not (q[-1] > 0 and paired_slack(g) >= 0 and min(secondary(g)) >= -1):
                        break
                    roots.append(root_certificate.root(g))
                    phases.append(root_certificate.centered(q, base))
                    q = differences(gaps(literal((*q, *(0,)*9), base), base))
                if len(set(roots)) != 8 or min(roots) != root_certificate.SIGNATURE:
                    continue
                assert q == first and center == min(phases)
                assert radius == max(abs(x) for phase in phases for x in phase)
                assert center[0] % 9 == base % 9
                identities[(radius, center)] += 1
        microscopic = {key for key in identities if key[0] <= 53}
        assert len(microscopic) == 7
        assert {center[0] % 9 for _, center in microscopic} == {3, 4, 6, 7}
        result = sum(identities.values()), sum(identities[key] for key in microscopic)
        assert result == expected
        sets.append(microscopic)
        print('Task 150 archive scale', scale, 'small-radius family hits', result[0],
              'of which seven fixed centered identities contribute', result[1])
    assert sets[0] == sets[1]
    for radius, center in sorted(sets[0]):
        print('  class', center[0] % 9, 'radius numerator', radius, 'centered identity', center)


if __name__ == '__main__':
    premise()
    fifteen()
    translations()
    outside_feeder()
    archived_small_radius()
