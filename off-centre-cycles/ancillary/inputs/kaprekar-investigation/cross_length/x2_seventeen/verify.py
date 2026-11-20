"""Exact affine certificates and independent literal-subtraction controls."""
import csv
import gzip
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import floor
from pathlib import Path

Vec = tuple[int, ...]
N = 17
M = 9
ROOT = (0, 1, 2, 4, 8)
COMPLEMENT = (3, 5, 6, 7)
SIGNATURE = (-1, 2, -1, 0)
SCALE = 100_000
ALPHA = Fraction(14, 135)
BETA = ALPHA + Fraction(1, 60)
ARCHIVE = Path(__file__).resolve().parent.parent / "x1_terminal_map"


@dataclass(frozen=True)
class Route:
    name: str
    word: str
    start: Vec
    rays: tuple[Vec, ...]
    change: Vec


COMMON = (0, 0, 26, 26)
FAR = (0, 0, 27, 0)
ROUTES = (
    Route("four positive", "++++", (5, 13, 107, 62),
          ((5, 5, 5, 5), (0, 12, 12, 12), COMMON, FAR), (-5, -5, -5, -5)),
    Route("one near", "+++-+++++++-++++", (0, 20, 114, 69),
          ((0, 12, 12, 12), COMMON, FAR), (0, -12, -12, -12)),
    Route("bridge", "+++-++++++--+++0++-++++-++0+++-0+++-++-++++-++-++++0",
          (0, 8, 102, 57), (COMMON, FAR), (0, -8, -31, -29)),
    Route("two near", "++-++++-++++++-++++-++++++-++++-++++++-++++-++++",
          (0, 0, 71, 28), (COMMON, FAR), (0, 0, -26, -26)),
    Route("three near", "++-+0++-++++-++0++-++++-0+++-+-0++0-++-+0++-++-+-++0",
          (0, 0, 45, 2), (FAR,), (0, 0, -27, 0)),
)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def step(x, sign):
    b, c, d, z = x
    if sign == "+":
        return (z - 3, b, c - 1, d - 1), 0
    if sign == "-":
        return (1 - z, b - 2 - 2*z, c - z, d - z), -1 - z
    assert sign == "0"
    return (z, b, c, d), 0


def branch(z):
    return "+" if z > 0 else "-" if z < 0 else "0"


def trace(x, word):
    states = [x]
    loss = 0
    for sign in word:
        x, h = step(x, sign)
        states.append(x)
        loss += h
    return states, loss


def spread(x):
    return max((0, *x)) - min((0, *x))


def check_route(route):
    base, loss = trace(route.start, route.word)
    rays = []
    for ray in route.rays:
        states, offset_loss = trace(add(route.start, ray), route.word)
        assert offset_loss == loss
        rays.append(list(map(sub, states, base)))
    strict_margin = []
    for i, sign in enumerate(route.word):
        z = base[i][-1]
        slopes = [ray[i][-1] for ray in rays]
        if sign == "0":
            assert z == 0 and all(d == 0 for d in slopes)
        else:
            direction = 1 if sign == "+" else -1
            strict_margin.append(direction*z)
            assert direction*z >= 1 and all(direction*d >= 0 for d in slopes), (route.name, i)
    assert base[-1] == add(route.start, route.change)
    assert tuple(ray[-1] for ray in rays) == route.rays
    bound = max(spread(x) + sum(spread(ray[i]) for ray in rays) for i, x in enumerate(base))
    print(route.name, len(route.word), min(strict_margin), loss, bound)
    return len(route.word), loss, bound


def tail():
    x = (0, 0, 18, 2)
    states = {}
    losses = []
    for _ in range(43):
        assert x not in states
        states[x] = len(states)
        x, h = step(x, branch(x[-1]))
        losses.append(h)
    pre = states[x]
    period = len(states) - pre
    cycle = tuple(states)[pre:]
    loss = sum(losses)
    bound = max(map(spread, states))
    assert (pre, period, min(cycle), loss, bound) == (35, 8, SIGNATURE, 10, 18)
    assert all(h == 0 for h in losses[pre:])
    print("tail", pre, period, cycle, loss, bound)


def fold(j):
    return min(2*j, N - 2*j)


def difference(word):
    digits = sorted(word, reverse=True)
    return tuple(digits[i] - digits[-1-i] for i in range(M-1))


def literal(word, base):
    digits = sorted(word)
    descending = ascending = 0
    for lo, hi in zip(digits, reversed(digits)):
        descending = base*descending + hi
        ascending = base*ascending + lo
    value = descending - ascending
    image = []
    for _ in range(N):
        value, digit = divmod(value, base)
        image.append(digit)
    assert value == 0
    return tuple(reversed(image))


def gaps(q, base):
    return (base - q[0], *(a-b for a, b in zip(q, q[1:])), q[-1])


def differences(g):
    return tuple(sum(g[j:]) for j in range(1, M))


def secondary(g):
    return tuple(sum(g[-j:]) - sum(g[:j]) for j in range(2, M//2 + 1))


def centered(q, base):
    return tuple(M*d - (M-1-i)*base for i, d in enumerate(q))


def root(g):
    return tuple(g[j] - g[0] for j in ROOT[1:])


def pairs(q, base, offsets=True):
    return tuple((d - int(offsets and i == M-2), base - q[-1-i] - int(offsets and i < M-2))
                 for i, d in enumerate(q))


def adjacent_slack(paired):
    return min(a - b for prev, curr in zip(paired, paired[1:]) for a in prev for b in curr)


def archived_witness():
    base = 1_111_013
    x = (105029, 210058, 315069, 170429, 275467, 163254, 28433, -105038)
    q0 = tuple((d + (M-1-i)*base)//M for i, d in enumerate(x))
    slope = tuple(196 if j in ROOT else 205 for j in range(M))
    increment = sum(slope)
    assert centered(q0, base) == x and increment == 1800
    assert all(slope[fold(j)] == d for j, d in enumerate(slope))
    assert adjacent_slack(pairs(differences(slope), increment, offsets=False)) == 178
    assert min(secondary(slope)) == 9 and slope[-1] == slope[0]
    assert increment - M*slope[0] == 36 and 450*x[0] >= M*base
    for t in (0, 1, 97):
        b = base + increment*t
        q = differences(add(gaps(q0, base), tuple(t*d for d in slope)))
        seen = []
        for _ in range(8):
            seen.append(q)
            q = difference(literal((*q, *(0,)*M), b))
        assert len(set(seen)) == 8 and q == seen[0]
        gs = [gaps(q, b) for q in seen]
        assert min(map(root, gs)) == SIGNATURE and all(g[-1] - g[0] >= -1 for g in gs)
        pair_slack = min(adjacent_slack(pairs(q, b)) for q in seen)
        assert pair_slack >= 55_018 + 178*t
        if t == 0:
            assert pair_slack == 55_018
            assert max(abs(v) for q in seen for v in centered(q, b)) == 315_096
            assert min(min(secondary(g)) for g in gs) == 24_137
    print("archived affine witness B=1111013+1800t, least period 8, radius >= 1/450")


def archived_counts():
    label = "8:-1,2,-1,0"
    with gzip.open(ARCHIVE / "tables/physical-root-families.tsv.gz", "rt") as handle:
        table = [r for r in csv.DictReader(handle, delimiter="\t")
                 if r["n"] == "17" and r["root_cycle"] == label]
    for scale, expected in ((2, 813), (3, 819)):
        with gzip.open(ARCHIVE / f"raw/n17-s{scale}.tsv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        hits = [r for r in rows if r["root_cycle"] == label]
        assert len(rows) == 2304 and len(hits) == expected
        assert all(int(r["period"]) == 8 for r in hits)
        assert all(1000*int(r["radius"]) > M*int(r["base"]) for r in hits)
        for aggregate in (False, True):
            assert expected == sum(int(r["count"]) for r in table
                                   if r["scale"] == str(scale) and (r["class"] == "all") == aggregate)
        print("archived panel", scale, "signature hits", len(hits), "/", len(rows))


def initial_root(p):
    n0, n1, p2, p3 = p
    n2, n3 = p2 + 1, p3 + 1
    v = 5*n0 + 12*n1
    return (5*n0, 8 + v, 49 + 27*n3 + 26*n2 + v, 31 + 26*n2 + v)


def initial_gaps(base, p, u):
    rel = (0, *initial_root(p))
    anchor = floor(ALPHA*base - Fraction(sum(rel), len(ROOT))) + u[0]
    g = [0] * M
    for j, d in zip(ROOT, rel):
        g[j] = anchor + d
    for j, d in zip(COMPLEMENT[:-1], u[1:]):
        g[j] = floor(BETA*base) + d
    g[COMPLEMENT[-1]] = base - sum(g)
    return tuple(g)


def lifted_step(g):
    z = g[-1] - g[0]
    e = (z > 0) - (z < 0)
    h = max(0, -1 - z)
    image = [0] * M
    for j, d in enumerate(g):
        image[fold(j)] = d
    image[0] += e - h
    image[1] += -2*e + h
    image[2] += e + h
    image[3] -= h
    return tuple(image), h


def word_with_gaps(g):
    left, right = zip(*((v//2, v - v//2) for v in g[1:]))
    word = [0]
    for d in (*left, *reversed(right)):
        word.append(word[-1] + d)
    return tuple(word)


def check_physical(base, p, u):
    limit = base//SCALE
    n0, n1, p2, p3 = p
    g = initial_gaps(base, p, u)
    word = word_with_gaps(g)
    assert len(set(word)) == N and 0 <= word[0] < word[-1] < base
    phases = (ROUTES[0].word*n0 + ROUTES[1].word*n1 + ROUTES[2].word
              + ROUTES[3].word*(p2+1) + ROUTES[4].word*(p3+1))
    start_cycle = len(phases) + 35
    profile = tuple((ALPHA if j in ROOT else BETA)*base for j in range(M))
    seen = {}
    loss = 0
    for i in range(start_cycle + 8):
        assert g not in seen
        seen[g] = i
        assert g == gaps(difference(word), base)
        assert 10*min(g) > base + 1 and min(secondary(g)) > 0
        assert max(abs(v-d) for v, d in zip(g, profile)) <= 191*limit
        assert 200*abs(centered(differences(g), base)[0]) > M*base
        if i < len(phases):
            assert phases[i] == branch(g[-1] - g[0])
        elif i == len(phases):
            assert root(g) == (0, 0, 18, 2)
        assert spread(root(g)) <= 179*limit
        g, h = lifted_step(g)
        word = literal(word, base)
        if i >= start_cycle:
            assert h == 0
        loss += h
    assert seen[g] == start_cycle and g == gaps(difference(word), base)
    assert loss <= 53*limit
    cycle = tuple(x for x, i in seen.items() if i >= start_cycle)
    assert min(map(root, cycle)) == SIGNATURE
    first = word
    for _ in range(8):
        word = literal(word, base)
    assert word == first


def physical_controls():
    assert 5*ALPHA + 4*BETA == 1
    error = Fraction(191, SCALE)
    assert ALPHA - error > Fraction(101, 1000)
    assert BETA - ALPHA - 8*error > 0
    assert Fraction(1, 9) - ALPHA - error > Fraction(1, 200)
    check_physical(100_000, (0, 0, 0, 0), (0, 0, 0, 0))
    count = 0
    for offset in range(18):
        for p in product((0, 1), repeat=4):
            u = tuple((i + offset) % 2 for i in range(4))
            check_physical(200_000 + offset, p, u)
            count += 1
    check_physical(10_000_003, (72, 86, 99, 97), (92, 45, 0, 99))
    print("literal-subtraction lifted itineraries", count + 2, "over all base residues mod 18")


if __name__ == "__main__":
    archived_witness()
    archived_counts()
    print("route, length, minimum strict sign slack, total loss, spread bound coefficient")
    results = tuple(map(check_route, ROUTES))
    assert results == ((4, 0, 177), (16, 4, 179), (52, 12, 155), (48, 18, 124), (52, 13, 72))
    tail()
    physical_controls()
