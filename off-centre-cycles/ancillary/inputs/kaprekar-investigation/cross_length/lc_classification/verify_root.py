"""Replay exact root-order certificates and check their affine lifts."""
import argparse
from fractions import Fraction
import json
from math import ceil, lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIFT_CONTROLS = {
    31: {},
    33: {},
    41: {},
    43: {3: 3, 7: 2},
    45: {3: 1, 5: 1, 9: 2, 15: 2},
    51: {3: 2, 5: 1, 9: 3, 17: 0},
    57: {3: 1, 5: 1, 19: -1},
    63: {},
    65: {},
}


def branch(z):
    return (z > 0) - (z < 0)


def step(x):
    z = x[-1]
    e, h = branch(z), max(0, -1 - z)
    return (z - 3*e + 2*h, x[0] + 2*h, *(q - e + h for q in x[1:-1])), h


def rotate(x, shift):
    shift %= len(x)
    return x[-shift:] + x[:-shift]


def spread(x):
    return max(0, *x) - min(0, *x)


def determinant(columns):
    rows = [[Fraction(x) for x in row] for row in zip(*columns)]
    determinant = Fraction(1)
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            determinant *= -1
        factor = rows[j][j]
        determinant *= factor
        for i in range(j + 1, len(rows)):
            ratio = rows[i][j] / factor
            rows[i][j:] = [a - ratio*b for a, b in zip(rows[i][j:], rows[j][j:])]
    return determinant


def verify(record):
    r = record["order"]
    source = tuple(record["source"])
    capture, stop = record["cycle"]
    assert len(source) == r and 0 <= capture < stop
    path, loss = [source], [0]
    for _ in range(stop):
        x, h = step(path[-1])
        path.append(x)
        loss.append(loss[-1] + h)
    assert path[capture] == path[stop] and loss[capture] == loss[stop]
    cycle = path[capture:stop]
    period = cycle[1:].index(cycle[0]) + 1 if cycle[0] in cycle[1:] else len(cycle)
    assert len(set(cycle[:period])) == period and all(x == cycle[i % period] for i, x in enumerate(cycle))
    assert all(z >= -(r + 1) for x in path[r:] for z in x)

    rays, words = [], []
    last = 0
    for start, end in record["insertions"]:
        assert last <= start < end <= capture and (end - start) % r == 0
        last = end
        ray = tuple(a-b for a,b in zip(path[start], path[end]))
        source_ray = rotate(ray, -start)
        # Nonzero slope is allowed only at a positive comparison.
        for t in range(end):
            if source_ray[(r-1-t) % r]:
                assert path[t][-1] > 0, (r, start, end, t)
        rays.append(source_ray)
        words.append((end - start, loss[end] - loss[start], spread(ray)))
    assert len(rays) == r and (det := abs(determinant(rays))) > 0

    # A larger lower coefficient on the last ray acting on each coordinate
    # protects the nonnegative aggregate perturbation until its own word.
    pivots = [max(i for i in range(r) if rays[i][j]) for j in range(r)]
    floors = []
    for i, ray in enumerate(rays):
        floor = 0
        for j, last_word in enumerate(pivots):
            if i != last_word:
                continue
            assert ray[j] > 0
            deficit = sum((floors[k]+1)*max(0,-rays[k][j]) for k in range(i))
            floor = max(floor, (deficit + ray[j] - 1)//ray[j])
        floors.append(floor)
    time_bound, loss_bound, spread_bound = [base + sum((floor+1)*word[j] for floor, word in zip(floors, words))
                                           for j, base in enumerate((stop, loss[-1], max(map(spread,path))))]
    print(f"r={r} determinant={det} period={period} last={stop} "
          f"time/L<={time_bound} loss/L<={loss_bound} spread/L<={spread_bound} coefficient_floors={floors}")
    return source, rays, tuple(cycle[:period]), max(loss_bound, spread_bound)


def literal(q, base):
    word = (*q, *(0,)*(len(q) + 1))
    desc = asc = 0
    for a, b in zip(word, reversed(word)):
        desc, asc = desc*base + a, asc*base + b
    num = desc - asc
    digits = []
    for _ in word:
        num, digit = divmod(num, base)
        digits.append(digit)
    assert num == 0
    digits.sort(reverse=True)
    return tuple(digits[j] - digits[-1-j] for j in range(len(q)))


def roots(n):
    k = n//2
    order = [1]
    while (new := min(2*order[-1], n - 2*order[-1])) != 1:
        order.append(new)
    assert order[-1] == k
    return (0, *order)


def gap_step(g, n):
    z = g[-1] - g[0]
    e, h = branch(z), max(0, -1 - z)
    nxt = [g[0], *(0,)*len(g[1:])]
    for j in range(1, len(g)):
        nxt[min(2*j, n-2*j)] = g[j]
    nxt[0] += e-h
    nxt[1] += -2*e+h
    nxt[2] += e+h
    nxt[3] -= h
    return tuple(nxt)


def physical_check(n, certificate):
    m = n//2 + 1
    j = roots(n)
    k = [x for x in range(m) if x not in j]
    source, rays, cycle, budget = certificate
    assert len(j) == len(source) + 1
    x = tuple(a + b for a,b in zip(source, rays[-1]))
    epsilon = Fraction(1, 4*m*(m+1))
    profile = [Fraction(1,m) - epsilon*len(k)/m + epsilon*(i in k) for i in range(m)]
    guards = [sum(profile[-i:]) - sum(profile[:i]) for i in range(2, m//2 + 1)]
    assert min(guards) > 0
    b = lcm(*(a.denominator for a in profile))
    base = b * (100 + 16*m*m*(budget+1))
    g = [int(a*base) for a in profile]
    for d, idx in zip(x, j[1:]):
        g[idx] += d
    g[k[-1]] -= sum(x)
    path = set()
    g = tuple(g)
    while g not in path:
        assert (m+1)*min(g) >= base+1
        assert all(sum(g[-i:]) >= sum(g[:i])-1 for i in range(2, m//2 + 1))
        path.add(g)
        q = tuple(sum(g[i:]) for i in range(1,m))
        nxt = gap_step(g, n)
        assert literal(q, base) == tuple(sum(nxt[i:]) for i in range(1,m))
        g = nxt
    assert tuple(g[i]-g[0] for i in j[1:]) in cycle
    print(f"matched n={n}, B={base}, {len(path)} literal whole-word updates")


def exact_lift(n, certificate):
    m = n//2+1
    j = roots(n)
    unvisited = set(range(m)) - set(j)
    weights = [0]*m
    while unvisited:
        orbit_rep = min(unvisited)
        i = orbit_rep
        while i in unvisited:
            unvisited.remove(i)
            weights[i] = LIFT_CONTROLS[n].get(orbit_rep, 1)
            i = min(2*i, n-2*i)
        assert i == orbit_rep
    center = Fraction(sum(weights), m)
    epsilon = Fraction(1,2*m*(m+1)*(max(weights)-min(weights)))
    pi = tuple(Fraction(1,m)+epsilon*(w-center) for w in weights)
    assert sum(pi) == 1 and all(pi[i] == pi[min(2*i, n-2*i)] for i in range(m))
    spare = min(pi) - Fraction(1,m+1)
    margin = min(sum(pi[-i:])-sum(pi[:i]) for i in range(2,m//2+1))
    amplitude = max(abs(sum(pi[i:])-Fraction(m-i,m)) for i in range(1,m))
    room = min(spare, margin/m, amplitude/m)
    assert room > 0
    error = m+2*certificate[3]+1
    denominator = ceil(4*error/room)
    assert Fraction(error, denominator) <= room/4
    print(f"strict n={n} order={len(j)-1} delta={amplitude/2} D={denominator} "
          f"eta={n}!/({2*(m+1)}^{m}*{denominator}^{m-1})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--literal", action="store_true")
    args = parser.parse_args()
    for line in (HERE / "root-certificates.jsonl").read_text().splitlines():
        record = json.loads(line)
        certificate = verify(record)
        for n in LIFT_CONTROLS:
            if len(roots(n)) == record["order"]+1:
                exact_lift(n, certificate)
        if args.literal and record["order"] in (5, 6):
            for n in (2**record["order"] - 1, 2**record["order"] + 1):
                physical_check(n, certificate)


if __name__ == "__main__":
    main()
