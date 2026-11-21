#!/usr/bin/env python3
"""Exact independent checks for Task 204. Usage: verify.py P9.tex P9-ancillary.zip."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path
import json
import re
import sys
from zipfile import ZipFile


def subtract(word, base):
    """Perform the column subtraction, with no gap update formulas."""
    ascending = sorted(word)
    result = []
    borrow = 0
    for high, low in zip(ascending, reversed(ascending)):
        borrow, digit = divmod(high - low + borrow, base)
        result.append(digit)
    assert borrow == 0
    return tuple(reversed(result))


def differences(word):
    ordered = sorted(word)
    return tuple(ordered[-i] - ordered[i - 1] for i in range(1, len(ordered)//2 + 1))


def gaps(word, base):
    q = (base,) + differences(word) + (0,)
    return tuple(a - b for a, b in zip(q, q[1:]))


def source_word(g):
    q = [sum(g[j:]) for j in range(1, len(g))]
    return tuple(q + [0]*len(g))


def perm(n):
    return (0,) + tuple(min(2*j, n - 2*j) for j in range(1, n//2 + 1))


def orbits(p):
    remaining = set(range(len(p)))
    output = []
    while remaining:
        start = j = min(remaining)
        orbit = []
        while j in remaining:
            orbit.append(j)
            remaining.remove(j)
            j = p[j]
        assert j == start
        output.append(tuple(orbit))
    return output


def permute(g, p):
    output = [0]*len(p)
    for j, target in enumerate(p):
        output[target] = g[j]
    return tuple(output)


def rows(g):
    output = []
    running = 0
    for i in range(len(g)//2):
        running += g[-i-1] - g[i]
        output.append(running)
    return tuple(output)


def paired(g):
    q = tuple(sum(g[i:]) for i in range(1, len(g)))
    if not q[-1]:
        return False
    base = sum(g)
    levels = [
        (x - (i == len(q) - 1), base - y - (i < len(q) - 1))
        for i, (x, y) in enumerate(zip(q, reversed(q)))
    ]
    return all(min(upper) >= max(lower) for upper, lower in zip(levels, levels[1:]))


def paired_model(g, p):
    d = rows(g)
    h = (0,) + tuple(max(0, -1 - a) for a in d)
    corrections = [0]*len(g)
    for i in range(len(d)):
        change = h[i + 1] - h[i]
        corrections[2*i:2*i + 2] = (-change, change)
    sign = (d[0] > 0) - (d[0] < 0)
    updated = list(permute(g, p))
    assert sum((a + b)**2 - a*a for a, b in zip(updated, corrections)) <= 0
    corrections[:3] = tuple(a + sign*b for a, b in zip(corrections[:3], (1, -2, 1)))
    return tuple(a + b for a, b in zip(updated, corrections))


def threshold_counts(word, base):
    return tuple(sum(d >= t for d in word) for t in range(1, base))


def repeat_digits(word, factor):
    return tuple(d for d in word for _ in range(factor))


def clock_checks(tex):
    control = 0
    for n in (7, 9, 11, 15, 17, 23, 25):
        p = perm(n)
        m = len(p)
        for a in range(-3, 4):
            for b in range(-4, 5):
                g = [256]*m
                g[1] += a
                g[-1] += b
                g[0] -= a + b
                r = len(orbits(p)[1])
                path = [tuple(g)]
                for _ in range(2*r):
                    base = sum(g)
                    assert paired(g)
                    successor = gaps(subtract(source_word(g), base), base)
                    assert paired_model(g, p) == successor
                    g = successor
                    path.append(g)
                H = max(max((max(0, -1-x) for x in rows(v)[1:]), default=0)
                        for v in path[:-1])
                D = r*(1+2*H) + 2
                mid = path[r]
                for j in range(r, 2*r):
                    assert max(abs(a-b) for a, b in zip(path[j+1], permute(path[j], p))) <= D
                rank_residual = sum(abs(mid[j]-mid[0]) for j in orbits(p)[1])
                work = sum(b*b-a*a for a, b in zip(mid, path[2*r]))
                assert work <= -2*rank_residual + 6*D*r*(r-1) + r*(8*(D-2)+6)
                control += 1
    print(f"Paired-update and two-turn work controls: {control} literal paths")

    p = perm(25)
    assert orbits(p)[2] == (5, 10) and len(orbits(p)[1]) == 10
    first = tuple(int(j == 5) for j in range(13))
    second = tuple(int(j == 10) for j in range(13))
    assert tuple(zip(rows(first)[1:], rows(second)[1:])) == ((0,0),(0,1),(0,1),(0,1),(-1,1))
    for base in range(100001, 100014):
        b = next(b for b in range(250, 263) if (base-2*b) % 13 == 0)
        a = (base-2*b)//13
        fixed = tuple(a + (b if j in (5, 10) else 0) for j in range(13))
        assert gaps(subtract(source_word(fixed), base), base) == fixed
        for copies in (1, 3):
            warm = list(fixed)
            warm[0] -= 10*copies
            for slot in orbits(p)[1]:
                warm[slot] += copies
            for _ in range(10*copies):
                assert paired(warm) and min(rows(warm)[1:]) >= -1
                warm = gaps(subtract(source_word(warm), base), base)
            assert warm == fixed
    print("25: ray on slots 5 and 10, 13 base classes checked with exact seeds and feeders")

    table = re.findall(r"^(\d+) & (\d+) & (\d+/\d+) & (\d+) \\\\", tex, re.MULTILINE)
    assert len(table) == 47
    witnesses = {}
    exceptions = {43:{3:3,7:2}, 45:{9:2,15:2}, 51:{3:2,9:3,17:0}, 57:{19:-1}}
    for length, order, peak, mass in table:
        n, r = int(length), int(order)
        p = perm(n)
        cycles = orbits(p)
        m = len(p)
        assert len(cycles[1]) == r
        weights = [0]*m
        for orbit in cycles[2:]:
            for j in orbit:
                weights[j] = exceptions.get(n, {}).get(orbit[0], 1)
        avg = Fraction(sum(weights), m)
        range_scale = avg - min(weights) if n in (45,91) else max(weights)-min(weights)
        pi = tuple(Fraction(1,m) + (w-avg)/(2*m*(m+1)*range_scale) for w in weights)
        assert permute(pi, p) == pi and min(pi) > Fraction(1,m+1)
        assert rows(pi)[0] == 0 and min(rows(pi)[1:]) > 0
        assert max(abs(sum(pi[j:])-Fraction(m-j,m)) for j in range(1,m))/2 == Fraction(peak)
        witnesses[n] = (r, pi, int(mass))
    assert len(witnesses) == 47
    print("P9 strict profiles: 47 all rational margins, correct peaks and orders")
    return witnesses


def rounded_seed(pi, base, fixed, p):
    m = len(pi)
    cycles = orbits(p)
    root = (0,) + cycles[1]
    values = [0]*m
    for j in root:
        values[j] = int(base*pi[0])
    if not fixed:
        eligible = [j for j in range(m) if j not in root]
        for j in eligible:
            values[j] = int(base*pi[j])
        values[eligible[-1]] += base-sum(values)
        return tuple(values)
    # Every complement orbit size divides r. An adjustment using a single
    # orbit has an integer solution by coprimality with r+1.
    size = len(cycles[2])
    for orbit in cycles[2:]:
        for j in orbit:
            values[j] = int(base*pi[j])
    remainder = base-sum(values)
    root_delta = (pow(len(root), -1, size)*remainder) % size if size > 1 else 0
    for j in root:
        values[j] += root_delta
    for j in cycles[2]:
        values[j] += (remainder-len(root)*root_delta)//size
    return tuple(values)


def root_step(x):
    z = x[-1]
    sign = (z > 0) - (z < 0)
    h = max(0, -1-z)
    return (z - 3*sign + 2*h, x[0] + 2*h,
            *(value - sign + h for value in x[1:-1])), h


def spread(vector):
    return max(0, max(vector))-min(0, min(vector))


def root_record(record, standard_rays=None):
    r = record["order"]
    source = tuple(record["source"])
    intervals = record["insertions"]
    start, end = record["cycle"]
    assert len(intervals) == r
    assert all(0 <= a < b <= start and (b-a) % r == 0 for a, b in intervals)
    assert all(intervals[i][1] <= intervals[i+1][0] for i in range(r-1))
    requested = {start, end, *(endpoint for interval in intervals for endpoint in interval)}
    seen = {}
    first_nonpositive = {}
    min_period = None
    x = source
    H, max_spread = 0, 0
    for t in range(end+1):
        if t in requested:
            seen[t] = x, H
        if start < t <= end and x == seen[start][0] and min_period is None:
            min_period = t-start
        max_spread = max(max_spread, spread(x))
        if t == end:
            break
        label = (r-1-t) % r
        if x[-1] <= 0:
            first_nonpositive.setdefault(label, t)
        x, h = root_step(x)
        H += h
    assert seen[start] == seen[end] and min_period is not None
    rays = []
    for a, b in intervals:
        difference = tuple(x-y for x, y in zip(seen[a][0], seen[b][0]))
        ray = tuple(difference[(j+a) % r] for j in range(r))
        assert min(ray) >= 0
        assert all(first_nonpositive.get(j, end) >= b for j in range(r) if ray[j])
        rays.append(ray)
    if standard_rays is None:
        sequence = list(first_nonpositive)
        assert len(sequence) == r
        assert all(rays[i][sequence[i]] > 0 and
                   all(rays[i][j] == 0 for j in sequence[:i]) for i in range(r))
    else:
        assert rays == standard_rays
    budget = (end+sum(b-a for a,b in intervals),
              H+sum(seen[b][1]-seen[a][1] for a,b in intervals),
              max_spread+sum(map(spread,rays)))
    return rays, min_period, budget


def retirement_checks(archive_path, witnesses):
    names = ("root-certificates.jsonl", "certificate.jsonl", "period132-root.jsonl")
    with ZipFile(archive_path) as archive:
        by_basename = {Path(name).name: name for name in archive.namelist()}
        records = [json.loads(line) for name in names
                   for line in archive.read(by_basename[name]).splitlines()]
    standard = {}
    for i, record in enumerate(records):
        shift = i == len(records)-1
        r = record["order"]
        rays, period, budget = root_record(record, standard[12][0] if shift else None)
        if not shift:
            standard[r] = rays, budget, record
        print(f"Retirement r={r}{' shifted' if shift else ''}: period {period}, "
              f"clock/loss/spread {budget}")
    assert set(standard) == set(range(3,14))
    for n, (r, pi, stored_D) in witnesses.items():
        _, budget, _ = standard[r]
        m = len(pi)
        E = m + 2*max(budget[1:]) + 1
        margin = min(min(pi)-Fraction(1,m+1), min(rows(pi)[1:])/m,
                     max(abs(sum(pi[j:])-Fraction(m-j,m)) for j in range(1,m))/m)
        D = 4*E/margin
        assert (D.numerator+D.denominator-1)//D.denominator == stored_D
    for n in (15, 17):
        p = perm(n)
        r, pi, _ = witnesses[n]
        m = len(pi)
        root = orbits(p)[1]
        complement = [i for i in range(1,m) if i not in root]
        assert len(root) == r and lcm(*(len(o) for o in orbits(p)[2:])) == (2 if n==15 else 4)
        rays, budget, record = standard[r]
        L = 3
        copies = tuple(i%L for i in range(r))
        x = tuple(a+sum(c*v[j] for c,v in zip(copies,rays))
                  for j,a in enumerate(record["source"]))
        for residue in range(m):
            base = 10**7 + (residue-10**7) % m
            for fixed_only in (False, True):
                seed = rounded_seed(pi, base, fixed_only, p)
                assert sum(seed) == base
                g = seed
                for _ in range(r):
                    assert paired(g) and min(rows(g)[1:]) >= 0
                    following = gaps(subtract(source_word(g), base), base)
                    assert following == permute(g,p)
                    g = following
                assert g == seed
            g = [int(base*v) for v in pi]
            anchor = int(base*pi[0]-Fraction(sum(x),r+1))+L-1
            for i, value in zip(root,x):
                g[i] = anchor+value
            g[0] = anchor
            free = complement[-1]
            g[free] += base-sum(g)
            seen = {}
            for t in range(budget[0]*L+1):
                g = tuple(g)
                assert paired(g) and min(rows(g)[1:]) >= -1
                if t == 0:
                    assert tuple(g[i]-g[0] for i in root) == x
                assert max(abs(Fraction(sum(g[j:]),base)-Fraction(m-j,m))
                           for j in range(1,m)) > Fraction(1,300)
                if g in seen:
                    break
                seen[g] = t
                g = gaps(subtract(source_word(g), base), base)
            else:
                raise AssertionError("lift failed to close")
    print("P9 logical budgets and 15/17 literal lifts in every base residue passed")


def prefix_checks():
    observations = 0
    max_fiber = 0
    for n in (7,9,11,13):
        for base in range(2,10):
            fibers = Counter()
            for q in combinations_with_replacement(range(base), n//2):
                source = tuple(reversed(q)) + (0,)*(n-n//2)
                output = subtract(source,base)
                target = differences(output)
                fibers[target] += 1
                h = sum(x>0 for x in q)
                if h:
                    lead = n-2*h
                    z = sorted(output,reverse=True)
                    assert z[:lead] == [base-1]*lead
                    eps = [z[lead+i]+z[-i-1]-(base-1) for i in range(h)]
                    assert max(map(abs,eps)) <= 1
                    recovered = [base-1]*lead + [None]*(2*h)
                    for i in range(h):
                        recovered[-i-1] = recovered[i]-target[i]
                        recovered[lead+i] = base-1+eps[i]-recovered[-i-1]
                    assert recovered == z
                observations += 1
            max_fiber = max(max_fiber,max(fibers.values()))
    print(f"Fixed-prefix inversion including shorter heights: {observations} source states, max fiber {max_fiber}")


def breach_checks():
    for n, representatives in ((229,(3,5)),(283,(3,5)),(1069,(3,7))):
        p = perm(n)
        cycles = orbits(p)
        assert len(cycles) == 4
        complement = tuple(next(o for o in cycles if a in o) for a in representatives)
        assert set(complement) == set(cycles[2:])
        indicators = [tuple(int(j in orbit) for j in range(len(p))) for orbit in complement]
        rs = tuple(zip(*(rows(v)[1:] for v in indicators)))
        det = lambda a,b: a[0]*b[1]-a[1]*b[0]
        certificate = None
        for tri in combinations(range(len(rs)),3):
            a,b,c = (rs[i] for i in tri)
            weights = (det(b,c),det(c,a),det(a,b))
            if max(weights)<0:
                weights = tuple(-w for w in weights)
            if min(weights)>0:
                common = gcd(*weights)
                weights = tuple(w//common for w in weights)
                assert tuple(sum(w*v[j] for w,v in zip(weights,(a,b,c)))
                             for j in range(2)) == (0,0)
                certificate = (tuple(i+2 for i in tri),(a,b,c),weights)
                break
        assert certificate is not None
        print(f"Zero-cone certificate n={n}, r={len(cycles[1])}: {certificate}")
    for scale in (7, 101, 1001):
        base = 71*scale
        dev = (-scale,-2*scale,-3*scale,-4*scale+12,scale-18,-18,
               -scale-6,-2*scale+6,-3*scale+18,2*scale-12,scale)
        numerators = [(12-i)*base+a for i,a in enumerate(dev,1)]
        assert all(num%12 == 0 for num in numerators)
        xs = tuple(num//12 for num in numerators)+(0,)*12
        g = gaps(xs,base)
        assert paired(g) and gaps(subtract(xs,base),base) == g
        assert max(max(0,-1-d) for d in rows(g)[1:]) == max(2,(scale-5)//2)
    print("23-digit full-orbit control: paired exterior fixed points, linear breach")


def replication_checks():
    scanned = commuting = 0
    for n in range(2,9):
        for base in range(2,8):
            for word in combinations_with_replacement(range(base),n):
                q = differences(word)
                h = sum(x>0 for x in q)
                predicted = h == 0 or q[0]+q[h-1] == base
                small = threshold_counts(subtract(word,base),base)
                for factor in (2,3):
                    large = threshold_counts(subtract(repeat_digits(word,factor),base),base)
                    diff = tuple(a-factor*b for a,b in zip(large,small))
                    expected = tuple((1-factor)*(int(t==base-q[0])-int(t==q[h-1]))
                                     if h else 0 for t in range(1,base))
                    assert diff == expected and (not any(diff)) == predicted
                scanned += 1
                commuting += predicted
    print(f"Replication factors 2 and 3: {scanned} multisets each, {commuting} commuting")

    checked_cycles = 0
    for n, base in ((3,4),(5,4),(5,5),(7,3),(7,5),(8,4),(9,4),(10,4)):
        seen = set()
        for start in combinations_with_replacement(range(base),n):
            orbit = {}
            current = start
            while current not in orbit and current not in seen:
                orbit[current] = len(orbit)
                current = tuple(sorted(subtract(current,base)))
            seen.update(orbit)
            if current not in orbit:
                continue
            cycle = list(orbit)[orbit[current]:]
            if len(cycle) <= 1:
                continue
            ranks = [tuple(q for q in differences(phase) if q) for phase in cycle]
            if not all(q and q[0]+q[-1] == base for q in ranks):
                continue
            lifted = [tuple(sorted(repeat_digits(phase,3))) for phase in cycle]
            words = [subtract(phase,base) for phase in lifted]
            assert len(set(words)) == len(cycle)
            assert all(subtract(words[i-1],base) == word for i,word in enumerate(words))
            checked_cycles += 1
    assert checked_cycles > 0
    print(f"Replicated nonfixed numerical cycles checked phasewise: {checked_cycles}")

    for N in (9,15,21,25,27,35,45,49,81,105):
        divisors = [d for d in range(3,N,2) if N%d == 0]
        m = (N+1)//2
        e = divisors[0]
        t = sum(j%e == 0 for j in range(1,m))
        v = tuple(int(j>0 and j%e==0) for j in range(m))
        assert rows(v)[1:] == tuple(int((i-1)%e >= (e-1)//2) for i in range(2,m//2+1))
        base = lcm(m, *((d+1)//2 for d in divisors), 2*m*m*(m+1))*2
        words = []
        for n in divisors+[N]:
            d = (n+1)//2
            a = base//d
            central = tuple(x for j in range(1,d) for x in (j*a-1,j*a))+(base-1,)
            word = subtract(repeat_digits(central,N//n),base)
            assert subtract(word,base) == word
            assert differences(word).count(0) == (N//n-1)//2
            words.append(word)
        pi = tuple(Fraction(1,m)+(Fraction(value)-Fraction(t,m))/(2*m*(m+1)) for value in v)
        assert all((base*u).denominator == 1 for u in pi)
        g = tuple(int(base*u) for u in pi)
        words.append(subtract(source_word(g),base))
        assert subtract(words[-1],base) == words[-1] and paired(g)
        assert len(set(words)) == len(divisors)+2
    print("Odd composite same-base divisor coexistence: 10 literal fixtures")


def main():
    tex = Path(sys.argv[1]).read_text()
    witnesses = clock_checks(tex)
    retirement_checks(sys.argv[2],witnesses)
    prefix_checks()
    breach_checks()
    replication_checks()


if __name__ == "__main__":
    main()
