"""Verify the prescribed cohort routes on the whole box using its 32 corners.

Usage in a private snapshot: python3 -B cohort_replay.py reconstruction.txt
The separately supplied itinerary is data. No discovery code is imported.
"""
import re
import sys
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from math import lcm, prod

from audit import CENTER, EPSILON, basis_form, centered, dot, trace


def parse(text):
    routes = [[], []]
    arm = 0
    for line in text.splitlines():
        if line == 'ODD':
            arm = 1
        elif line.startswith('step '):
            routes[arm].append(int(line.split()[1]))
        elif line.startswith('repeat '):
            match = re.fullmatch(
                r'repeat \(([^)]*)\) at \d+ delta \(([^)]*)\) numerator \(([^)]*)\)'
                r' divisor (\d+) remainder (\d+) count \d+', line)
            assert match
            w, d, u = match.groups()[:3]
            routes[arm].append((tuple(map(int,w.split(','))), tuple(map(int,d.split(','))),
                                tuple(map(Q,u.split(','))), int(match[4]), int(match[5])))
    assert all(sum(isinstance(t,tuple) for t in route) == 6 for route in routes)
    return routes


def translate(x, delta, n):
    return tuple(u+n*d for u,d in zip(x,delta))


def word_step(x, word, slacks):
    for code in word:
        digits = [-6]
        for i,u in enumerate(x):
            a, b = u-6*(i == 4), -x[4-i]-6*(i < 4)
            if code & (1 << i):
                digits.extend((a,b))
                slacks.append(a-b)
            else:
                digits.extend((b,a))
                slacks.append(b-a-6)
        x = tuple(digits[i]-digits[-1-i] for i in range(5))
    return x


def replay(route, start, arm):
    x = translate(start,(1,2,3,4,-1),arm)
    elapsed = 0
    slacks = []
    for token in route:
        if isinstance(token,int):
            x = word_step(x,(token,),slacks)
            elapsed += 1
            continue
        word, delta, row, modulus, residue = token
        count = 1+(dot(row,(*start,1))-residue)/modulus
        slacks.append(count-1)
        assert word_step(x,word,slacks) == translate(x,delta,1)
        assert word_step(translate(x,delta,1),word,[]) == translate(x,delta,2)
        last = translate(x,delta,count-1)
        x = translate(x,delta,count)
        assert word_step(last,word,slacks) == x
        elapsed += len(word)*count
    slacks.append(count-30)
    endpoint = ((12,12,36,42,18),(31,50,75,82,47))[arm]
    assert x == endpoint
    spot = translate(endpoint,delta,-30)
    mu, period = trace(spot, 271140+arm)
    assert (mu,period) == ((658,12),(619,24))[arm]
    clock = elapsed - 30*len(word)+mu
    expected = (('3/16','55/48','-5/8','-55/48','95/48','39'),
                ('205/1152','715/1152','-65/192','-715/1152','1235/1152','209/24'))[arm]
    assert clock == dot(tuple(map(Q,expected)),(*start,1))
    return slacks


def gap_index(routes):
    rows, moduli = [], []
    seed = centered(CENTER)
    for route in routes:
        for token in route:
            if isinstance(token,int):
                continue
            _,_,form,m,r = token
            assert (dot(form,(*seed,1))-r) % m == 0
            slope = [-6*sum(form[j:5])/m for j in range(5)]
            d = lcm(*(u.denominator for u in slope))
            moduli.append(d)
            rows.append(tuple(int(u*d) for u in slope))
    diagonal = [[d*int(i == j) for j in range(len(moduli))]
                for i,d in enumerate(moduli)]
    h = basis_form((*diagonal,*zip(*rows)))
    return Q(prod(moduli),prod(h[i][i] for i in range(len(moduli))))


if __name__ == '__main__':
    routes = parse(Path(sys.argv[1]).read_text())
    assert gap_index(routes) == 3317760
    for arm, route in enumerate(routes):
        for signs in product((-1,1),repeat=5):
            y = tuple(Q(c,180000)+s*EPSILON for s,c in zip(signs,CENTER))
            unit = centered((*y,1-sum(y)))
            low, next_ = (replay(route,tuple(b*u for u in unit),arm) for b in (10**6,10**6+1))
            assert len(low) == len(next_)
            assert all(0 <= a <= b for a,b in zip(low,next_))
        print(f'arm {arm}: 32 corners, {len(low)} affine guard/count inequalities at the cutoff and on its outgoing ray')
    print('all route identities, both tail clocks, 12 seed residues and joint index 3317760 pass')
