"""Independent checks of equations (2), (7)–(9), and (37), using the Python standard library."""
from collections import Counter, deque
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def cross(p, q, r):
    return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])


def vertices(forms):
    """Find the convex feasible polygon for a*x+b*y+c >= 0 by line intersections."""
    points = set()
    for (a,b,c), (d,e,f) in combinations(forms, 2):
        det = a*e-b*d
        if det:
            x, y = (b*f-c*e)/det, (c*d-a*f)/det
            if all(i*x+j*y+k >= 0 for i,j,k in forms):
                points.add((x,y))
    def half(ordered):
        part = []
        for p in ordered:
            while len(part) > 1 and cross(part[-2], part[-1], p) <= 0:
                part.pop()
            part.append(p)
        return part[:-1]
    ordered = sorted(points)
    return half(ordered) + half(reversed(ordered))


def integral(poly):
    """Four-point triangle rule, exact through degree three (the density is cubic)."""
    total = Q(0)
    for i in range(1, len(poly)-1):
        triangle = [poly[0], poly[i], poly[i+1]]
        def density(weights):
            x, y = (sum(w * p[j] for w, p in zip(weights, triangle)) for j in (0,1))
            return 120*(1-x)*y*(x-y)
        value = -Q(9,16)*density([Q(1,3)]*3)
        for j in range(3):
            value += Q(25,48)*density([Q(3 if k == j else 1,5) for k in range(3)])
        total += cross(*triangle)*value/2
    return total


def cdf():
    # Inside the event, |A|, |C| <= r, hence U=(2+A+C)/3 <= 53/75.
    # V <= U, so (V-1/2)/(1-U) < 1. Only shear counts zero and one can contribute.
    r, half = Q(3,50), Q(1,2)
    def guards(x, y, cell, n):
        regions = [(half-x,), (x-half, 1-x-y), (x+y-1, half-y), (y-half,)]
        U, V = [(1-y,1-x-y),(1-y,x-y),(x,x-y),(x,x+y-1)][cell]
        w = V - n*(1-U)
        A, C = U+w-1, 2*U-w-1
        shear = (half-V,) if n == 0 else (V-half,half-w)
        return (y, x-y, 1-x, *regions[cell], *shear, r-A, r+A, r-C, r+C, r-A+C, r+A-C)
    pieces = []
    for cell in range(4):
        for n in (0,1):
            evaluations = [guards(Q(x),Q(y),cell,n) for x,y in [(0,0),(1,0),(0,1)]]
            forms = [(i-c,j-c,c) for c,i,j in zip(*evaluations)]
            pieces.append(integral(vertices(forms)))
    value = sum(pieces)
    require(value == Q(713903,39062500) and 1-value == Q(38348597,39062500), 'Equation (2)')
    return {'cdf': str(value), 'tail': str(1-value), 'cell_masses': list(map(str,pieces))}


def index(u, v):
    return u*(u+1)//2 - 1 + v


def successor(B, u, v):
    # A representative is the descending word (u,v,0,0,0). Subtract its reversal
    # as positional integers, then unpack exactly five digits, including leading zeros.
    word = (u,v,0,0,0)
    desc = asc = 0
    for a,b in zip(word, reversed(word)):
        desc, asc = desc*B+a, asc*B+b
    value = desc-asc
    digits = []
    for _ in word:
        value, d = divmod(value,B)
        digits.append(d)
    digits.sort()
    return index(digits[4]-digits[0], digits[3]-digits[1])


def weight(B,u,v):
    if v == 0:
        return (B-u)*(20*u-10)
    if u == v:
        return (B-u)*(30*u-10)
    return (B-u)*(120*v*(u-v)-20)


def graph(B):
    pairs = [(u,v) for u in range(1,B) for v in range(u+1)]
    nxt = [successor(B,*p) for p in pairs]
    degree = [0]*len(nxt)
    for j in nxt:
        degree[j] += 1
    leaves = deque(i for i,d in enumerate(degree) if d == 0)
    peeled = []
    while leaves:
        i = leaves.popleft()
        peeled.append(i)
        j = nxt[i]
        degree[j] -= 1
        if degree[j] == 0:
            leaves.append(j)
    h = [0]*len(nxt)
    two = [nxt[i] != i and nxt[nxt[i]] == i for i in range(len(nxt))]
    for i in reversed(peeled):
        h[i] = h[nxt[i]]+1
        two[i] = two[nxt[i]]
    # One numerical cycle word per cyclic pair. All other words with this pair
    # have clock h+1. Every pair has weight > 1 in this range.
    mass = [weight(B,*p) for p in pairs]
    require(min(mass)>1 and sum(mass)==B**5-B, f'Input counts at {B}')
    maximum = max(h)+1
    stop = (B-1)//2
    n = (1+isqrt(1+8*(B-3-stop)))//2
    witness = h[index(B-1,stop+n*(n-1)//2)]+1
    ceil_sqrt = 1+isqrt(B-1)
    require(3*(witness+5)-B >= ceil_sqrt and witness <= maximum
            and 3*(maximum-60)-B <= isqrt(B), f'Maximum bounds at {B}')
    hist = Counter()
    basin_pairs = []
    for i,p in enumerate(pairs):
        if two[i]:
            cycle_word = int(h[i] == 0)
            hist[h[i]+1] += mass[i]-cycle_word
            hist[0] += cycle_word
            basin_pairs.append(p)
    if B%2 == 0:
        m = B//2
        require(set(basin_pairs)=={(m,m-1),(m+1,0),(m,0)}, f'Two-cycle basin at {B}')
        require(hist=={0:2,1:35*B*B-75*B-12,2:5*B*B-5*B}
                and sum(hist.values())==40*B*B-80*B-10, f'Two-cycle histogram at {B}')
    else:
        require(not basin_pairs, f'Unexpected two-cycle at {B}')
    return {'base':B,'maximum':maximum,'witness_time':witness,'two_cycle_histogram':dict(hist)}


if __name__ == '__main__':
    values = cdf()
    rows = [graph(b) for b in range(30,201)]
    out = Path(__file__).resolve().parent / 'results'
    out.mkdir(exist_ok=True)
    (out/'headline_results.json').write_text(json.dumps({**values,'bases':rows},indent=2)+'\n')
    print(f"F(1/25) = {values['cdf']}, tail = {values['tail']}.")
    print('All 171 weighted graphs passed, including 86 even-base basins and all 171 witnesses.')
