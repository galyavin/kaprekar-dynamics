"""Exact predecessor reconstruction for eleven-digit difference states."""

import itertools
from fractions import Fraction

def digits_of(q, B):
    """digits (unsorted, low->high) of desc-asc computed by position-wise borrow."""
    v = [-q[j] for j in range(5)] + [0] + [q[10 - j] for j in range(6, 11)]
    d = []; carry = 0
    for j in range(11):
        x = v[j] + carry
        if x < 0: d.append(x + B); carry = -1
        else: d.append(x); carry = 0
    assert carry == 0
    return d

def F(q, B):
    d = sorted(digits_of(q, B), reverse=True)
    return tuple(d[i] - d[10 - i] for i in range(5))

def Fint(q, B):  # independent big-int version
    L = sum(q[i-1] * (B**(11 - i) - B**(i - 1)) for i in range(1, 6))
    d = []
    for _ in range(11): d.append(L % B); L //= B
    assert L == 0
    d.sort(reverse=True)
    return tuple(d[i] - d[10 - i] for i in range(5))

def gaps(q, B): return (B-q[0], q[0]-q[1], q[1]-q[2], q[2]-q[3], q[3]-q[4], q[4])

def matchings(items):
    if not items: yield []; return
    a = items[0]
    seen = set()
    for j in range(1, len(items)):
        if items[j] in seen: continue
        seen.add(items[j])
        rest = items[1:j] + items[j+1:]
        for m in matchings(rest): yield [(a, items[j])] + m

def q_from_core(core, h, B):
    """all q (h positive entries) whose pair-multiset equals core (multiset of 2h digits)"""
    C = B - 1; res = set()
    items = sorted(core, reverse=True)
    if h == 1:
        for x in set(items):
            q1 = x + 1
            res.add((q1,0,0,0,0))
        return res
    for mt in matchings(items):
        sums = [a+b for a, b in mt]
        if sorted(sums) != sorted([B, C-1] + [C]*(h-2)): continue
        pB = mt[sums.index(B)]; pm = mt[sums.index(C-1)]
        mids = [p for p, s in zip(mt, sums) if s == C]
        for q1 in set(pB):
            for qh0 in set(pm):
                for choice in itertools.product(*[set(p) for p in mids]):
                    qs = sorted([q1, qh0 + 1] + list(choice), reverse=True)
                    res.add(tuple(qs) + (0,)*(5-h))
    return res


_cache = {}
def _system(h):
    """per-h: matrix A (rows: h mirror, 5 target), reduced to solution map R (n x m) and consistency map Z (k x m)"""
    if h in _cache: return _cache[h]
    n = 2*h; nC = 11 - n
    A = []; Bc = []   # Bc: for each row, coefficient dict over [C, r_0..r_{h-1}, t_1..t_5]  (b = vector in R^{1+h+5})
    for i in range(h):
        row = [0]*n; row[i] = 1; row[n-1-i] += 1; A.append(row)
        v = [0]*(1+h+5); v[0] = 1; v[1+i] = 1; Bc.append(v)
    for i in range(1, 6):
        row = [0]*n; v = [0]*(1+h+5); v[1+h+i-1] = 1
        for pos, sgn in ((i, 1), (12-i, -1)):
            if pos <= nC: v[0] -= sgn
            else: row[pos-nC-1] += sgn
        A.append(row); Bc.append(v)
    M = [[Fraction(x) for x in a] + [Fraction(x) for x in b] for a, b in zip(A, Bc)]
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: raise RuntimeError("underdetermined h=%d" % h)
        M[r], M[p] = M[p], M[r]; M[r] = [x/M[r][c] for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f*y for x, y in zip(M[i], M[r])]
        r += 1
    R = [row[n:] for row in M[:n]]; Z = [row[n:] for row in M[n:]]
    _cache[h] = (R, Z); return R, Z

def preds(t, B, rmax=2):
    C = B - 1; out = set()
    for h in range(1, 6):
        R, Z = _system(h); n = 2*h
        for r in itertools.product(range(-rmax, rmax+1), repeat=h):
            if sum(r) != 0: continue
            b = [C] + list(r) + list(t)
            if any(sum(z*x for z, x in zip(row, b)) != 0 for row in Z): continue
            sol = [sum(z*x for z, x in zip(row, b)) for row in R]
            if any(x.denominator != 1 for x in sol): continue
            d = [int(x) for x in sol]
            if any(d[i] < d[i+1] for i in range(n-1)) or d[0] > C or d[-1] < 0: continue
            for q in q_from_core(d, h, B):
                if q[0] <= C and q[h-1] > 0 and F(q, B) == t:
                    assert Fint(q, B) == t
                    out.add(q)
    return out
