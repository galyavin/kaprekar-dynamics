"""Exact checks of the unpaired periodic basin certificates.

Run from the repository root with a Python interpreter containing SymPy.
Branch geometry and Fourier--Motzkin cone elimination are regenerated here.
The four earlier normal-map proofs and paired-cell proofs are separate inputs.
"""

import json, gzip, itertools, functools, math, collections
import sympy as sp
from fractions import Fraction as Q
from pathlib import Path

def atom(row):
    d = math.gcd(*row)
    return tuple((x // d for x in row)) if d else tuple(row)

def neg(r):
    return tuple((-x for x in r))

def sub(r, v):
    return tuple((x - y for (x, y) in zip(r, v)))

def add(r, v):
    return tuple((x + y for (x, y) in zip(r, v)))

def tran(u, M):
    return tuple((sum((x * y for (x, y) in zip(u, c))) for c in zip(*M)))

def matrix(A, B):
    return tuple((tran(a, B) for a in A))

def matvec(A, v):
    return tuple((sum((x * y for (x, y) in zip(a, v))) for a in A))
I = tuple((tuple((int(i == j) for j in range(5))) for i in range(5)))

def can(u):
    return u if next((x for x in u if x)) > 0 else neg(u)

def strict(rr):
    return feasible(tuple(sorted(set((atom(r) for r in rr)))))

@functools.lru_cache(maxsize=None)
def feasible(rows):
    """Eliminate the variable with the fewest positive/negative pairs from a strict homogeneous system. Opposite signs combine with positive multipliers, and a strict zero row is impossible."""
    if not rows:
        return True
    n = len(rows[0])
    if not n or (0,) * n in rows:
        return False
    S = set(rows)
    if any((neg(r) in S for r in S)):
        return False
    i = min(range(n), key=lambda j: sum((r[j] > 0 for r in rows)) * sum((r[j] < 0 for r in rows)))
    plus = [r for r in rows if r[i] > 0]
    minus = [r for r in rows if r[i] < 0]
    rem = {r[:i] + r[i + 1:] for r in rows if not r[i]}
    for p in plus:
        for m in minus:
            v = tuple((p[i] * y - m[i] * x for (x, y) in zip(p, m)))
            rem.add(atom(v[:i] + v[i + 1:]))
    return feasible(tuple(sorted(rem)))

def family(row):
    """Regenerate all adjacent signed sorting branches from the rational reference cycle. Compute phase tangents from pulled-back tie equations and the derivative of the periodicity equation."""
    q = tuple(map(Q, row['witness']))
    p = row['period']
    A = I
    eq = []
    phases = []
    prefix = []
    points = []
    for step in range(p):
        points.append(q)
        prefix.append(A)
        mag = [abs(x - Q(1, 2)) for x in q]
        levels = sorted(set(mag), reverse=True)
        groups = [tuple((i for i in range(5) if mag[i] == u)) for u in levels]
        signs = tuple((1 if x > Q(1, 2) else -1 if x < Q(1, 2) else 0 for x in q))
        branches = []
        for perm in itertools.product(*(itertools.permutations(g) for g in groups)):
            order = sum(perm, ())
            for orient in itertools.product(*((signs[i],) if signs[i] else (-1, 1) for i in order)):
                r = {i: tuple((s * x for x in I[i])) for (i, s) in zip(order, orient)}
                guards = tuple((sub(r[i], r[j]) for g in perm for (i, j) in zip(g, g[1:])))
                if levels[-1] == 0:
                    guards += (r[order[-1]],)
                u = tuple((r[i] for i in order))
                H = (u[0],) + tuple((add(a, b) for (a, b) in zip(u, u[1:])))
                branches.append((H, guards))
        phases.append(branches)
        (H, L) = branches[0]
        for g in L:
            eq.append(tran(g, A))
        A = matrix(H, A)
        v = [abs(x - Q(1, 2)) for x in q]
        v.sort(reverse=True)
        q = (Q(1, 2) + v[0],) + tuple((a + b for (a, b) in zip(v, v[1:])))
    assert q == points[0]
    eq.extend((sub(r, s) for (r, s) in zip(A, I)))
    tangent = [tuple(map(Q, v)) for v in sp.Matrix(eq).nullspace()]
    assert len(tangent) == row['dimension']
    return (phases, [tuple((matvec(A, t) for t in tangent)) for A in prefix], points)

def check(row, data):
    """Check every backward cone transition, including feasibility, orientation up to sign, and inclusion in a stored state. A reversed strict guard detects failed inclusion because its source cone is open."""
    (ph, T, ps) = family(row)
    cover = collections.defaultdict(list)
    for (j, v, C) in data['cover']:
        v = tuple(v)
        C = tuple(map(tuple, C))
        assert any(v) and v == can(v)
        assert all((not sum((a * b for (a, b) in zip(v, t))) for t in T[j]))
        assert strict(C)
        cover[j, v].append(C)
    assert () in cover[0, tuple(data['seed'])]
    for ((j, v), Cs) in list(cover.items()):
        h = (j - 1) % len(ph)
        for C in Cs:
            for (A, L) in ph[h]:
                H = L + tuple((tran(g, A) for g in C))
                if not strict(H):
                    continue
                w = can(tran(v, A))
                other = cover.get((h, w), [])
                assert any((all((g in H or not strict(H + (neg(g),)) for g in D)) for D in other))

def mixed(rr, ww):
    rows = {atom(r): False for r in ww}
    rows.update(((atom(r), True) for r in rr))
    return mixed_elim(tuple(sorted(rows.items())))

@functools.lru_cache(maxsize=None)
def mixed_elim(rows):
    """Eliminate variables from a cone of strict and weak homogeneous inequalities. Tag each derived row strict if at least one summand is strict."""
    if not rows:
        return True
    n = len(rows[0][0])
    S = dict(rows)
    if S.get((0,) * n):
        return False
    S.pop((0,) * n, None)
    if not S:
        return True
    if any(((s or S.get(neg(r), False)) and neg(r) in S for (r, s) in S.items())):
        return False
    i = min(range(n), key=lambda j: sum((r[j] > 0 for r in S)) * sum((r[j] < 0 for r in S)))
    rem = {r[:i] + r[i + 1:]: s for (r, s) in S.items() if not r[i]}
    for (p, sp_) in S.items():
        if p[i] <= 0:
            continue
        for (m, sm) in S.items():
            if m[i] >= 0:
                continue
            v = tuple((p[i] * y - m[i] * x for (x, y) in zip(p, m)))
            v = atom(v[:i] + v[i + 1:])
            rem[v] = rem.get(v, False) or sp_ or sm
    return mixed_elim(tuple(sorted(rem.items())))
G = (neg(I[0]),) + tuple((sub(a, b) for (a, b) in zip(I, I[1:]))) + (I[-1],)

def gaps(q):
    a = (1,) + tuple(q) + (0,)
    return sub(a[:-1], a[1:])

def dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)))

def order_profiles(ps, T):
    """Enumerate endpoints and open intervals of the within-cell changes of gap order. Sorting-cell endpoints themselves belong to other catalogue rows."""
    if not T[0]:
        return [ps]
    assert len(T[0]) == 1
    cuts = []
    guards = []
    for (q, (t,)) in zip(ps, T):
        gg = gaps(q)
        vv = matvec(G, t)
        guards.extend(zip(gg, vv))
        for (i, j) in itertools.combinations(range(6), 2):
            if vv[i] != vv[j]:
                cuts.append((gg[j] - gg[i]) / (vv[i] - vv[j]))
        mags = []
        for (x, v) in zip(q, t):
            a = x - Q(1, 2)
            if not a:
                assert not v
            mag = (abs(a), v if a > 0 else -v)
            if a:
                guards.append(mag)
            mags.append(mag)
        mags.sort(reverse=True)
        for (x, y) in zip(mags, mags[1:]):
            (a, b) = sub(x, y)
            if a:
                guards.append((a, b))
            else:
                assert not b
    lo = max((-a / b for (a, b) in guards if b > 0))
    hi = min((-a / b for (a, b) in guards if b < 0))
    knots = sorted({lo, hi, *(c for c in cuts if lo < c < hi)})
    pts = knots[1:-1] + [(a + b) / 2 for (a, b) in zip(knots, knots[1:])]
    return [[tuple((a + x * b for (a, b) in zip(q, t[0]))) for (q, t) in zip(ps, T)] for x in pts]

def major(row, d):
    """Enforce all limiting majorization, anchor, and sorting inequalities on the same finite prefix. Every branch must be empty by the recorded depth."""
    (ph, T, ps) = family(row)
    profiles = order_profiles(ps, T)
    depth = max((len(p['counts']) - 1 for p in d['profiles']))
    assert len(d['profiles']) == len(profiles)
    for zs in profiles:
        rules = []
        for (z, vs) in zip(zs, T):
            lines = [neg(I[0])]
            gs = gaps(z)
            for v in sorted(set(gs))[:-1]:
                lines.append(neg(tuple((sum((g[i] for (g, h) in zip(G, gs) if h <= v)) for i in range(5)))))
            rules.append((tuple((r for r in lines if all((not dot(r, t) for t in vs)))), tuple((r for r in lines if any((dot(r, t) for t in vs))))))
        stack = [(0, I, (), ())]
        while stack:
            (h, A, R, W) = stack.pop()
            j = h % len(ph)
            (r, w) = rules[j]
            R += tuple((tran(v, A) for v in r))
            W += tuple((tran(v, A) for v in w))
            if not mixed(R, W):
                continue
            assert h < depth
            for (H, L) in ph[j]:
                RR = R + tuple((tran(v, A) for v in L))
                if mixed(RR, W):
                    stack.append((h + 1, matrix(H, A), RR, W))
    return (depth, len(profiles))

def normal(row, d):
    """Test every possible uniquely maximizing signed covector and require all the next indicated return norms to be smaller. Every competitor word must be infeasible."""
    (ph, T, z) = family(row)
    L = [tuple(v) for v in d['norm']]
    assert all((not dot(r, t) for r in L for t in T[0]))
    assert sp.Matrix(L).rank() == 5 - len(T[0])
    signed = set(L + [neg(v) for v in L])
    k = len(d['counts']) - 1
    for l in signed:
        stack = [(0, I, tuple((sub(l, s) for s in signed if s != l)))]
        while stack:
            (h, A, C) = stack.pop()
            if h and h % len(ph) == 0:
                C += tuple((sub(l, tran(r, A)) for r in signed))
            if not strict(C):
                continue
            assert h < k * len(ph)
            for (H, L) in ph[h % len(ph)]:
                stack.append((h + 1, matrix(H, A), C + tuple((tran(v, A) for v in L))))
    return k

def unguarded(records, all_guarded, base):
    """Independently close 54+17 unrestricted covector orbits, check every normal, and reconcile the full unpaired catalogue against all three evidence families."""
    import re, ast
    old = Path('docs/drafts/research-direction/astra-strategy-13/basins.out').read_text()
    seeds = {int(k): (ast.literal_eval(v), int(h)) for (k, v, h) in re.findall('^(\\d+) dimension \\d+ period \\d+ row (\\([^)]*\\)) orbit size up to sign (\\d+)', old, re.M)}
    extra = ast.literal_eval(next((x.split(': ', 1)[1] for x in (base / 'e1b_guarded/verification.log').read_text().splitlines() if x.startswith('Additional closed orbits:'))))
    assert len(seeds) == 54 and len(extra) == 17 and (not set(seeds) & set(extra))
    seeds.update({k: (u, h) for (k, (u, h, _)) in extra.items()})
    for (k, (u, h)) in seeds.items():
        row = records[k]
        (phases, T, _) = family(row)
        period = len(phases)
        phase_count = len(row['faces'])
        orbit = {(0, tuple(u))}
        todo = list(orbit)
        while todo:
            (j, v) = todo.pop()
            prev = (j - 1) % phase_count
            assert all((not dot(v, t) for t in T[j % period]))
            for (A, _) in phases[prev % period]:
                item = (prev, can(tran(v, A)))
                if item not in orbit:
                    orbit.add(item)
                    todo.append(item)
            assert len(orbit) <= h
        assert len(orbit) == h
    earlier = {147, 176, 177, 179}
    all_ids = set(seeds) | earlier | {int(k) for k in all_guarded}
    assert len(all_ids) == 71 + 4 + 86
    assert all_ids == {i for (i, r) in enumerate(records) if not r['paired']}
    assert len(records) == 191
    print('Unrestricted covector orbits', len(seeds), 'four earlier exclusions', sorted(earlier))
if __name__ == '__main__':
    base = Path('kaprekar-investigation/eleven_access')
    cats = json.loads((base / 'entry_pi/catalogue.json').read_text())
    dat = json.load(gzip.open(base / 'e1b_guarded/certificate.json.gz', 'rt'))
    n = 0
    m = 0
    for (k, d) in dat.items():
        if d.get('method') != 'covector':
            continue
        check(cats[int(k)], d)
        n += 1
        m += len(d['cover'])
    assert (n, m) == (74, 8130)
    print('Guard-cone certificates and states', n, m)
    for (k, d) in dat.items():
        m = d.get('method')
        r = cats[int(k)]
        if m == 'majorization':
            print('empty by step / profiles', k, major(r, d))
        if m == 'noncontraction':
            print('noncontraction within returns', k, normal(r, d))
    unguarded(cats, dat, base)
