"""Theorem 1 tests (independent code).
(a) real: g_1(f x) < g_1(x) (majorization), exact rationals, k=1..7;
(b) circle inequality (eq. circle) for general C<=L, exact rationals;
(c) integer: S_j(g_B(F_B q)) <= S_j(g_B(q)) + 2 for q_k>0, exhaustive small B
    (n=3,5,7,9) and random large B (n up to 15); n=9 also against +1;
(d) eventual periods of f on random rationals: <= (k+1)!, and for k=4 with
    positive gaps period in {1,3} with limit in the core or in E.
"""
import random, math, itertools, sys
from fractions import Fraction as Fr
from km import *

random.seed(20261117)
ok = True


def S_list(g):
    s = sorted(g, reverse=True)
    out, a = [], 0
    for t in s:
        a += t
        out.append(a)
    return out


# (a)
viol = 0
tested = 0
for k in range(1, 8):
    for _ in range(3000):
        D = random.choice([7, 60, 997, 10 ** 6])
        x = sorted([Fr(random.randint(0, D), D) for _ in range(k)], reverse=True)
        if random.random() < 0.2:  # force ties / boundary values
            x[random.randrange(k)] = x[0]
            x = sorted(x, reverse=True)
        tested += 1
        if not majorized(g1(f_real(tuple(x))), g1(tuple(x))):
            viol += 1
print(f"(a) real majorization g1(fx) < g1(x): {tested} exact rational inputs, k=1..7, violations={viol}")
ok &= viol == 0


# (b) circle inequality with C <= L
def DC(x, C):
    vals = [C] + list(x) + [C - t for t in x]
    d = sorted(vals, reverse=True)
    n = len(d)
    return [d[i] - d[n - 1 - i] for i in range(len(x))]


def gL(x, L):
    return [L - x[0]] + [x[i] - x[i + 1] for i in range(len(x) - 1)] + [x[-1]]


viol = 0
tested = 0
for k in range(1, 7):
    for _ in range(3000):
        D = random.choice([12, 1000, 10 ** 5])
        L = Fr(random.randint(1, D), D)
        C = L * Fr(random.randint(0, D), D)
        x = sorted([C * Fr(random.randint(0, D), D) for _ in range(k)], reverse=True)
        tested += 1
        if not majorized(gL(DC(x, C), L), gL(x, L)):
            viol += 1
print(f"(b) circle inequality g_L(D_C x) < g_L(x), C<=L: {tested} inputs, k=1..6, violations={viol}")
ok &= viol == 0

# (c) integer +2 (and +1 for n=9)
for n, Bmax in ((3, 60), (5, 40), (7, 22), (9, 16)):
    k = n // 2
    mx = -99
    cnt = 0
    for B in range(2, Bmax + 1):
        for q in itertools.combinations_with_replacement(range(B - 1, 0, -1), k):
            q = tuple(sorted(q, reverse=True))
            g0 = S_list(gaps(q, B))
            g1_ = S_list(gaps(F(q, B, n), B))
            mx = max(mx, max(a - b for a, b in zip(g1_, g0)))
            cnt += 1
    print(f"(c) n={n}: exhaustive B=2..{Bmax}, {cnt} positive sources, max_j [S_j(g(Fq))-S_j(g(q))] = {mx}")
    ok &= mx <= (1 if n == 9 else 2)
for n in (3, 5, 7, 9, 11, 13, 15):
    k = n // 2
    mx = -99
    for _ in range(20000):
        B = random.choice([random.randint(2, 50), random.randint(50, 5000), random.randint(10 ** 5, 10 ** 9)])
        q = tuple(sorted([random.randint(1, B - 1) if B > 1 else 1 for _ in range(k)], reverse=True))
        if random.random() < 0.3:  # near-ties
            j = random.randrange(k)
            q = tuple(sorted([min(B - 1, max(1, t + random.randint(-1, 1))) if i == j else t for i, t in enumerate(q)], reverse=True))
        g0 = S_list(gaps(q, B))
        g1_ = S_list(gaps(F(q, B, n), B))
        mx = max(mx, max(a - b for a, b in zip(g1_, g0)))
    print(f"(c) n={n}: 20000 random sources (B up to 1e9), max excess = {mx}")
    ok &= mx <= (1 if n == 9 else 2)

# (d) eventual periods of the real map on rationals
for k in range(1, 7):
    pers = {}
    for _ in range(300 if k < 6 else 100):
        D = random.choice([97, 1000, 65537])
        x = tuple(sorted([Fr(random.randint(1, D - 1), D) for _ in range(k)], reverse=True))
        seen = {}
        t = 0
        while x not in seen:
            seen[x] = t
            x = f_real(x)
            t += 1
        p = t - seen[x]
        pers[p] = pers.get(p, 0) + 1
    bad = [p for p in pers if p > math.factorial(k + 1)]
    print(f"(d) k={k}: eventual periods {dict(sorted(pers.items()))}  (bound (k+1)!={math.factorial(k+1)}) {'OK' if not bad else 'FAIL'}")
    ok &= not bad

# k=4 specifics: positive gaps -> period 1 or 3, limit in core or E
from km import Hmap, in_core, Cmap
import z3  # noqa (not used; keep model pure python)


def in_E_py(r):
    a, b, c, d = r
    if not (Fr(5, 7) <= a <= Fr(3, 4)):
        return False
    p1 = b == 2 * a - 1 and c == 2 * a - 1 and 0 < d <= 1 - a
    p2 = c == 3 * a - 2 and d == 3 - 4 * a and 1 - a <= b < 2 - 2 * a
    p3 = b == a and c + d == 2 - 2 * a and 0 <= d < 1 - a
    return p1 or p2 or p3


stats = {}
entry = []
for _ in range(4000):
    D = random.choice([1009, 10007, 10 ** 6 + 3])
    xs = sorted(random.sample(range(1, D), 4), reverse=True)
    x = tuple(Fr(t, D) for t in xs)
    r = tuple(sorted([abs(2 * t - 1) for t in x], reverse=True))
    seen = {}
    t = 0
    first_core = None
    while r not in seen:
        seen[r] = t
        if first_core is None and in_core(r, strict=True):
            first_core = t
        r = Hmap(r)
        t += 1
    p = t - seen[r]
    where = "core-int" if in_core(r, strict=True) else ("core-bdry" if in_core(r) else ("E" if in_E_py(r) else "OTHER"))
    stats[(p, where)] = stats.get((p, where), 0) + 1
    if first_core is not None:
        entry.append(first_core)
print("(d) k=4 positive gaps, 4000 random rationals: (period, limit location) ->", dict(sorted(stats.items())))
print(f"    entry time e into int(core): max={max(entry)}, mean={sum(entry)/len(entry):.2f}, histogram up to 10:",
      [entry.count(i) for i in range(11)])
ok &= all(p in (1, 3) and w != "OTHER" for (p, w) in stats)
print("THEOREM 1/2-REAL TESTS OK" if ok else "SOME TEST FAILED")
