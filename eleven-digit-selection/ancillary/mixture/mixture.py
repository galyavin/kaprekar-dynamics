"""Independent exact analysis of the six-residue lane mixture.

Inputs are only the definitions and the published tables of proof.md
(transcribed below); the route verification code is separate.

  mixture.py tables              exact table checks, distinctness, ramp, m(U)
  mixture.py samples FILE...     lane attribution of kap.c outputs
  mixture.py box FILE...         threshold-base offset check (all lattice points)
"""
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

NUM = [26293, 27557, 28481, 40447, 29865, 27357]
C = [F(n, 180000) for n in NUM]
EPS = F(1, 10**6)
OFFSET_BOUND = 144

# (v_a, v_b): {beta: [(w or None, mass), ...]}   transcribed from proof.md
TABLE = {
    (F(16, 5), F(2)): {0: [(F(24, 5), F(2651, 5184))],
                       2: [(F(23, 5), F(2651, 6480)), (F(24, 5), F(2651, 25920))],
                       4: [(F(24, 5), F(2651, 12960)), (F(47, 10), F(2651, 8640))]},
    (F(33, 10), F(21, 10)): {0: [(F(22, 5), F(535, 3024)), (F(24, 5), F(535, 4032))],
                             2: [(F(24, 5), F(535, 1728))],
                             4: [(F(24, 5), F(535, 1728))]},
    (F(18, 5), F(12, 5)): {0: [(F(22, 5), F(29, 162))],
                           2: [(F(24, 5), F(29, 162))],
                           4: [(F(24, 5), F(29, 162))],
                           1: [(F(24, 5), F(55, 1296))],
                           3: [(F(22, 5), F(55, 1296))],
                           5: [(F(24, 5), F(55, 1296))]},
    (F(3), F(9, 5)): {1: [(F(24, 5), F(2461, 10368)), (F(114, 25), F(2683, 31104)), (F(117, 25), F(67, 15552))],
                      3: [(F(22, 5), F(737, 5184)), (F(24, 5), F(107, 576))],
                      5: [(F(24, 5), F(9395, 31104)), (F(51, 10), F(805, 31104))]},
    (F(18, 5), F(9, 5)): {1: [(None, F(79, 864))], 3: [(None, F(79, 864))], 5: [(None, F(79, 864))]},
    (F(21, 5), F(12, 5)): {1: [(F(24, 5), F(1571, 5184))],
                           3: [(F(22, 5), F(1571, 13824)), (F(24, 5), F(7855, 41472))],
                           5: [(F(24, 5), F(7855, 41472)), (F(51, 10), F(1571, 13824))]},
    (F(21, 5), F(27, 10)): {1: [(F(24, 5), F(1219, 11664)), (F(114, 25), F(6095, 46656))],
                            3: [(F(24, 5), F(1219, 5184))],
                            5: [(F(24, 5), F(1219, 5184))]},
}
SHEET = {0: {(F(16, 5), F(2)): F(2651, 5184), (F(33, 10), F(21, 10)): F(535, 1728), (F(18, 5), F(12, 5)): F(29, 162)},
         1: {(F(3), F(9, 5)): F(425, 1296), (F(18, 5), F(9, 5)): F(79, 864), (F(18, 5), F(12, 5)): F(55, 1296),
             (F(21, 5), F(12, 5)): F(1571, 5184), (F(21, 5), F(27, 10)): F(1219, 5184)}}


def kappa(va, vb, w):
    if va == 2 * vb:
        assert w is None
        return 1 / (2 * vb)
    return 1 / (2 * vb) + (1 - va / (2 * vb)) / w


def lanes(beta):
    out = []
    for (va, vb), d in TABLE.items():
        for w, p in d.get(beta, []):
            out.append(((va, vb, w), kappa(va, vb, w), p))
    return out


# --- affine forms in the free gaps y_0..y_4 -------------------------------
# x_i = i - 6 sum_{j<i} y_j  (1<=i<=5): represent an affine form as (const, [c0..c4])
def x_form(i):
    return (F(i), [F(-6) if j < i else F(0) for j in range(5)])


def lin(*terms):
    const = F(0)
    co = [F(0)] * 5
    for a, (c0, cs) in terms:
        const += a * c0
        co = [u + a * v for u, v in zip(co, cs)]
    return (const, co)


A_FORM = lin((F(1, 6), x_form(1)), (F(55, 6), x_form(2)), (F(-30, 6), x_form(3)),
             (F(-55, 6), x_form(4)), (F(95, 6), x_form(5)))
P_FORM = lin((F(1, 6), x_form(1)), (F(-1, 6), A_FORM))


def clock_form(k):
    return lin((F(1), P_FORM), (k, A_FORM))


def ev(form, y):
    c0, cs = form
    return c0 + sum(a * b for a, b in zip(cs, y))


def range_on_box(form):
    """exact min, max of an affine form over the closed box |y_i-c_i|<=EPS"""
    mid = ev(form, C[:5])
    rad = EPS * sum(abs(a) for a in form[1])
    return mid - rad, mid + rad


def m_box(a):
    p = F(1)
    for i in range(5):
        p *= C[i]
    return math.factorial(11) * (2 * a) ** 5 * p * (C[5] - a * a / 3 * sum(1 / C[i] for i in range(5)))


def simpson_mass(a):
    """tensor Simpson rule, exact for degree<=3 per variable (independent of m_box)"""
    nodes = [(-a, F(1, 3)), (F(0), F(4, 3)), (a, F(1, 3))]   # weights * (a) per dim
    tot = F(0)

    def rec(i, y, wt):
        nonlocal tot
        if i == 5:
            y5 = 1 - sum(y)
            p = y5
            for v in y:
                p *= v
            tot += wt * p
            return
        for d, w in nodes:
            rec(i + 1, y + [C[i] + d], wt * w * a)
    rec(0, [], F(1))
    return math.factorial(11) * tot


def chi2_sf(x, k):
    """upper tail of chi-square with k dof (regularized gamma Q(k/2, x/2))"""
    s, z = k / 2.0, x / 2.0
    if z <= 0:
        return 1.0
    if z < s + 1:
        term = 1.0 / s
        tot = term
        n = 1
        while abs(term) > 1e-15 * abs(tot):
            term *= z / (s + n)
            tot += term
            n += 1
        return max(0.0, 1.0 - tot * math.exp(-z + s * math.log(z) - math.lgamma(s)))
    # continued fraction (Lentz)
    b = z + 1 - s
    c = 1e300
    d = 1 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - s)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-15:
            break
    return math.exp(-z + s * math.log(z) - math.lgamma(s)) * h


def fmt(q):
    return f"{q} ({float(q):.9g})"


def tables():
    print("== table sums and sheet marginals")
    for beta in range(6):
        L = lanes(beta)
        tot = sum(p for _, _, p in L)
        sheet = defaultdict(F)
        for (va, vb, w), k, p in L:
            sheet[(va, vb)] += p
        ok = dict(sheet) == SHEET[beta % 2]
        print(f"beta {beta}: {len(L)} lanes, total {tot}, sheet marginal matches stated parity masses: {ok}")
        assert tot == 1 and ok
    print("\n== clock coefficients kappa and clock forms t = P + kappa A on closed U")
    amin, amax = range_on_box(A_FORM)
    pmin, pmax = range_on_box(P_FORM)
    print(f"A on closed U: [{float(amin):.9f}, {float(amax):.9f}]  (A>0: {amin > 0})")
    print(f"P on closed U: [{float(pmin):.9f}, {float(pmax):.9f}]")
    print(f"A coefficients in y: {[str(v) for v in A_FORM[1]]}, const {A_FORM[0]}")
    kap = {}
    for beta in range(6):
        kap[beta] = Counter()
        for lane, k, p in lanes(beta):
            kap[beta][k] += p
    for beta in range(6):
        for lane, k, p in lanes(beta):
            tmin, tmax = range_on_box(clock_form(k))
            others = [b for b in range(6) if b != beta and b % 2 == beta % 2 and k in kap[b]]
            print(f"beta {beta} lane {tuple(str(v) for v in lane)} kappa {k} = {float(k):.6f} "
                  f"mass {p} t in [{float(tmin):.7f},{float(tmax):.7f}] same-parity residues sharing kappa: {others}")
    print("\n== pairwise distinctness of the six limits")
    # separating test: kappa present (positive mass) in b1, absent in b2.
    for b1 in range(6):
        for b2 in range(b1 + 1, 6):
            found = None
            for (x, y) in ((b1, b2), (b2, b1)):
                for k, p in sorted(kap[x].items(), key=lambda kv: -kv[1]):
                    if k not in kap[y]:
                        sep = min(abs(k - k2) for k2 in kap[y]) * amin   # min graph separation on closed U
                        found = (x, y, k, p, sep)
                        break
                if found:
                    break
            x, y, k, p, sep = found
            print(f"({b1},{b2}): h(y,t)=max(0,1-|t-(P+{k}A)|/delta), delta<{float(sep):.3e}: "
                  f"E_{x}[h] >= {p} = {float(p):.4f}, E_{y}[h] = 0")
    print("\n== fixed descending ramp psi(t)=clamp((131/5000-t)/(131/5000-13/500))")
    lo_r, hi_r = F(13, 500), F(131, 5000)
    expect = {}
    for beta in range(6):
        e = F(0)
        for lane, k, p in lanes(beta):
            tmin, tmax = range_on_box(clock_form(k))
            if tmax <= lo_r:
                e += p
            elif tmin >= hi_r:
                pass
            else:
                raise SystemExit(f"lane {lane} straddles the ramp")
        expect[beta] = e
        print(f"beta {beta}: E[psi] = {e}")
    fast_max = max(range_on_box(clock_form(k))[1] for b in range(6) for _, k, _ in lanes(b)
                   if range_on_box(clock_form(k))[1] <= lo_r)
    slow_min = min(range_on_box(clock_form(k))[0] for b in range(6) for _, k, _ in lanes(b)
                   if range_on_box(clock_form(k))[0] >= hi_r)
    print(f"largest fast-lane time {float(fast_max):.7f}; smallest other-lane time {float(slow_min):.7f}")
    contrast = expect[1] - expect[0]
    assert all(expect[b] == expect[b % 2] for b in range(6))
    print(f"odd-minus-even conditional contrast = {fmt(contrast)}")
    print("\n== m(U)")
    m = m_box(EPS)
    s = simpson_mass(EPS)
    print(f"closed form m(eps) = {float(m):.12e}; tensor Simpson (exact) = {float(s):.12e}; equal: {m == s}")
    # bump phi = clamp(3-4||y-c||/eps): layer cake I_phi = int_0^1 m((3-s)eps/4) ds; m is a polynomial in a
    # integrate exactly: m(a) = K a^5 (c5 - a^2 S/3) with K = 11! 32 prod c
    p = F(1)
    for i in range(5):
        p *= C[i]
    Kc = math.factorial(11) * 32 * p
    S = sum(1 / C[i] for i in range(5))
    # int_0^1 ((3-s)eps/4)^n ds = (eps/4)^n (3^{n+1}-2^{n+1})/(n+1)
    def mom(n):
        return (EPS / 4) ** n * F(3 ** (n + 1) - 2 ** (n + 1), n + 1)
    I_phi = Kc * (C[5] * mom(5) - S / 3 * mom(7))
    print(f"I_phi = {float(I_phi):.12e}; (155/288) I_phi = {float(contrast * I_phi):.12e}")
    # brute check of I_phi by midpoint quadrature in s
    n = 2000
    I_num = sum(m_box((3 - (j + F(1, 2)) / n) * EPS / 4) for j in range(n)) / n
    print(f"I_phi midpoint check (n={n}): {float(I_num):.12e}")


def parse(files):
    rows = []
    for fn in files:
        with open(fn) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                v = list(map(int, line.split()))
                rows.append(v)
    return rows


def in_open_box(B, G):
    return all(abs(F(G[i], B) - C[i]) < EPS for i in range(5)) and sum(G) == B


def attribute(B, G, T):
    beta = B % 6
    X = [i * B - 6 * sum(G[:i]) for i in range(1, 6)]
    AB = F(X[0] + 55 * X[1] - 30 * X[2] - 55 * X[3] + 95 * X[4], 6)
    PB = (X[0] - AB) / 6
    res = []
    for lane, k, p in lanes(beta):
        off = T - (PB + k * AB)
        res.append((abs(off), off, lane, k))
    res.sort(key=lambda r: r[0])
    return res


def samples(files):
    rows = parse(files)
    by = defaultdict(list)
    for r in rows:
        by[(r[0] % 6, r[0])].append(r)
    agg = defaultdict(Counter)
    aggn = Counter()
    bad = 0
    offsets = defaultdict(set)
    maxoff = 0
    lam_seen = Counter()
    rel_err = 0.0
    for (beta, B), rs in sorted(by.items()):
        cnt = Counter()
        for r in rs:
            G, T, lam = r[1:7], r[7], r[8]
            assert in_open_box(B, G), (B, G)
            lam_seen[lam] += 1
            res = attribute(B, G, T)
            hits = [x for x in res if x[0] < OFFSET_BOUND]
            if len(hits) != 1:
                bad += 1
                print(f"UNATTRIBUTED B={B} G={G} T={T} nearest {[(float(a), b) for a, _, b, _ in res[:2]]}")
                continue
            a, off, lane, k = hits[0]
            gap2 = res[1][0] if len(res) > 1 else None
            cnt[lane] += 1
            offsets[(beta, lane)].add(off)
            maxoff = max(maxoff, a)
            y = [F(g, B) for g in G]
            rel_err = max(rel_err, abs(float(F(T, B) - ev(clock_form(k), y))))
        n = sum(cnt.values())
        agg_line = []
        chi = 0.0
        L = lanes(beta)
        for lane, k, p in L:
            o = cnt[lane]
            e = n * p
            chi += (o - e) ** 2 / e
            agg[beta][lane] += o
            agg_line.append(f"{o}/{float(e):.1f}")
        aggn[beta] += n
        extra = set(cnt) - {l for l, _, _ in L}
        print(f"beta {beta} B {B}: n={n} observed/expected {' '.join(agg_line)} "
              f"chi2={chi:.2f} df={len(L)-1} p={chi2_sf(chi, len(L)-1):.3f}{' EXTRA '+str(extra) if extra else ''}")
    print("\n== pooled over bases")
    for beta in range(6):
        if not aggn[beta]:
            continue
        n = aggn[beta]
        chi = 0.0
        L = lanes(beta)
        parts = []
        for lane, k, p in L:
            o = agg[beta][lane]
            e = n * p
            chi += (o - e) ** 2 / e
            z = (o - e) / math.sqrt(n * float(p) * (1 - float(p)))
            parts.append(f"{tuple(str(v) for v in lane)}: {o/n:.4f} vs {float(p):.4f} (z={z:+.2f})")
        print(f"beta {beta}: n={n} chi2={chi:.2f} df={len(L)-1} p={chi2_sf(chi, len(L)-1):.3f}")
        for s in parts:
            print("   " + s)
    print(f"\nunattributed samples: {bad}; max |T - B t_lane(G/B)| = {float(maxoff):.4f}; "
          f"max |T/B - t_lane| = {rel_err:.3e}; cycle lengths seen {dict(lam_seen)}")
    print("lane independence of position in U (pooled per beta; halves split on y_i - c_i, i=0..4):")
    for beta in range(6):
        rs = [r for r in rows if r[0] % 6 == beta]
        if not rs:
            continue
        L = [l for l, _, _ in lanes(beta)]
        worst = 1.0
        for i in range(5):
            tab = {h: Counter() for h in (0, 1)}
            for r in rs:
                B, G, T = r[0], r[1:7], r[7]
                hits = [x for x in attribute(B, G, T) if x[0] < OFFSET_BOUND]
                if len(hits) != 1:
                    continue
                tab[int(F(G[i], B) > C[i])][hits[0][2]] += 1
            n0, n1 = sum(tab[0].values()), sum(tab[1].values())
            chi = 0.0
            df = 0
            for l in L:
                tot = tab[0][l] + tab[1][l]
                if tot == 0:
                    continue
                df += 1
                for h, nh in ((0, n0), (1, n1)):
                    e = tot * nh / (n0 + n1)
                    chi += (tab[h][l] - e) ** 2 / e
            worst = min(worst, chi2_sf(chi, df - 1))
        print(f"   beta {beta}: smallest p over the five half-splits = {worst:.3f}")
    print("distinct exact offsets per lane (count, min, max):")
    for key in sorted(offsets, key=lambda k: (k[0], str(k[1]))):
        s = offsets[key]
        print(f"   beta {key[0]} {tuple(str(v) for v in key[1])}: {len(s)} values in "
              f"[{float(min(s)):.4f}, {float(max(s)):.4f}]")
    return offsets


def box(files):
    rows = parse(files)
    nb = 0
    worst = 0
    amb = 0
    perbeta = Counter()
    for r in rows:
        B, G, T = r[0], r[1:7], r[7]
        assert in_open_box(B, G)
        res = attribute(B, G, T)
        hits = [x for x in res if x[0] < OFFSET_BOUND]
        perbeta[B % 6] += 1
        if not hits:
            nb += 1
            print(f"NO LANE within {OFFSET_BOUND}: B={B} G={G} T={T} nearest {float(res[0][0]):.3f}")
        else:
            worst = max(worst, res[0][0])
            if len(hits) > 1:
                amb += 1
    print(f"lattice points: {len(rows)} by residue {dict(sorted(perbeta.items()))}; "
          f"without any lane within {OFFSET_BOUND}: {nb}; ambiguous (>1 lane within bound): {amb}; "
          f"max nearest |offset| = {float(worst):.4f}")


if __name__ == "__main__":
    if sys.argv[1] == "tables":
        tables()
    elif sys.argv[1] == "samples":
        samples(sys.argv[2:])
    elif sys.argv[1] == "box":
        box(sys.argv[2:])
