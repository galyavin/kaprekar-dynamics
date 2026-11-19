#!/usr/bin/env python3
"""Independent odd-cohort checks.

Written from the gap and route definitions.

Step 1 (black-box lattice).  At an even base B = 3360 mod 3840 near 4.2e9,
take P = G* (mod 3840) inside the closed rectangle (the witness
residues mod 3840 lie in the cohort).  Run the centered map (./cen) at P and
at P + 3840 e_i (i < 5), on the even arm (B, G) and the odd arm (B+1, G+e4).
Each orbit is compressed into translation words; the repeat counts of the
words are read off.  On the cohort these counts are affine in G, so their
slopes are (difference)/3840.  The cohort lattice is
    L = {v in Z^5 : every repeat-count slope . v is an integer}
(integrality of the affine repeat counts), and its
index is computed by imposing one row at a time.  Extra points test that the
counts are affine (predicted vs observed).

Step 2 (clocks).  Every run's T is compared with mu+1 for proof.md's exact
cohort clocks mu_e, mu_o (in even centered X).  Then cohort points P + L in
the HALF-size closed rectangle U' at several bases are simulated with the
independent whole-word program ./ww on both arms; we check T_e = mu_e+1,
T_o = mu_o+1, T_o/(B+1) <= 13/500, T_e/B > 131/5000.

Step 3 (mass).  m' from the index, by exact Simpson cubature (exact for this
integrand) and by the closed form, and a lattice count of P + L in U'.

Usage: cohort_check.py CEN WW STEP
"""
import random
import subprocess
import sys
from fractions import Fraction as F
from math import ceil, floor, gcd, factorial

GS = (26293, 27557, 28481, 40447, 29865, 27357)
C = [F(g, 180000) for g in GS]
EPS = F(1, 10**6)
MU_E = ([F(3, 16), F(55, 48), F(-5, 8), F(-55, 48), F(95, 48)], F(39))
MU_O = ([F(205, 1152), F(715, 1152), F(-65, 192), F(-715, 1152), F(1235, 1152)], F(209, 24))


def Xof(B, G):
    s, X = 0, []
    for i in range(5):
        s += G[i]
        X.append((i + 1) * B - 6 * s)
    return X


def clock(form, X):
    a, k = form
    return sum(coefficient * xi for coefficient, xi in zip(a, X)) + k


def box(B, eps):
    return ([ceil((C[i] - eps) * B) for i in range(5)], [floor((C[i] + eps) * B) for i in range(5)])


def run_cen(cen, jobs):
    inp = "".join(f"{B} {' '.join(map(str, G))}\n" for B, G in jobs)
    out = subprocess.run([cen], input=inp, capture_output=True, text=True, check=True).stdout
    res, cur = [], None
    for line in out.splitlines():
        f = line.split()
        if f[0] == "T":
            cur = dict(T=int(f[1]), MU=int(f[3]), LAM=int(f[5]), toks=[])
        elif f[0] == "END":
            res.append(cur)
        elif f[0] == "FAIL":
            raise SystemExit(line)
        else:
            cur["toks"].append((int(f[0]), int(f[1]), f[2], f[3]))
    assert len(res) == len(jobs)
    return res


def run_ww(ww, jobs, workers=5):
    from concurrent.futures import ThreadPoolExecutor
    chunks = [jobs[k::workers] for k in range(workers)]

    def go(ch):
        if not ch:
            return []
        inp = "".join(f"{B} {' '.join(map(str, G))}\n" for B, G in ch)
        out = subprocess.run([ww], input=inp, capture_output=True, text=True, check=True).stdout
        r = []
        for line in out.splitlines():
            v = list(map(int, line.split()))
            r.append(((v[0], tuple(v[1:6])), v[6], v[7], v[8]))
        return r

    with ThreadPoolExecutor(workers) as ex:
        res = {}
        for part in ex.map(go, chunks):
            for key, T, MUD, LAM in part:
                res[key] = (T, MUD, LAM)
    return res


def long_words(toks):
    """(p, delta, word, count) for translation words; skeleton = all tokens
    with counts of long words blanked."""
    lw = [(p, d, w, c) for p, c, d, w in toks if d != "(0,0,0,0,0)"]
    skel = tuple((p, d, w, (c if d == "(0,0,0,0,0)" else None)) for p, c, d, w in toks)
    return lw, skel


def impose(basis, row, mod):
    """basis: list of 5 integer column vectors; impose row.v = 0 (mod mod)."""
    vals = [sum(r * b for r, b in zip(row, v)) % mod for v in basis]
    basis = [list(v) for v in basis]
    # unimodular column ops: gcd-reduce vals into position 0
    for j in range(1, len(basis)):
        while vals[j] != 0:
            q = vals[0] // vals[j] if vals[j] else 0
            basis[0] = [a - q * b for a, b in zip(basis[0], basis[j])]
            vals[0] = vals[0] - q * vals[j]
            basis[0], basis[j] = basis[j], basis[0]
            vals[0], vals[j] = vals[j], vals[0]
    g = gcd(vals[0], mod)
    f = mod // g
    basis[0] = [f * a for a in basis[0]]
    return basis, f


def hnf_lower(basis):
    """column-style lower-triangular HNF of a full-rank 5x5 lattice basis."""
    cols = [list(v) for v in basis]
    n = 5
    H = []
    for i in range(n):
        # among remaining cols, gcd-reduce coordinate i
        rem = cols
        while True:
            nz = [c for c in rem if c[i] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda c: abs(c[i]))
            piv = nz[0]
            new = [piv]
            for c in rem:
                if c is piv:
                    continue
                if c[i] != 0:
                    q = c[i] // piv[i]
                    c = [a - q * b for a, b in zip(c, piv)]
                new.append(c)
            rem = new
        piv = [c for c in rem if c[i] != 0][0]
        if piv[i] < 0:
            piv = [-a for a in piv]
        H.append(piv)
        cols = [c for c in rem if c[i] == 0]
    # reduce below-diagonal entries
    for i in range(n):
        for j in range(i):
            q = H[j][i] // H[i][i]
            H[j] = [a - q * b for a, b in zip(H[j], H[i])]
    return H


def coset_points(P, H, lo, hi):
    """all points P + sum a_j H[j] in the box (H lower-triangular columns)."""
    out = []

    def rec(i, v):
        if i == 5:
            out.append(tuple(v))
            return
        # coordinate i: v_i + a * H[i][i] within [lo_i, hi_i]
        h = H[i][i]
        amin = ceil(F(lo[i] - v[i], h))
        amax = floor(F(hi[i] - v[i], h))
        for a in range(amin, amax + 1):
            rec(i + 1, [x + a * y for x, y in zip(v, H[i])])

    # shift P to be reduced first: start from P itself
    rec(0, list(P))
    return out


def step1(cen, ww):
    B = 3840 * 1093750 + 3360
    assert B % 3840 == 3360
    lo, hi = box(B, EPS)
    P = [lo[i] + ((GS[i] - lo[i]) % 3840) for i in range(5)]
    assert all(lo[i] <= P[i] and P[i] + 3840 <= hi[i] for i in range(5))
    pts = [P] + [[P[j] + (3840 if j == i else 0) for j in range(5)] for i in range(5)]
    checks = [[P[0] + 3840, P[1] + 3840, P[2], P[3], P[4]],
              [P[0], P[1], P[2] + 3840, P[3] + 3840, P[4] + 3840]]
    allpts = pts + checks
    print(f"STEP1 B={B} (mod 3840 = {B % 3840}); P={P}; width per gap = {hi[0]-lo[0]}")
    ev = run_cen(cen, [(B, G) for G in allpts])
    od = run_cen(cen, [(B + 1, G[:4] + [G[4] + 1]) for G in allpts])
    rows = []
    for arm, res, form in (("even", ev, MU_E), ("odd", od, MU_O)):
        lws = [long_words(r["toks"]) for r in res]
        skel0 = lws[0][1]
        same = all(s == skel0 for _, s in lws)
        print(f" {arm}: identical token skeleton at all 8 points: {same}; tokens={len(skel0)}; "
              f"long words={len(lws[0][0])}")
        for G, r in zip(allpts, res):
            X = Xof(B, G)
            pred = clock(form, X) + 1
            ok = (pred == r["T"])
            print(f"   G={G} T={r['T']} mu+1(formula)={pred} match={ok} period={r['LAM']} "
                  f"T/B={r['T'] / (B + (arm == 'odd')):.9f}")
        n0 = [c for _, _, _, c in lws[0][0]]
        print(f"   long words (p, delta, count at P):")
        for (p, d, w, c) in lws[0][0]:
            print(f"     p={p} delta={d} count={c}")
        slopes = []
        for k in range(len(n0)):
            s = [F(lws[1 + i][0][k][3] - n0[k], 3840) for i in range(5)]
            slopes.append(s)
            print(f"   word {k}: count slope per unit G_i = {[str(x) for x in s]}")
        # affine prediction at check points
        for ci, G in enumerate(checks):
            v = [G[j] - P[j] for j in range(5)]
            obs = [c for _, _, _, c in lws[6 + ci][0]]
            pr = [n0[k] + sum(slopes[k][j] * v[j] for j in range(5)) for k in range(len(n0))]
            print(f"   affine check point {ci}: predicted {[str(x) for x in pr]} observed {obs} "
                  f"match={all(F(o) == p for o, p in zip(obs, pr))}")
        rows += slopes
    # lattice: impose each slope row (scaled to integers) on Z^5
    basis = [[1 if i == j else 0 for j in range(5)] for i in range(5)]
    factors = []
    for s in rows:
        D = 1
        for x in s:
            D = D * x.denominator // gcd(D, x.denominator)
        row = [int(x * D) for x in s]
        basis, f = impose(basis, row, D)
        factors.append(f)
    index = 1
    for f in factors:
        index *= f
    print(f" per-word index factors (even words then odd words, in route order): {factors}")
    print(f" joint cohort lattice index = {index}   (certified index 3317760)")
    H = hnf_lower(basis)
    det = 1
    for i in range(5):
        det *= H[i][i]
    print(f" HNF diagonal {[H[i][i] for i in range(5)]}, |det| = {abs(det)}")
    print(f" HNF basis columns: {H}")
    # 3840 Z^5 inside L?
    inside = all(all(int(sum(r * b for r, b in zip(s_row, [3840 if j == i else 0 for j in range(5)]))) == sum(r * b for r, b in zip(s_row, [3840 if j == i else 0 for j in range(5)])) for s_row in rows) for i in range(5))
    print(f" 3840 Z^5 contained in L: {inside}")
    # whole-word confirmation at P on both arms
    r = run_ww(ww, [(B, tuple(P)), (B + 1, tuple(P[:4] + [P[4] + 1]))], workers=2)
    print(f" whole-word at P: even T={r[(B, tuple(P))][0]} (centered {ev[0]['T']}), "
          f"odd T={r[(B + 1, tuple(P[:4] + [P[4] + 1]))][0]} (centered {od[0]['T']})")
    return P, H, index


def coset_sample(P, H, lo, hi, n, rng):
    """n random points of P + L in the box, choosing a_0, a_1, ... in turn."""
    out = []
    while len(out) < n:
        v = list(P)
        ok = True
        for i in range(5):
            h = H[i][i]
            amin = ceil(F(lo[i] - v[i], h))
            amax = floor(F(hi[i] - v[i], h))
            if amax < amin:
                ok = False
                break
            a = rng.randint(amin, amax)
            v = [x + a * y for x, y in zip(v, H[i])]
        if ok:
            out.append(tuple(v))
    return out


def step2(ww, P, H, specs):
    rng = random.Random(57057)
    for B, nmax in specs:
        assert B % 3840 == 3360
        lo, hi = box(B, EPS / 2)
        # translate P into the box region: P is a coset rep valid for every B in the progression
        expected = (EPS * B) ** 5 / 3317760
        if expected > 10**5:
            pts = coset_sample(P, H, lo, hi, nmax, rng)
            npts = f"~{float(expected):.3g} (sampled)"
        else:
            pts = coset_points(P, H, lo, hi)
            npts = len(pts)
            if nmax > 0 and npts > nmax:
                pts = rng.sample(pts, nmax)
        jobs = []
        for G in pts:
            jobs.append((B, tuple(G)))
            jobs.append((B + 1, tuple(list(G[:4]) + [G[4] + 1])))
        res = run_ww(ww, jobs)
        bad_e = bad_o = ramp_o = ramp_e = 0
        maxo = None
        mine = None
        for G in pts:
            X = Xof(B, list(G))
            Te = res[(B, tuple(G))][0]
            To = res[(B + 1, tuple(list(G[:4]) + [G[4] + 1]))][0]
            bad_e += Te != clock(MU_E, X) + 1
            bad_o += To != clock(MU_O, X) + 1
            ramp_o += not (F(To, B + 1) <= F(13, 500))
            ramp_e += not (F(Te, B) > F(131, 5000))
            to = F(To, B + 1)
            te = F(Te, B)
            maxo = to if maxo is None or to > maxo else maxo
            mine = te if mine is None or te < mine else mine
        print(f"STEP2 B={B}: cohort points in closed U' = {npts}; simulated {len(pts)} (both arms)")
        if pts:
            print(f"   even T != mu_e+1: {bad_e}; odd T != mu_o+1: {bad_o}; "
                  f"odd T/(B+1) > 13/500: {ramp_o}; even T/B <= 131/5000: {ramp_e}")
            print(f"   max odd T/(B+1) = {float(maxo):.9f}; min even T/B = {float(mine):.9f}")
        sys.stdout.flush()


def step3(index, P, H):
    a = EPS / 2
    pc = 1
    for i in range(5):
        pc *= C[i]
    closed = F(factorial(11), index) * (2 * a) ** 5 * pc * (C[5] - a * a / 3 * sum(1 / C[i] for i in range(5)))
    # Simpson cubature in 5 dims on [c_i - a, c_i + a]: exact for degree <= 3 per variable
    nodes = [(-a, a / 3), (F(0), 4 * a / 3), (a, a / 3)]
    tot = F(0)
    import itertools
    for combo in itertools.product(nodes, repeat=5):
        y = [C[i] + combo[i][0] for i in range(5)]
        w = 1
        for t in combo:
            w *= t[1]
        y5 = 1 - sum(y)
        prod = y5
        for v in y:
            prod *= v
        tot += w * prod
    simpson = F(factorial(11), index) * tot
    print(f"STEP3 m' closed form = {float(closed):.12e}; Simpson (exact) = {float(simpson):.12e}; equal: {closed == simpson}")
    print(f"   m' exact = {closed}")
    print(f"   m' > 0: {closed > 0}")
    # lattice count check of P + L in U' at a base: count * weight vs m'
    for B in (3840 * 26041 + 3360, 3840 * 52083 + 3360):
        lo, hi = box(B, a)
        pts = coset_points(P, H, lo, hi)
        wsum = F(0)
        for G in pts:
            G5 = B - sum(G)
            w = factorial(11) * G[0]
            for g in list(G[1:]) + [G5]:
                w *= (g - 1)
            wsum += w
        # include odd translation: weights at B+1 of G+e4 (distinct-digit lower bound)
        wodd = F(0)
        for G in pts:
            Gp = list(G[:4]) + [G[4] + 1]
            G5 = B + 1 - sum(Gp)
            w = factorial(11) * Gp[0]
            for g in Gp[1:] + [G5]:
                w *= (g - 1)
            wodd += w
        pe = wsum / (B ** 11 - B)
        po = wodd / ((B + 1) ** 11 - (B + 1))
        print(f"   B={B}: cohort points in closed U' = {len(pts)} (expected ~ {float((EPS * B) ** 5 / index):.1f}); "
              f"distinct-digit mass even = {float(pe):.6e}, odd partner = {float(po):.6e}, ratio to m' = {float(po / closed):.4f}")


def main():
    cen, ww, which = sys.argv[1], sys.argv[2], sys.argv[3]
    P, H, index = step1(cen, ww)
    if which in ("2", "all"):
        step2(ww, P, H, [(3840 * 7813 + 3360, -1), (3840 * 13021 + 3360, -1),
                         (3840 * 26041 + 3360, 300), (3840 * 260416 + 3360, 40)])
    if which in ("3", "all"):
        step3(index, P, H)


if __name__ == "__main__":
    main()
