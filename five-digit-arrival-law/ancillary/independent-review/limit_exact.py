"""Exact rational evaluation of the limiting law of the paper, written by the
reviewer from the paper's definitions (13)-(16) and (25)-(29) only.

Method (own implementation): convex half-plane clipping with Fractions and
exact integration of polynomials over triangles via
    int_{simplex} s^i t^j ds dt = i! j! / (i+j+2)!.
(The paper's implementation uses Green's theorem on polygon edges instead.)

usage: python limit_exact.py [N_for_mean]
"""
import math
import sys
from fractions import Fraction as Fr

HALF = Fr(1, 2)


# ---------------- polynomials in x,y as {(i,j): coeff} -----------------------
def padd(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + c * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            k = (a + c, b + d)
            r[k] = r.get(k, 0) + v * w
    return {k: v for k, v in r.items() if v != 0}


def aff(cx, cy, c0):
    """affine form cx*x + cy*y + c0 as polynomial"""
    return {k: v for k, v in {(1, 0): Fr(cx), (0, 1): Fr(cy), (0, 0): Fr(c0)}.items() if v != 0}


X, Y, ONE = aff(1, 0, 0), aff(0, 1, 0), aff(0, 0, 1)
RHO = pmul(pmul(pmul({(0, 0): Fr(120)}, padd(ONE, X, -1)), Y), padd(X, Y, -1))

FACT = [math.factorial(i) for i in range(40)]


def integrate_triangle(poly, P0, P1, P2):
    (x0, y0), (x1, y1), (x2, y2) = P0, P1, P2
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if det == 0:
        return Fr(0)
    # x = x0 + s*(x1-x0) + t*(x2-x0), y = y0 + s*(y1-y0) + t*(y2-y0)
    xs = {(0, 0): x0, (1, 0): x1 - x0, (0, 1): x2 - x0}
    ys = {(0, 0): y0, (1, 0): y1 - y0, (0, 1): y2 - y0}
    xs = {k: v for k, v in xs.items() if v != 0}
    ys = {k: v for k, v in ys.items() if v != 0}
    maxi = max(i for i, _ in poly) if poly else 0
    maxj = max(j for _, j in poly) if poly else 0
    xp = [{(0, 0): Fr(1)}]
    for _ in range(maxi):
        xp.append(pmul(xp[-1], xs))
    yp = [{(0, 0): Fr(1)}]
    for _ in range(maxj):
        yp.append(pmul(yp[-1], ys))
    tot = Fr(0)
    for (i, j), c in poly.items():
        st = pmul(xp[i], yp[j])
        for (a, b), w in st.items():
            tot += c * w * Fr(FACT[a] * FACT[b], FACT[a + b + 2])
    return tot * abs(det)


def integrate_polygon(poly, pts):
    tot = Fr(0)
    for k in range(1, len(pts) - 1):
        tot += integrate_triangle(poly, pts[0], pts[k], pts[k + 1])
    return tot


def clip(pts, form):
    """keep form >= 0, form = (cx, cy, c0)"""
    cx, cy, c0 = form
    out = []
    n = len(pts)
    for k in range(n):
        P, Q = pts[k], pts[(k + 1) % n]
        fp = cx * P[0] + cy * P[1] + c0
        fq = cx * Q[0] + cy * Q[1] + c0
        if fp >= 0:
            out.append(P)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    # remove consecutive duplicates
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    if len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res if len(res) >= 3 else []


def neg(f):
    return (-f[0], -f[1], -f[2])


def lin(f, g, a=1, b=1):
    return (a * f[0] + b * g[0], a * f[1] + b * g[1], a * f[2] + b * g[2])


# affine forms in (x,y) are triples (cx, cy, c0)
TRI = [(Fr(0), Fr(1), Fr(0)), (Fr(1), Fr(-1), Fr(0)), (Fr(-1), Fr(0), Fr(1))]  # y>0, x-y>0, 1-x>0
SQUARE = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(1), Fr(1)), (Fr(0), Fr(1))]

# initial cells of (13): guards and (U, V) as affine forms
CELLS = [
    ("x<1/2", [(Fr(-1), Fr(0), HALF)], (Fr(0), Fr(-1), Fr(1)), (Fr(-1), Fr(-1), Fr(1))),
    ("x>=1/2,x+y<1", [(Fr(1), Fr(0), -HALF), (Fr(-1), Fr(-1), Fr(1))],
     (Fr(0), Fr(-1), Fr(1)), (Fr(1), Fr(-1), Fr(0))),
    ("x+y>=1,y<1/2", [(Fr(1), Fr(1), Fr(-1)), (Fr(0), Fr(-1), HALF)],
     (Fr(1), Fr(0), Fr(0)), (Fr(1), Fr(-1), Fr(0))),
    ("y>=1/2", [(Fr(0), Fr(1), -HALF)], (Fr(1), Fr(0), Fr(0)), (Fr(1), Fr(1), Fr(-1))),
]


def base_polygon(guards):
    pts = SQUARE
    for f in TRI + guards:
        pts = clip(pts, f)
        if not pts:
            break
    return pts


def pieces(N):
    """Yield (polygon, A, C) over retained region delta >= eps, eps=1/(2(N+1)),
    for each initial cell and shear count n in 0..N+1 (n=N+1 must be empty)."""
    eps = Fr(1, 2 * (N + 1))
    keep = [(Fr(-1), Fr(0), 1 - eps), (Fr(0), Fr(1), -eps)]   # 1-x>=eps, y>=eps
    for name, guards, Uf, Vf in CELLS:
        P0 = base_polygon(guards + keep)
        if not P0:
            continue
        delta = lin((Fr(0), Fr(0), Fr(1)), Uf, 1, -1)            # 1-U
        for n in range(0, N + 2):
            if n == 0:
                g = [lin((Fr(0), Fr(0), HALF), Vf, 1, -1)]        # 1/2 - V >= 0
                w = Vf
            else:
                g = [lin(lin(Vf, delta, 1, -(n - 1)), (Fr(0), Fr(0), HALF), 1, -1),   # V-(n-1)d-1/2 > 0
                     lin(lin((Fr(0), Fr(0), HALF), Vf, 1, -1), delta, 1, n)]           # 1/2-V+n d >= 0
                w = lin(Vf, delta, 1, -n)
            P = P0
            for f in g:
                P = clip(P, f)
                if not P:
                    break
            if not P:
                continue
            A = lin(lin(Uf, w), (Fr(0), Fr(0), Fr(1)), 1, -1)      # U + w - 1
            C = lin(lin(Uf, w, 2, -1), (Fr(0), Fr(0), Fr(1)), 1, -1)  # 2U - w - 1
            yield name, n, P, A, C


def R(e):
    return 20 * e**2 - 40 * e**3 + 32 * e**5


def cdf(t):
    t = Fr(t)
    assert 0 < t < Fr(1, 3)
    N = max(1, math.ceil(1 / (2 * (Fr(1, 3) - t))) - 1)
    eps = Fr(1, 2 * (N + 1))
    assert eps <= Fr(1, 3) - t and eps <= Fr(1, 4)
    r = Fr(3, 2) * t
    total_mass = Fr(0)
    F = Fr(0)
    nmax_used = 0
    for name, n, P, A, C in pieces(N):
        m = integrate_polygon(RHO, P)
        if m:
            nmax_used = max(nmax_used, n)
        total_mass += m
        Q = P
        AC = lin(A, C, 1, -1)
        for f in [lin((0, 0, r), A, 1, -1), lin((0, 0, r), A, 1, 1),
                  lin((0, 0, r), C, 1, -1), lin((0, 0, r), C, 1, 1),
                  lin((0, 0, r), AC, 1, -1), lin((0, 0, r), AC, 1, 1)]:
            Q = clip(Q, f)
            if not Q:
                break
        if Q:
            F += integrate_polygon(RHO, Q)
    return F, N, total_mass, 1 - R(eps), nmax_used


def mean_enclosure(N, k=1):
    eps = Fr(1, 2 * (N + 1))
    S = Fr(0)
    mass = Fr(0)
    for name, n, P, A, C in pieces(N):
        mass += integrate_polygon(RHO, P)
        # split by signs of A, C, A-C; range = max(0,A,C) - min(0,A,C)
        AC = lin(A, C, 1, -1)
        for sA in (1, -1):
            for sC in (1, -1):
                for sD in (1, -1):
                    Q = P
                    for f, sg in ((A, sA), (C, sC), (AC, sD)):
                        Q = clip(Q, f if sg > 0 else neg(f))
                        if not Q:
                            break
                    if not Q:
                        continue
                    vals = {"0": (Fr(0), Fr(0), Fr(0)), "A": A, "C": C}
                    # determine max/min from signs
                    order_ok = True
                    # A>=0? C>=0? A>=C?
                    mx = "A" if (sA > 0 and sD > 0) else ("C" if (sC > 0 and sD < 0) else "0")
                    mn = "A" if (sA < 0 and sD < 0) else ("C" if (sC < 0 and sD > 0) else "0")
                    rng = lin(vals[mx], vals[mn], 1, -1)
                    tau = {kk: Fr(2, 3) * vv for kk, vv in aff(*rng).items()}
                    integrand = RHO
                    for _ in range(k):
                        integrand = pmul(integrand, tau)
                    S += integrate_polygon(integrand, Q)
    # omitted region delta<eps: tau in [1/3-delta, 1/3-2delta/3], P(delta<s)=R(s)
    # R'(s) = 40 s - 120 s^2 + 160 s^4
    def omit(c1, c2):  # int_0^eps (c1 + c2 s)^k R'(s) ds  exact
        tot = Fr(0)
        base = {0: Fr(c1), 1: Fr(c2)}
        p = {0: Fr(1)}
        for _ in range(k):
            q = {}
            for a, va in p.items():
                for b, vb in base.items():
                    q[a + b] = q.get(a + b, 0) + va * vb
            p = q
        Rp = {1: Fr(40), 2: Fr(-120), 4: Fr(160)}
        for a, va in p.items():
            for b, vb in Rp.items():
                tot += va * vb * eps ** (a + b + 1) / (a + b + 1)
        return tot
    lo = S + omit(Fr(1, 3), Fr(-1))
    hi = S + omit(Fr(1, 3), Fr(-2, 3))
    return lo, hi, mass, 1 - R(eps)


if __name__ == "__main__":
    import time
    for t in (Fr(1, 25), Fr(1, 10), Fr(1, 5), Fr(3, 10)):
        t0 = time.time()
        F, N, mass, expect_mass, nmax = cdf(t)
        print(f"F({t}) = {F} = {float(F):.12f}   [cutoff N={N}, retained mass "
              f"== 1-R(eps): {mass == expect_mass}, max shear count used {nmax}] "
              f"({time.time()-t0:.1f}s)")
        if t == Fr(1, 25):
            print(f"   paper (2): 713903/39062500; equal: {F == Fr(713903, 39062500)};"
                  f" tail 1-F = {1-F} = {float(1-F):.10f}")
    Nm = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    t0 = time.time()
    lo, hi, mass, em = mean_enclosure(Nm)
    print(f"E[L] enclosure with cutoff N={Nm}: [{float(lo):.12f}, {float(hi):.12f}] "
          f"width {float(hi-lo):.3e}; retained mass ok {mass == em} ({time.time()-t0:.1f}s)")
    plo, phi = Fr(1827011465, 10**10), Fr(1827011669, 10**10)
    print(f"   paper (5): ({float(plo)}, {float(phi)}); paper interval contains ours: "
          f"{plo < lo and hi < phi}; intervals overlap: {lo < phi and plo < hi}")
