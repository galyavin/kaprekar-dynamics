"""Lemma A.2 second half: terminal lifts and the discrepancy interval, own code.

For each cycle phase p=(m,S), B mod 4 = r4 (parity-compatible), j=1..D:
  start at u_end = ceil(3B0/4)+256-j, B0 = 100000+r4, follow F by integer subtraction
  (model.F_direct) until the first repetition; record h_end (first cyclic index).
Translation lemma certificate: under (q,B) -> (q+t(3,2,1), B+4t) the borrow digits move with
slopes (3,2,1,4,3,2,1); along every state of the path we check w>=1, that the slopes are
non-increasing in the (value, slope)-lexicographic order (so the sorted order persists for all
real t>=0), and that the symmetric slope differences are (3,2,1).
Also re-derives the discrepancy 3k + r e + h_end - r(u_end - 3B/4) over all initial offsets.
"""
from fractions import Fraction as Fr
import math, sys
import model
from edge import build_graph, analyse

SLOPES = (3, 2, 1, 4, 3, 2, 1)


def certify_translation(q, B):
    u, v, w = q
    if w < 1:
        return False
    d = model.borrow_digits(q, B)
    idx = sorted(range(7), key=lambda i: (-d[i], -SLOPES[i]))
    sl = [SLOPES[i] for i in idx]
    if any(sl[i] < sl[i + 1] for i in range(6)):
        return False
    return (sl[0] - sl[6], sl[1] - sl[5], sl[2] - sl[4]) == (3, 2, 1)


def terminal(p, r4, j, check_extra=()):
    m, S = p
    B0 = 100000 + r4
    HB = -((-3 * B0) // 4) + 256
    ue = HB - j
    q = (ue, (B0 + m) // 2, B0 - ue + S)
    assert 2 * q[1] - B0 == m and q[0] + q[2] - B0 == S
    path = [q]
    seen = {q: 0}
    while True:
        x = model.F_direct(path[-1], B0)
        assert x == model.F(path[-1], B0)
        if x in seen:
            hend = seen[x]
            break
        seen[x] = len(path)
        path.append(x)
    cert = all(certify_translation(x, B0) for x in path)
    cyc = path[hend:]
    cyc_ok = all(Fr(x[0]) - Fr(3 * B0, 4) < 128 for x in cyc)
    uoff = Fr(ue) - Fr(3 * B0, 4)
    # brute-force translation spot check at a few far bases
    spot = True
    for t in check_extra:
        B1 = B0 + 4 * t
        q1 = (q[0] + 3 * t, q[1] + 2 * t, q[2] + t)
        y = q1
        seen1 = {y: 0}
        k = 0
        while True:
            y = model.F(y, B1)
            k += 1
            if y in seen1:
                h1 = seen1[y]
                break
            seen1[y] = k
        spot &= (h1 == hend)
    return hend, cert, cyc_ok, uoff, spot, len(path)


if __name__ == "__main__":
    total_fam = 0
    allok = True
    for beta in (0, 1):
        start, rule, fails = build_graph(beta)
        cycles, maxk, maxdisp, maxpartial, entries = analyse(beta, start, rule)
        H = {}
        for key, (cyc, ell, D) in cycles.items():
            for p in cyc:
                for r4 in (beta, beta + 2):
                    for j in range(1, D + 1):
                        res = terminal(p, r4, j, check_extra=(1, 7, 250001) if j == 1 else ())
                        H[(p, r4, j)] = res
                        total_fam += 1
                        if not (res[1] and res[2] and res[4]):
                            allok = False
                            print("  FAIL family", beta, p, r4, j, res)
        lo, hi = None, None
        arg = None
        for (n0, k, e, key, entry) in entries:
            cyc, ell, D = cycles[key]
            r = Fr(3 * ell, D)
            for r4 in (beta, beta + 2):
                for j in range(1, D + 1):
                    hend, cert, cyc_ok, uoff, spot, L = H[(entry, r4, j)]
                    disc = 3 * k + r * e + hend - r * uoff
                    if lo is None or disc < lo:
                        lo = disc
                    if hi is None or disc > hi:
                        hi = disc
                        arg = (n0, k, e, entry, r4, j, hend, uoff)
        maxh = max(v[0] for v in H.values())
        print("beta=%d families=%d  translation certificates=%s  terminal cycles u-3B/4<128: %s  far-base spot checks: %s"
              % (beta, len(H), all(v[1] for v in H.values()), all(v[2] for v in H.values()), all(v[4] for v in H.values())))
        print("   discrepancy range [%s, %s] = [%.4f, %.4f]; max h_end=%d; argmax %s" % (lo, hi, float(lo), float(hi), maxh, arg))
    print("total endpoint families:", total_fam, "ALL OK" if allok else "FAILURES")
