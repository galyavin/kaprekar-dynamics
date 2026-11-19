#!/usr/bin/env python3
"""Test the four leading preterminal clock forms against
independent whole-word clocks.

For every lattice point of the closed rectangle at an even base (exhaustive),
T is computed by ./ww.  For each of the four forms f_k (with the certificate's
constant range [k_min, k_max] for the preterminal time), r_k = T - f_k(X).
A point is explained by form k when k_min < r_k <= k_max + 30*20 + 700
(30 returns of length <= 20 between the preterminal state and the exit,
plus a terminal tail).  We report how many points each form explains and
whether every point is explained by at least one form; the minimum leading
values are recomputed by exact corner enumeration.

Usage: forms_check.py WW B [B ...]
"""
import sys
from fractions import Fraction as F
from itertools import product

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from even_check import all_points, run_ww          # noqa: E402  (own code)

FORMS = [  # (row, k_min, k_max) from exclusion.txt
    ((F(145, 792), F(65, 72), F(-65, 132), F(-65, 72), F(1235, 792)), F(-16099, 33), F(-15289, 33)),
    ((F(125, 672), F(715, 672), F(-65, 112), F(-715, 672), F(1235, 672)), F(-36471, 56), F(-68501, 112)),
    ((F(115, 616), F(185, 168), F(-185, 308), F(-185, 168), F(3515, 1848)), F(-76877, 154), F(-143069, 308)),
    ((F(3, 16), F(55, 48), F(-5, 8), F(-55, 48), F(95, 48)), F(-2645, 4), F(-4917, 8)),
]
GS = (26293, 27557, 28481, 40447, 29865, 27357)
C = [F(g, 180000) for g in GS]
EPS = F(1, 10**6)


def corner_min(row):
    best = None
    for signs in product((-1, 1), repeat=5):
        y = [C[i] + signs[i] * EPS for i in range(5)]
        s, val = F(0), F(0)
        for i in range(5):
            s += y[i]
            val += row[i] * ((i + 1) - 6 * s)
        best = val if best is None or val < best else best
    return best


def main():
    ww = sys.argv[1]
    for k, (row, a, b) in enumerate(FORMS):
        print(f"form {k}: corner-enumerated min on closed rectangle = {corner_min(row)} = {float(corner_min(row)):.9f}")
    for B in map(int, sys.argv[2:]):
        pts = all_points(B)
        res = run_ww(ww, [(B, p) for p in pts], 4)
        counts = [0] * 4
        unexplained = 0
        multi = 0
        rmin = [None] * 4
        for _, G, T, _, _ in res:
            X = [(i + 1) * B - 6 * sum(G[:i + 1]) for i in range(5)]
            hits = []
            for k, (row, a, b) in enumerate(FORMS):
                r = T - sum(u * x for u, x in zip(row, X))
                if a < r <= b + 600 + 700:
                    hits.append(k)
                    rmin[k] = r if rmin[k] is None or r < rmin[k] else rmin[k]
            if not hits:
                unexplained += 1
            if len(hits) > 1:
                multi += 1
            for k in hits:
                counts[k] += 1
        print(f"B={B}: points={len(res)} explained per form={counts} unexplained={unexplained} "
              f"explained by >1 form={multi}")
        print("   min observed T - f_k(X) per form: " + ", ".join(
            f"{float(r):.2f} (claimed constant >= {float(FORMS[k][1]):.2f})" if r is not None else "-"
            for k, r in enumerate(rmin)))


if __name__ == "__main__":
    main()
