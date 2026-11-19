#!/usr/bin/env python3
"""Spot-check the even exclusion tree.

Written AFTER the pre-read hashes (PRE-READ-SHA256) were recorded.

Terminal classes are drawn by random descent through the certificate's own tree
(their macro/remainders/repeat, imported read-only from a private snapshot):
at every node one remainder child is chosen uniformly, so classes are drawn
with probability equal to their lattice density (a random initial lattice
point's class).  exclusion.txt records only aggregate data, so the classes
themselves must be regenerated.

For each drawn class, independently of the certificate's arithmetic:
  * the leading minimum of its preterminal time form over the closed
    rectangle is recomputed by exact enumeration of the 32 corners, and
    compared with the certificate's interval and with 131/5000;
  * lattice points of the class inside the closed rectangle are found at
    two even bases; for each, ./spot (centered map) reports the first time
    the orbit visits the class's constant preterminal state, and the first
    cyclic time; ./ww (whole-word subtraction) reports T.  We check
      visit time == claimed preterminal time form (exact),
      visit time < first cyclic difference time,
      residual(preterminal state) > 504,
      T_ww == T_centered, T > claimed form, T/B > 131/5000, and
      T - form(X) equal at both points (constant terminal tail).

Usage: terminal_check.py SNAPSHOT_DIR SPOT WW N SEED [DUMP.jsonl.gz]
With DUMP (from dump_tree.py, a full rerun of the certificate's tree that records
every leaf), N classes are drawn uniformly from all leaves instead.
"""
import json
import random
import subprocess
import sys
from fractions import Fraction as F
from itertools import product
from math import ceil, floor

sys.path.insert(0, sys.argv[1])
import exclude as ex                       # noqa: E402  (primary certificate, read-only)
from affine_routes import (VARIABLES, ZERO, normalized_interval,   # noqa: E402
                           repeat)
from primitive import step                 # noqa: E402

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cohort_check import hnf_lower, coset_sample, box   # noqa: E402  (own code)

GS = (26293, 27557, 28481, 40447, 29865, 27357)
C = [F(g, 180000) for g in GS]
EPS = F(1, 10**6)


def descend(rng):
    state, elapsed = VARIABLES, ZERO
    x = (0,) * 5
    basis = tuple(tuple(6 * int(i == j) for j in range(5)) for i in range(5))
    hist = []
    for depth in range(ex.MAX_BLOCKS + 1):
        for _ in range(ex.MAX_SHORT):
            word, delta, modulus, num = ex.macro(state)
            if modulus:
                break
            state, elapsed = step(state, word[0])[0], elapsed + 1
        else:
            raise ArithmeticError("no return")
        children = list(ex.remainders(num, modulus, x, basis))
        rem, point, nb = rng.choice(children)
        assert ex.positive((num - rem) / modulus)
        nxt, time = repeat(state, elapsed, (word, delta, modulus, num), rem)
        hist.append((len(word), delta, len(children)))
        if all(not any(s.row) for s in nxt):
            endpoint = tuple(int(s.constant) for s in nxt)
            spot = tuple(e - ex.PRETERMINAL_RETURNS * d for e, d in zip(endpoint, delta))
            cut = time - ex.PRETERMINAL_RETURNS * len(word)
            return dict(point=point, basis=nb, cut=cut, endpoint=endpoint, spot=spot,
                        hist=hist, word_len=len(word), delta=delta)
        state, elapsed, x, basis = nxt, time, point, nb
    raise ArithmeticError("too deep")


def corner_min(row):
    """exact min over the closed y-box of row . X/B, X_i/B = i - 6 sum_{j<i} y_j."""
    best = None
    for signs in product((-1, 1), repeat=5):
        y = [C[i] + signs[i] * EPS for i in range(5)]
        val = F(0)
        s = F(0)
        for i in range(5):
            s += y[i]
            val += row[i] * ((i + 1) - 6 * s)
        best = val if best is None or val < best else best
    return best


def to_gaps(B, X):
    G = [F(B - X[0], 6)] + [F(B + X[i] - X[i + 1], 6) for i in range(4)]
    assert all(g.denominator == 1 for g in G)
    return [int(g) for g in G]


def gap_basis(basis):
    out = []
    for b in basis:
        v = [F(-b[0], 6)] + [F(b[i] - b[i + 1], 6) for i in range(4)]
        assert all(t.denominator == 1 for t in v)
        out.append([int(t) for t in v])
    return out


def residual(s):
    g = (-s[0], s[0] - s[1], s[1] - s[2], s[2] - s[3], s[3] - s[4], s[4])
    return -min(g)


def coset_sample_once(G0, H, lo, hi, rng):
    v = list(G0)
    for i in range(5):
        h = H[i][i]
        amin = ceil(F(lo[i] - v[i], h))
        amax = floor(F(hi[i] - v[i], h))
        if amax < amin:
            return None
        a = rng.randint(amin, amax)
        v = [p + a * q for p, q in zip(v, H[i])]
    return v


def run_par(binary, jobs, workers=5):
    """run line jobs in parallel chunks, preserving order."""
    from concurrent.futures import ThreadPoolExecutor
    idx = list(range(len(jobs)))
    chunks = [idx[k::workers] for k in range(workers)]

    def go(ch):
        if not ch:
            return []
        out = subprocess.run([binary], input="".join(jobs[i] for i in ch),
                             capture_output=True, text=True, check=True).stdout.splitlines()
        assert len(out) == len(ch)
        return list(zip(ch, out))

    res = [None] * len(jobs)
    with ThreadPoolExecutor(workers) as ex_:
        for part in ex_.map(go, chunks):
            for i, line in part:
                res[i] = line
    return res


def load_dump(path, n, rng):
    """uniform sample of n leaf classes from a full-tree dump (dump_tree.py)."""
    import gzip
    from affine_routes import Form
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    out = []
    for line in rng.sample(lines, n):
        d = json.loads(line)
        cut = Form(tuple(F(x) for x in d["cut"][:5]), F(d["cut"][5]))
        out.append(dict(point=tuple(d["point"]), basis=tuple(tuple(b) for b in d["basis"]), cut=cut,
                        endpoint=tuple(d["endpoint"]), spot=tuple(d["spot"]), hist=[(h[0], tuple(h[1]), 0) for h in d["hist"]],
                        word_len=d["word_len"], delta=tuple(d["delta"])))
    return out, len(lines)


def main():
    snap, spot_bin, ww_bin, n, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    rng = random.Random(seed)
    if len(sys.argv) > 6:
        classes, total = load_dump(sys.argv[6], n, rng)
        print(f"uniform sample of {n} classes from full-tree dump with {total} leaves")
    else:
        classes = [descend(rng) for _ in range(n)]
        print(f"{n} classes drawn by density-weighted random descent")
    jobs_spot, jobs_ww, meta = [], [], []
    rows_seen = {}
    for k, cls in enumerate(classes):
        row = cls["cut"].row
        cm = corner_min(row)
        lo_t, _ = normalized_interval(row)
        rows_seen.setdefault(tuple(row), [cm, lo_t, 0])
        rows_seen[tuple(row)][2] += 1
        H = hnf_lower(gap_basis(cls["basis"]))
        index = 1
        for i in range(5):
            index *= H[i][i]
        hmax = max(H[i][i] for i in range(5))
        need = max(F(2 * hmax), F(4) * F(index) ** F(1, 5) if index < 10**40 else F(10**9)) / (2 * EPS)
        B1 = 6 * ceil(max(need, F(10**7)) / 6)
        pts = []
        for B in (B1, 6 * ceil(F(17, 10) * B1 / 6)):
            lo, hi = box(B, EPS)
            p = None
            for _ in range(200):
                p = coset_sample_once(to_gaps(B, cls["point"]), H, lo, hi, rng)
                if p is not None:
                    break
            if p is None:
                pts.append(None)
                continue
            X = [(i + 1) * B - 6 * sum(p[:i + 1]) for i in range(5)]
            cutv = cls["cut"].evaluate(X)
            assert cutv.denominator == 1
            pts.append((B, p, int(cutv)))
            jobs_spot.append(f"{B} {' '.join(map(str, p))} {' '.join(map(str, cls['spot']))}\n")
            jobs_ww.append(f"{B} {' '.join(map(str, p))}\n")
        meta.append((cls, cm, lo_t, index, pts))
    out_s = run_par(spot_bin, jobs_spot)
    out_w = run_par(ww_bin, jobs_ww)
    js = 0
    tot = dict(classes=n, points=0, missing=0, visit_eq=0, before_cycle=0, resid_ok=0, ww_eq=0,
               above_form=0, above_ramp=0, tail_const=0, tail_pairs=0, cm_eq=0, cm_ramp=0)
    worst = None
    for cls, cm, lo_t, index, pts in meta:
        tot["cm_eq"] += (cm == lo_t)
        tot["cm_ramp"] += (cm > F(131, 5000))
        tails = []
        for p in pts:
            if p is None:
                tot["missing"] += 1
                continue
            B, G, cutv = p
            fs = list(map(int, out_s[js].split()))
            fw = list(map(int, out_w[js].split()))
            js += 1
            hit, mu, Tc = fs[6], fs[7], fs[8]
            Tw = fw[6]
            tot["points"] += 1
            tot["visit_eq"] += (hit == cutv)
            tot["before_cycle"] += (0 <= hit < mu)
            tot["resid_ok"] += residual(cls["spot"]) > 504
            tot["ww_eq"] += (Tw == Tc)
            tot["above_form"] += (Tw > cutv)
            tot["above_ramp"] += (F(Tw, B) > F(131, 5000))
            tails.append(Tw - cutv)
            slack = F(Tw, B) - F(131, 5000)
            if worst is None or slack < worst[0]:
                worst = (slack, B, G, Tw, cutv)
            print(json.dumps(dict(B=B, G=G, index=index, hist=[(h[0], list(h[1]), h[2]) for h in cls["hist"]],
                                  spot=cls["spot"], form=[str(r) for r in cls["cut"].row] + [str(cls["cut"].constant)],
                                  claimed=cutv, visit=hit, mu=mu, T_ww=Tw, T_cen=Tc,
                                  corner_min=str(cm), certificate_lower=str(lo_t))))
        if len(tails) == 2:
            tot["tail_pairs"] += 1
            tot["tail_const"] += (tails[0] == tails[1])
    print("SUMMARY", json.dumps(tot))
    print("worst T/B - 131/5000:", float(worst[0]), "at B", worst[1], "G", worst[2], "T", worst[3], "claimed form", worst[4])
    for row, (cm, lo_t, cnt) in sorted(rows_seen.items(), key=lambda r: r[1][0]):
        print("form", [str(r) for r in row], "drawn", cnt, "corner min", cm, float(cm), "certificate", lo_t, "equal", cm == lo_t)


if __name__ == "__main__":
    main()
