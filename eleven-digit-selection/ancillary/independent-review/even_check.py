#!/usr/bin/env python3
"""Independent even-arm exclusion check.

Written from proof.md definitions only.  For each chosen even base B, list
lattice points G (integer, G5 = B - sum) of the CLOSED rectangle
|G_i/B - c_i| <= 1e-6 (i < 5), compute the whole-word first cyclic word
time T_B with ./ww, and compare exactly (rationals) with
    131/5000   and   9451/337500 - 16099/(33 B).

Usage: even_check.py WW_BINARY MODE [args]   (modes below)
"""
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from math import ceil, floor

GS = (26293, 27557, 28481, 40447, 29865, 27357)
C = [F(g, 180000) for g in GS]
EPS = F(1, 10**6)
RAMP_HI = F(131, 5000)
LEAD = F(9451, 337500)
OFF = F(16099, 33)


def box(B, eps=EPS):
    lo = [ceil((C[i] - eps) * B) for i in range(5)]
    hi = [floor((C[i] + eps) * B) for i in range(5)]
    return lo, hi


def all_points(B):
    lo, hi = box(B)
    pts = [[]]
    for i in range(5):
        if hi[i] < lo[i]:
            return []
        pts = [p + [g] for p in pts for g in range(lo[i], hi[i] + 1)]
    return pts


def corners(B):
    lo, hi = box(B)
    if any(hi[i] < lo[i] for i in range(5)):
        return []
    pts = [[]]
    for i in range(5):
        pts = [p + [g] for p in pts for g in sorted({lo[i], hi[i]})]
    return pts


def random_points(B, n, rng):
    lo, hi = box(B)
    return [[rng.randint(lo[i], hi[i]) for i in range(5)] for _ in range(n)]


def run_ww(ww, jobs, workers):
    """jobs: list of (B, G[5]); returns list of (B, G, T, MUD, LAM)."""
    chunks = [jobs[k::workers] for k in range(workers)]

    def go(chunk):
        if not chunk:
            return []
        inp = "".join(f"{B} {' '.join(map(str, G))}\n" for B, G in chunk)
        out = subprocess.run([ww], input=inp, capture_output=True, text=True, check=True).stdout
        res = []
        for line in out.splitlines():
            f = line.split()
            if "BAD" in f:
                raise SystemExit("BAD input " + line)
            v = list(map(int, f))
            res.append((v[0], v[1:6], v[6], v[7], v[8]))
        if len(res) != len(chunk):
            raise SystemExit("missing output")
        return res

    with ThreadPoolExecutor(workers) as ex:
        out = []
        for r in ex.map(go, chunks):
            out.extend(r)
    return out


def analyse(results, tag):
    """Return summary dict; exact comparisons."""
    n = len(results)
    fails_ramp = fails_bound = mud_bad = 0
    min_tau = None
    min_slack = None
    arg_tau = arg_slack = None
    lanes = {}
    for B, G, T, MUD, LAM in results:
        if MUD != T - 1:
            mud_bad += 1
        tau = F(T, B)
        bound = LEAD - OFF / B
        if not tau > RAMP_HI:
            fails_ramp += 1
        if not tau > bound:
            fails_bound += 1
        slack = T - (LEAD * B - OFF)          # in steps; bound needs slack > 0
        if min_tau is None or tau < min_tau:
            min_tau, arg_tau = tau, (B, G, T)
        if min_slack is None or slack < min_slack:
            min_slack, arg_slack = slack, (B, G, T)
        key = round(float(tau), 4)
        lanes[key] = lanes.get(key, 0) + 1
    return dict(tag=tag, n=n, fails_ramp=fails_ramp, fails_bound=fails_bound,
                mud_bad=mud_bad, min_tau=min_tau, arg_tau=arg_tau,
                min_slack=min_slack, arg_slack=arg_slack, lanes=lanes)


def report(s):
    if s["n"] == 0:
        print(f"{s['tag']}: no lattice points")
        return
    print(f"{s['tag']}: points={s['n']} fail(T/B<=131/5000)={s['fails_ramp']} "
          f"fail(T/B<=bound)={s['fails_bound']} MUD!=T-1:{s['mud_bad']}")
    print(f"   min T/B={float(s['min_tau']):.9f} at B={s['arg_tau'][0]} G={s['arg_tau'][1]} T={s['arg_tau'][2]}")
    print(f"   min slack T-(9451/337500 B-16099/33)={float(s['min_slack']):.3f} at B={s['arg_slack'][0]} "
          f"G={s['arg_slack'][1]} T={s['arg_slack'][2]}")
    print("   T/B histogram (4 d.p.): " + ", ".join(f"{k:.4f}:{v}" for k, v in sorted(s['lanes'].items())))
    sys.stdout.flush()


def main():
    ww, mode = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3])
    rng = random.Random(20261118)
    if mode == "threshold":                     # all B = 0 mod 6 in [lo, hi], every lattice point
        lo, hi = int(sys.argv[4]), int(sys.argv[5])
        b0 = lo + (-lo) % 6
        jobs = []
        nb = nonempty = 0
        for B in range(b0, hi + 1, 6):
            pts = all_points(B)
            nb += 1
            nonempty += bool(pts)
            jobs.extend((B, p) for p in pts)
        print(f"threshold scan B=0 mod 6 in [{lo},{hi}]: bases={nb} nonempty={nonempty} points={len(jobs)}")
        res = run_ww(ww, jobs, workers)
        report(analyse(res, f"all B=0 mod 6 in [{lo},{hi}]"))
        prog = [r for r in res if r[0] % 3840 == 3360]
        report(analyse(prog, f"  of which B=3360 mod 3840"))
        # smallest bases with points
        firsts = sorted({r[0] for r in res})[:5]
        for B in firsts:
            report(analyse([r for r in res if r[0] == B], f"  B={B}"))
    elif mode == "exhaustive":                  # listed bases, every lattice point
        for B in map(int, sys.argv[4:]):
            pts = all_points(B)
            res = run_ww(ww, [(B, p) for p in pts], workers)
            report(analyse(res, f"exhaustive B={B} (B mod 3840 = {B % 3840}, mod 6 = {B % 6})"))
    elif mode == "sample":                      # listed B:n:c ; n random points (+ all 32 corners if c=1)
        for spec in sys.argv[4:]:
            B, n, c = map(int, spec.split(":"))
            pts = (corners(B) if c else []) + random_points(B, n, rng)
            res = run_ww(ww, [(B, p) for p in pts], workers)
            report(analyse(res, f"sample B={B} (mod 3840 = {B % 3840}, mod 6 = {B % 6}) "
                                f"{'32 corners + ' if c else ''}{n} uniform random lattice points"))
    else:
        raise SystemExit("mode?")


if __name__ == "__main__":
    main()
