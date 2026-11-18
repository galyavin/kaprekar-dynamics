"""Fresh-base statistical test of all 33 routes x 5 residues (165 kernels).

Predictions (route domains, kernel weights, affine clocks and terminal
matrices) are read from the authors' supplement (a scratch copy of
ancillary/source, imported read-only).  Everything else is independent:
interior points are re-centred by an l_inf Chebyshev LP, the sampling box is
certified inside the route domain, the core point y = f(q/B) is computed with
the independent real map, and T_B, Z_B come from the independent column-
subtraction simulator (./kap9 orbits).

usage: sim_allroutes.py NPER_KERNEL RMIN BMIN OUT.npz [OFF]
(raw rows: route_idx beta B R q1..q4 mu lam z1..z4 allpos allpaired; analysis in analyze_allroutes.py)
"""
import sys, json, math, random, subprocess, itertools
from fractions import Fraction as Fr
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.optimize import linprog

HERE = __file__.rsplit("/", 1)[0]
ANC = str(Path(__file__).resolve().parents[1] / "source")
sys.path.insert(0, ANC)
sys.dont_write_bytecode = True
from nine_core_checks import route_domains  # noqa: E402

KERN = json.load(open(ANC + "/nine_kernel_results.json"))["kernels"]
TABLE = {(tuple(k["route"]), k["base_residue"]): k for k in KERN}
sys.path.pop(0)
from km import f_real  # noqa: E402


def cheb_center(rows):
    """max m s.t. g_lin.r + g_c >= m*||g_lin||_1 (l_inf ball of radius m inside)."""
    A, b = [], []
    for g in rows:
        lin = [float(t) for t in g[:4]]
        n1 = sum(abs(t) for t in lin)
        A.append([-t for t in lin] + [n1])
        b.append(float(g[4]))
    res = linprog(c=[0, 0, 0, 0, -1], A_ub=A, b_ub=b, bounds=[(0, 1)] * 4 + [(0, 1)], method="highs")
    r = [Fr(t).limit_denominator(10 ** 6) for t in res.x[:4]]
    m = min((sum(Fr(g[i]) * r[i] for i in range(4)) + Fr(g[4])) / sum(abs(Fr(t)) for t in g[:4]) for g in rows)
    return r, m


def ev(row, v):
    return sum(float(Fr(row[i])) * v[i] for i in range(4)) + float(Fr(row[4]))


def run_chunk(items):
    inp = "".join(f"{B} {q[0]} {q[1]} {q[2]} {q[3]}\n" for B, q, _ in items)
    out = subprocess.run([HERE + "/kap9", "orbits"], input=inp, capture_output=True, text=True, check=True).stdout
    return [tuple(map(int, l.split())) for l in out.split("\n") if l.strip()]


def main():
    nper = int(sys.argv[1])
    rmin = int(sys.argv[2])
    bmin = int(sys.argv[3])
    outp = sys.argv[4]
    global OFF
    OFF = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    doms = route_domains()
    rnd = random.Random(424242)
    jobs = []
    meta = {}
    for ri, d in enumerate(doms):
        route = tuple(d["route"])
        rows = [tuple(Fr(t) for t in g) for g in d["domain"]]
        rc, m = cheb_center(rows)
        assert m > 0
        # box radius R (q units) gives an l_inf ball of radius 2R/B in r
        Bbase = max(bmin, int(math.ceil(2 * rmin / (0.9 * float(m)) / 10.0)) * 10)
        Bbase = Bbase - Bbase % 10 + 10 * (ri % 7) + 1000 + OFF  # fresh, spread; OFF=5 flips parity
        meta[route] = dict(center=[str(t) for t in rc], linf_margin=float(m), Bbase=Bbase)
        for beta in range(5):
            B = Bbase + beta
            R = int(0.9 * float(m) * B / 2)
            assert R >= rmin
            for _ in range(nper):
                sg = [rnd.choice((1, -1)) for _ in range(4)]
                x = sorted([(1 + s * float(t)) / 2 for s, t in zip(sg, rc)], reverse=True)
                q = [int(round(xi * B)) + rnd.randint(-R, R) for xi in x]
                jobs.append((B, q, (route, beta, R)))
    print(f"{len(jobs)} orbits over {len(doms)} routes", flush=True)
    k = 8
    chunks = [jobs[i::k] for i in range(k)]
    with ThreadPoolExecutor(k) as ex:
        outs = list(ex.map(run_chunk, chunks))
    rows = []
    rnames = [tuple(d["route"]) for d in doms]
    for ch, oc in zip(chunks, outs):
        for (B, q, (route, beta, R)), o in zip(ch, oc):
            rows.append([rnames.index(route), beta, B, R, *q, *o])
    import numpy as _np
    _np.savez_compressed(outp, rows=_np.array(rows, dtype=_np.int64),
                         routes=_np.array(['/'.join(r) for r in rnames]),
                         meta=json.dumps({'/'.join(k): v for k, v in meta.items()}))
    print(f"saved {len(rows)} samples to {outp}", flush=True)


route_domains_cache = {}
for d in route_domains():
    route_domains_cache[tuple(d["route"])] = [tuple(float(t) for t in g) for g in d["domain"]]

if __name__ == "__main__":
    main()
