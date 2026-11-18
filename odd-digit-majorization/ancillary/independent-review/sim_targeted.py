"""Targeted replication of selected kernels (same construction as
sim_allroutes.py: Chebyshev centre, certified box, 16 branches), with an
explicit base and a fresh seed.  Output uses the same raw format, so
analyze_allroutes.py applies unchanged.

usage: sim_targeted.py OUT.npz SEED N "route|beta|B" ["route|beta|B" ...]
"""
import sys, json, random
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import sim_allroutes as S


def main():
    outp, seed, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    specs = [a.split("|") for a in sys.argv[4:]]
    doms = S.route_domains()
    rnames = [tuple(d["route"]) for d in doms]
    rnd = random.Random(seed)
    jobs = []
    for rname, beta, B in specs:
        route = tuple(rname.split("/"))
        beta, B = int(beta), int(B)
        assert B % 5 == beta
        d = doms[rnames.index(route)]
        rows = [tuple(S.Fr(t) for t in g) for g in d["domain"]]
        rc, m = S.cheb_center(rows)
        R = int(0.9 * float(m) * B / 2)
        for _ in range(n):
            sg = [rnd.choice((1, -1)) for _ in range(4)]
            x = sorted([(1 + s * float(t)) / 2 for s, t in zip(sg, rc)], reverse=True)
            q = [int(round(xi * B)) + rnd.randint(-R, R) for xi in x]
            jobs.append((B, q, (route, beta, R)))
    k = 8
    chunks = [jobs[i::k] for i in range(k)]
    with ThreadPoolExecutor(k) as ex:
        outs = list(ex.map(S.run_chunk, chunks))
    out = []
    for ch, oc in zip(chunks, outs):
        for (B, q, (route, beta, R)), o in zip(ch, oc):
            out.append([rnames.index(route), beta, B, R, *q, *o])
    np.savez_compressed(outp, rows=np.array(out, dtype=np.int64), routes=np.array(['/'.join(r) for r in rnames]),
                        meta=json.dumps(dict(specs=specs, seed=seed, n=n)))
    print(f"saved {len(out)} samples to {outp}", flush=True)


if __name__ == "__main__":
    main()
