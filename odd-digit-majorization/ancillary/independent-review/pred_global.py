"""Predicted global statistics of mu_beta from the authors' supplement
(nine_global.leading_target + route domains + kernel table, scratch copy,
read-only), integrated by Monte Carlo over X ~ rho_9 (sorted differences of
nine i.i.d. uniform digits at base 10^7).

Outputs per beta: E[T], P(T<=t) on a grid, P(terminal = Q(4/5)), E[s_terminal].
usage: pred_global.py NSAMP SEED OUT.json
"""
import sys, json, random
from fractions import Fraction as Fr
from pathlib import Path

ANC = str(Path(__file__).resolve().parents[1] / "source")
sys.path.insert(0, ANC)
sys.dont_write_bytecode = True
from nine_global import leading_target  # noqa
from nine_entry import inverse_C, evaluate  # noqa
from nine_collar import route_table  # noqa

GRID = [0.05 * i for i in range(1, 21)]


def main():
    n = int(sys.argv[1])
    rnd = random.Random(int(sys.argv[2]))
    B0 = 10 ** 7
    domains, table = route_table()
    doms = [(tuple(s["route"]), s["domain"]) for s in domains]
    acc = {b: dict(ET=0.0, ET2=0.0, cdf=[0.0] * len(GRID), P45=0.0, Es=0.0) for b in range(5)}
    used = rej = 0
    routes_seen = {}
    prefix_hist = {}
    while used < n:
        d = sorted((rnd.randrange(B0) for _ in range(9)), reverse=True)
        x = tuple(Fr(d[i] - d[8 - i], B0) for i in range(4))
        try:
            y, prefix = leading_target(x, max_steps=20000)
        except Exception:
            rej += 1
            continue
        r = inverse_C(y)
        rs = [rt for rt, dom in doms if min(evaluate(g, (*r, 1)) for g in dom) > 0]
        if len(rs) != 1:
            rej += 1
            continue
        used += 1
        rt = rs[0]
        routes_seen[rt] = routes_seen.get(rt, 0) + 1
        prefix_hist[prefix] = prefix_hist.get(prefix, 0) + 1
        yf = (*y, 1)
        for b in range(5):
            rec = table[rt, b]
            for a in rec["atoms"]:
                w = float(Fr(a["weight"]))
                t = float(evaluate(tuple(map(Fr, a["clock"])), yf))
                z = [float(evaluate(tuple(map(Fr, g)), yf)) for g in a["terminal"]]
                s = z[0]
                A = acc[b]
                A["ET"] += w * t
                A["ET2"] += w * t * t
                for i, gq in enumerate(GRID):
                    if t <= gq:
                        A["cdf"][i] += w
                if abs(s - 0.8) < 1e-12:
                    A["P45"] += w
                A["Es"] += w * s
    out = {}
    for b in range(5):
        A = acc[b]
        out[b] = dict(ET=A["ET"] / used, ET2=A["ET2"] / used, cdf=[c / used for c in A["cdf"]], P45=A["P45"] / used,
                      Es=A["Es"] / used)
        print(f"beta={b}: E[T/B]={out[b]['ET']:.5f} P(Z=Q(4/5))={out[b]['P45']:.5f} E[s]={out[b]['Es']:.6f}", flush=True)
    print(f"used={used} rejected={rej}; prefix-length histogram (steps): {dict(sorted(prefix_hist.items())[:12])} ... max={max(prefix_hist)}")
    json.dump(dict(n=used, rejected=rej, grid=GRID, pred={str(b): v for b, v in out.items()},
                   routes={'/'.join(k): v for k, v in routes_seen.items()}), open(sys.argv[3], "w"), indent=1)


if __name__ == "__main__":
    main()
