"""Analyse raw all-route samples against the supplement's 165 kernels.

For every sample: y = f(q/B) by the independent real map; r = C^{-1}(y) must
lie strictly in the route domain; observed (T/B, Z/B) with T = mu+1 and Z the
successor of the first cyclic difference.  For each atom l of the kernel the
error e_l = max(|T/B - t_l(y)|, max_j |Z_j/B - U_l(y)_j|).

Atoms are resolved only up to the finite-base error bound E = CERR/B:
  * a sample is MATCHED if some atom has e_l <= E, else UNMATCHED (reported);
  * atoms l, l' are merged into one group whenever, at some sample of the
    kernel, both are within E of the observation or their predicted
    observables are within 2E of each other.  Distinct groups are therefore
    > 2E apart at every sample, so the group of a matched sample is unique.
Group frequencies are compared with summed weights (chi-square).

usage: analyze_allroutes.py RAW.npz [RAW2.npz] [CERR]
"""
import sys, json, math
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np
from scipy.stats import chi2 as C2

HERE = __file__.rsplit("/", 1)[0]
ANC = str(Path(__file__).resolve().parents[1] / "source")
sys.path.insert(0, ANC)
sys.dont_write_bytecode = True
from nine_core_checks import route_domains  # noqa
KERN = json.load(open(ANC + "/nine_kernel_results.json"))["kernels"]
sys.path.pop(0)

TABLE = {('/'.join(k["route"]), k["base_residue"]): k for k in KERN}
DOM = {'/'.join(d["route"]): np.array([[float(Fr(t)) for t in g] for g in d["domain"]]) for d in route_domains()}


def f_real_np(x):
    """vectorised leading map: x (n,4) -> y (n,4)"""
    n = x.shape[0]
    vals = np.concatenate([np.ones((n, 1)), x, 1 - x], axis=1)
    d = -np.sort(-vals, axis=1)
    return np.stack([d[:, i] - d[:, 8 - i] for i in range(4)], axis=1)


def analyse(files, cerr):
    rows = []
    for f in files:
        D = np.load(f, allow_pickle=False)
        names = list(D["routes"])
        R = D["rows"]
        rows.append((names, R, f))
    results = {}
    for names, R, f in rows:
        for ri, rname in enumerate(names):
            for beta in range(5):
                sel = R[(R[:, 0] == ri) & (R[:, 1] == beta)]
                if len(sel) == 0:
                    continue
                B = sel[:, 2].astype(float)
                q = sel[:, 4:8].astype(float)
                mu = sel[:, 8]
                z = sel[:, 10:14].astype(float)
                paired = sel[:, 15]
                y = f_real_np(q / B[:, None])
                a = 2 * y[:, 0] - 1
                b = 2 * y[:, 1] - a
                c = 2 * y[:, 2] - b
                dd = 2 * y[:, 3] - c
                rr = np.stack([a, b, c, dd, np.ones_like(a)], axis=1)
                inside = (rr @ DOM[rname].T).min(axis=1) > 0
                kern = TABLE[(rname, beta)]
                yb = np.concatenate([y, np.ones((len(y), 1))], axis=1)
                tpred = np.stack([yb @ np.array([float(Fr(t)) for t in at["clock"]]) for at in kern["atoms"]], axis=1)
                zpred = np.stack([yb @ np.array([[float(Fr(t)) for t in g] for g in at["terminal"]]).T
                                  for at in kern["atoms"]], axis=1)  # (n, L, 4)
                tobs = (mu + 1) / B
                zobs = z / B[:, None]
                err = np.maximum(np.abs(tobs[:, None] - tpred), np.abs(zobs[:, None, :] - zpred).max(axis=2))
                E = cerr / B
                L = err.shape[1]
                # union-find
                parent = list(range(L))

                def find(i):
                    while parent[i] != i:
                        parent[i] = parent[parent[i]]
                        i = parent[i]
                    return i

                within = err <= E[:, None]
                for i in range(L):
                    for j in range(i + 1, L):
                        sep = np.maximum(np.abs(tpred[:, i] - tpred[:, j]), np.abs(zpred[:, i] - zpred[:, j]).max(axis=1))
                        if np.any((within[:, i] & within[:, j]) | (sep <= 2 * E)):
                            parent[find(i)] = find(j)
                groups = sorted({find(i) for i in range(L)})
                gidx = {g: k for k, g in enumerate(groups)}
                matched = within.any(axis=1) & inside
                best = err.argmin(axis=1)
                counts = np.zeros(len(groups), dtype=int)
                for s in np.where(matched)[0]:
                    counts[gidx[find(best[s])]] += 1
                w = np.zeros(len(groups))
                for l, at in enumerate(kern["atoms"]):
                    w[gidx[find(l)]] += float(Fr(at["weight"]))
                key = (rname, beta)
                r_ = results.setdefault(key, dict(counts=np.zeros(len(groups), dtype=int), w=w, n_unmatched=0, n_outside=0,
                                                  n_notpaired=0, groups=[[l for l in range(L) if find(l) == g] for g in groups],
                                                  natoms=L, maxerr=0.0, B=[]))
                if len(r_["counts"]) != len(counts):  # different grouping in another file: keep the coarser one per file
                    r_.setdefault("extra", []).append((f, counts.tolist(), w.tolist()))
                else:
                    r_["counts"] += counts
                r_["n_unmatched"] += int((~within.any(axis=1) & inside).sum())
                r_["n_outside"] += int((~inside).sum())
                r_["n_notpaired"] += int((paired == 0).sum())
                r_["maxerr"] = max(r_["maxerr"], float((err.min(axis=1)[inside] * B[inside]).max()))
                r_["B"].append(int(B[0]))
    return results


def report(results, label):
    tot = 0.0
    df = 0
    out = []
    worst = []
    for (rname, beta), r_ in sorted(results.items()):
        n = int(r_["counts"].sum())
        w = r_["w"]
        dfk = len(w) - 1
        chi = float(((r_["counts"] - n * w) ** 2 / (n * w)).sum()) if dfk > 0 else 0.0
        tot += chi
        df += dfk
        p = C2.sf(chi, dfk) if dfk > 0 else 1.0
        worst.append((p, rname, beta, r_["counts"].tolist(), [round(x, 4) for x in w], r_["groups"]))
        out.append(dict(route=rname, beta=beta, B=r_["B"], n=n, natoms=r_["natoms"], groups=r_["groups"],
                        counts=r_["counts"].tolist(), weights=w.tolist(), chi2=chi, df=dfk, p=p,
                        unmatched=r_["n_unmatched"], outside=r_["n_outside"], notpaired=r_["n_notpaired"],
                        max_err_steps=r_["maxerr"], extra=r_.get("extra")))
    nun = sum(o["unmatched"] for o in out)
    nout = sum(o["outside"] for o in out)
    npd = sum(o["notpaired"] for o in out)
    ntot = sum(o["n"] for o in out)
    merged = sum(1 for o in out if len(o["groups"]) < o["natoms"])
    print(f"[{label}] kernels={len(out)} matched={ntot} unmatched={nun} outside-box={nout} not-paired={npd} "
          f"max matched error={max(o['max_err_steps'] for o in out):.1f} steps; kernels with merged atoms={merged}")
    print(f"[{label}] total chi2={tot:.1f} df={df} p={C2.sf(tot, df):.3f}")
    worst.sort()
    for p, rn, b, c, w, g in worst[:5]:
        print(f"   lowest p={p:.4f} {rn} beta={b} counts={c} weights={w} groups={g}")
    return out, tot, df


if __name__ == "__main__":
    args = sys.argv[1:]
    cerr = 120.0
    if args and not args[-1].endswith(".npz"):
        cerr = float(args.pop())
    allout = {}
    for f in args:
        res = analyse([f], cerr)
        out, tot, df = report(res, f.rsplit("/", 1)[-1])
        allout[f] = dict(kernels=out, chi2=tot, df=df, cerr=cerr)
    if len(args) == 2:
        # homogeneity between the two base families, on the coarser common grouping
        A = {(o["route"], o["beta"]): o for o in allout[args[0]]["kernels"]}
        Bk = {(o["route"], o["beta"]): o for o in allout[args[1]]["kernels"]}
        ch = 0.0
        d = 0
        for k in A:
            a, b = A[k], Bk[k]
            if a["groups"] != b["groups"] or a["df"] == 0:
                continue
            na, nb = a["n"], b["n"]
            for ca, cb in zip(a["counts"], b["counts"]):
                t = ca + cb
                if t:
                    ea, eb = t * na / (na + nb), t * nb / (na + nb)
                    ch += (ca - ea) ** 2 / ea + (cb - eb) ** 2 / eb
            d += a["df"]
        print(f"homogeneity of the two base families (same B mod 5, opposite parity): chi2={ch:.1f} df={d} p={C2.sf(ch, d):.3f}")
        allout["homogeneity"] = dict(chi2=ch, df=d)
    json.dump(allout, open(HERE + "/allroutes_analysis.json", "w"), indent=1, default=str)
