"""Compare uniform-string Monte Carlo at fresh bases with the predicted limits
(merged pred_part*.json), and test 5-periodicity across bases."""
import glob, json, math, sys
import numpy as np
from scipy.stats import ks_2samp

HERE = __file__.rsplit("/", 1)[0]
parts = [json.load(open(f)) for f in sorted(glob.glob(HERE + "/pred_part*.json"))]
N = sum(p["n"] for p in parts)
GRID = parts[0]["grid"]
pred = {}
for b in range(5):
    pred[b] = {k: sum(p["pred"][str(b)][k] * p["n"] for p in parts) / N for k in ("ET", "ET2", "P45", "Es")}
    pred[b]["cdf"] = [sum(p["pred"][str(b)]["cdf"][i] * p["n"] for p in parts) / N for i in range(len(GRID))]
# MC standard error of predicted P45 and ET from the 8 independent parts
for b in range(5):
    v = [p["pred"][str(b)]["P45"] for p in parts]
    e = [p["pred"][str(b)]["ET"] for p in parts]
    pred[b]["P45_se"] = np.std(v, ddof=1) / math.sqrt(len(v))
    pred[b]["ET_se"] = np.std(e, ddof=1) / math.sqrt(len(e))
print(f"prediction: {N} limit-law samples from the supplement's kernels")
for b in range(5):
    print(f"  beta={b}: P(Z=Q(4/5))={pred[b]['P45']:.4f}+-{pred[b]['P45_se']:.4f}  E[T/B]={pred[b]['ET']:.5f}+-{pred[b]['ET_se']:.5f}")
# predicted class differences with common random numbers
for b in range(1, 5):
    d = [p["pred"][str(b)]["P45"] - p["pred"]["0"]["P45"] for p in parts]
    print(f"  pred P45(beta={b})-P45(0) = {np.mean(d):+.4f} +- {np.std(d,ddof=1)/math.sqrt(len(d)):.4f}")


def seg_s(Zn):
    a = np.array([1.0, 2.0, 3.0, -1.0])
    bb = np.array([0.0, -1.0, -2.0, 1.0])
    s = ((Zn - bb) * a).sum(axis=1) / float((a * a).sum())
    sc = np.clip(s, 0.8, 5 / 6)
    P = sc[:, None] * a[None, :] + bb[None, :]
    return sc, np.max(np.abs(Zn - P), axis=1)


obs = {}
rows = []
for f in sorted(glob.glob(HERE + "/uniform_B*.npy")):
    B = int(f.split("_B")[1].split(".")[0])
    A = np.load(f)
    T = A[:, 0] / B
    s, dist = seg_s(A[:, 1:5] / B)
    p45 = np.mean(np.abs(s - 0.8) * B <= 3)
    n = len(T)
    se = math.sqrt(p45 * (1 - p45) / n)
    b = B % 5
    zP = (p45 - pred[b]["P45"]) / math.sqrt(se ** 2 + pred[b]["P45_se"] ** 2)
    zT = (T.mean() - pred[b]["ET"]) / math.sqrt(T.var() / n + pred[b]["ET_se"] ** 2)
    cdf = [np.mean(T <= g) for g in GRID]
    ks = max(abs(c - pc) for c, pc in zip(cdf, pred[b]["cdf"]))
    obs[B] = (T, p45, n)
    rows.append((b, B, n, p45, se, pred[b]["P45"], zP, T.mean(), pred[b]["ET"], zT, ks, dist.max() * B))
print("\nbeta  B       n       obs P45  (se)     pred P45  z      obs E[T/B]  pred E[T/B]  z      sup|CDF-pred|  max segdist*B")
for r in sorted(rows):
    print(f"{r[0]}  {r[1]:6d} {r[2]:7d}  {r[3]:.4f} ({r[4]:.4f})  {r[5]:.4f}  {r[6]:+.2f}  {r[7]:.5f}     {r[8]:.5f}     {r[9]:+.2f}  {r[10]:.4f}         {r[11]:.1f}")
print("\nsame residue, different parity/mod 3 (should agree):")
Bs = sorted(obs)
for B in Bs:
    for B2 in Bs:
        if B2 > B and (B2 - B) % 5 == 0:
            T1, p1, n1 = obs[B]
            T2, p2, n2 = obs[B2]
            z = (p1 - p2) / math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
            ks = ks_2samp(T1, T2)
            print(f"  B={B} vs {B2}: dP45={p1-p2:+.4f} z={z:+.2f}; KS(T/B) D={ks.statistic:.4f} p={ks.pvalue:.3f}")
print("\ndifferent residues (pooled per class):")
pool = {}
for B, (T, p, n) in obs.items():
    pool.setdefault(B % 5, []).append((p, n))
pp = {b: (sum(p * n for p, n in v) / sum(n for p, n in v), sum(n for p, n in v)) for b, v in pool.items()}
for b in range(5):
    for b2 in range(b + 1, 5):
        (p1, n1), (p2, n2) = pp[b], pp[b2]
        z = (p1 - p2) / math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
        print(f"  beta {b} vs {b2}: obs dP45={p1-p2:+.4f} (z={z:+.1f})  pred dP45={pred[b]['P45']-pred[b2]['P45']:+.4f}")
