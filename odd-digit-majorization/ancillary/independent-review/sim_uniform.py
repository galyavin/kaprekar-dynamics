"""Fresh-base Monte Carlo of the nine-digit joint law for uniformly random
non-constant strings N (direct digit subtraction, ./kap9 uniform).

Per base: T_B/B summary, terminal Z_B/B: distance to the segment
Q(s)=(s,2s-1,3s-2,1-s), 4/5<=s<=5/6, the mass at Q(4/5) (|s-4/5|*B<=3),
and the full samples saved for cross-base comparisons.

usage: sim_uniform.py NSAMP B1 [B2 ...]
"""
import sys, subprocess, json, math
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = __file__.rsplit("/", 1)[0]


def run(B, n, seed):
    out = subprocess.run([HERE + "/kap9", "uniform", str(B), str(n), str(seed)], capture_output=True, text=True,
                         check=True).stdout
    return np.array([list(map(int, l.split())) for l in out.split("\n") if l.strip()], dtype=np.int64)


def seg_dist(Zn):
    """sup-distance of normalised terminal points to the segment Q([4/5,5/6]);
    least squares s then clip."""
    a = np.array([1.0, 2.0, 3.0, -1.0])
    b = np.array([0.0, -1.0, -2.0, 1.0])
    s = ((Zn - b) * a).sum(axis=1) / float((a * a).sum())
    s = np.clip(s, 0.8, 5 / 6)
    P = s[:, None] * a[None, :] + b[None, :]
    return np.max(np.abs(Zn - P), axis=1), s


def main():
    n = int(sys.argv[1])
    bases = [int(b) for b in sys.argv[2:]]
    k = 8
    for B in bases:
        with ThreadPoolExecutor(k) as ex:
            parts = list(ex.map(lambda i: run(B, n // k, 1000 * B + i), range(k)))
        A = np.vstack(parts)
        np.save(HERE + f"/uniform_B{B}.npy", A)
        T = A[:, 0] / B
        Zn = A[:, 1:5] / B
        d, s = seg_dist(Zn)
        at45 = np.abs(s - 0.8) * B <= 3
        print(f"B={B} beta={B%5} n={len(T)} meanT/B={T.mean():.5f} sdT/B={T.std():.4f} "
              f"q10/50/90={np.quantile(T,0.1):.4f}/{np.quantile(T,0.5):.4f}/{np.quantile(T,0.9):.4f} maxT/B={T.max():.4f} "
              f"max seg-dist*B={d.max()*B:.1f} P(terminal at Q(4/5))={at45.mean():.4f} mean s={s.mean():.5f}", flush=True)


if __name__ == "__main__":
    main()
