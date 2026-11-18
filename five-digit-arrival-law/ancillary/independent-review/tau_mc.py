"""Plain Monte Carlo of the paper's limiting variable tau(x,y) of (13)-(15),
sampling (x,y) as outer/next-outer gaps of five sorted uniforms (density (1)).
Cross-checks limit_exact.py.  usage: python tau_mc.py nsamples"""
import sys
import numpy as np
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
rng = np.random.default_rng(424242)
ts = [1/25, 0.1, 0.2, 0.3]
hits = np.zeros(4); s1 = 0.0; done = 0
while done < n:
    m = min(2_000_000, n - done)
    p = np.sort(rng.random((m, 5)), axis=1)
    x = p[:, 4] - p[:, 0]; y = p[:, 3] - p[:, 1]
    U = np.where(x < 0.5, 1 - y, np.where(x + y < 1, 1 - y, x))
    V = np.where(x < 0.5, 1 - x - y, np.where(x + y < 1, x - y, np.where(y < 0.5, x - y, x + y - 1)))
    d = 1 - U
    nn = np.maximum(0, np.ceil((V - 0.5) / d))
    w = V - nn * d
    A = U + w - 1; C = 2 * U - w - 1
    z = np.zeros_like(A)
    tau = (2 / 3) * (np.maximum(np.maximum(z, A), C) - np.minimum(np.minimum(z, A), C))
    for k, t in enumerate(ts):
        hits[k] += np.count_nonzero(tau <= t)
    s1 += tau.sum(); done += m
for k, t in enumerate(ts):
    p = hits[k] / n
    print(f"MC F({t}) = {p:.6f} +- {np.sqrt(p*(1-p)/n):.6f}")
print(f"MC E[L] = {s1/n:.6f}")
