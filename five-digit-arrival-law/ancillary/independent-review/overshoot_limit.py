"""Limit predictions of Theorem 5, equation (12), by numerical quadrature.
usage: python overshoot_limit.py alpha beta gamma zeta"""
import sys
import numpy as np
al, be, ga, ze = map(float, sys.argv[1:5])
n = 2000
lam = al + (np.arange(n) + 0.5) * (be - al) / n
z = ga + (np.arange(n) + 0.5) * (ze - ga) / n
L, Z = np.meshgrid(lam, z, indexing="ij")
dens = 120 * L * Z * (1 - Z)
cell = (be - al) * (ze - ga) / n**2
M = dens.sum() * cell
Mexact = (be**2 - al**2) / 2 * 120 * ((ze**2 - ga**2) / 2 - (ze**3 - ga**3) / 3)
S = np.sqrt(L**2 + 2 * Z - 1)
EY = (dens * (S / 6 - L)).sum() * cell / M
print(f"M = {Mexact:.6f} (quadrature {M:.6f}); limit B*P(E_B) -> M")
print(f"limit E[Y|E] = {EY:.5f}")
for y in (-1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25):
    p = np.clip((y + L) / (S / 3), 0, 1)
    print(f"limit P(Y<={y:.2f}|E) = {(dens * p).sum() * cell / M:.5f}")
