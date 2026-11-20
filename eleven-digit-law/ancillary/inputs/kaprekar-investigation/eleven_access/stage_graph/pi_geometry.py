"""Exact geometry, first-face weights and entrance remainders."""
from collections import Counter
from fractions import Fraction as Q
import json
from math import factorial, gcd, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sb_exhaustion"))
from geometry import Poly, det, dot, factorial_product  # noqa: E402


def minus(a, b):
    return tuple(u-v for u, v in zip(a, b))


def within_interlevel():
    base = (1, 1, 1, 1, 1, 1, 0)
    names = list(range(6))
    offsets = [0]*6
    by_row = {}
    up = set()
    for phase in range(5):
        q = [tuple(int(j in names[i+1:]) for j in range(6))
             + (sum(offsets[i+1:]),) for i in range(5)]
        a = [row[:6]+(row[6]-int(i == 4),) for i, row in enumerate(q)]
        comp = [minus(base, q[4-i]) for i in range(5)]
        b = [row[:6]+(row[6]-int(i < 4),) for i, row in enumerate(comp)]
        for i in range(5):
            row = minus(a[i], b[i])
            if row[:-1] not in by_row or row[-1] < by_row[row[:-1]][-1]:
                by_row[row[:-1]] = (phase, i+1, row[-1])
        for i in range(4):
            up.add(minus(b[i], a[i+1])[:-1])
        perm = (0, 5, 1, 4, 2, 3)
        u = (1, -2, 1, 0, 0, 0)
        names = [names[i] for i in perm]
        offsets = [offsets[i]+u[j] for j, i in enumerate(perm)]
    assert names == list(range(6)) and offsets == [5,-1,-1,-1,-1,-1]
    assert set(dot(r, offsets) for r in by_row) == {-6}
    assert set(dot(r, offsets) for r in up) == {5}
    return by_row, up


def dirichlet_mass(poly):
    if poly is None:
        return Q(0)
    total = Q(0)
    for sx in poly.simplices():
        vertices = [poly.vertices[i] for i in sx]
        den = lcm(*(x.denominator for v in vertices for x in v))
        terms = {(0,)*6: 1}
        for gap in range(6):
            extra = Counter()
            for powers, coeff in terms.items():
                for i, vert in enumerate(vertices):
                    p = powers[:i]+(powers[i]+1,)+powers[i+1:]
                    if vert[gap]:
                        extra[p] += coeff*int(den*vert[gap])
            terms = extra
        # The normalization 11! cancels the total degree-6 simplex factor.
        total += abs(det(vertices))*Q(sum(c*factorial_product(p) for p,c in terms.items()), den**6)
    return total


def relative_volume(poly):
    if poly is None:
        return Q(0)
    return sum(abs(det([poly.vertices[i] for i in sx])) for sx in poly.simplices())


def run():
    desc, inter = within_interlevel()
    poly = Poly.root().restrict([*sorted(inter), *sorted(desc)])
    full_vol = relative_volume(poly)
    full_mass = dirichlet_mass(poly)
    assert dirichlet_mass(Poly.root()) == 1
    rows = []
    for a, (phase, pair, c) in sorted(desc.items()):
        winner = poly.restrict(minus(b,a) for b in sorted(desc) if b != a)
        volume = relative_volume(winner)
        gap_mass = dirichlet_mass(winner)
        stride = gcd(6, *(a[i]-a[5] for i in range(5)))
        rows.append({
            "phase": phase, "pair": pair, "gap_row": a, "offset": c,
            "relative_simplex_volume": volume, "uniform_pi_mass": volume/full_vol,
            "original_direct_pi_mass": gap_mass,
            "by_base_mod_12": {
                r: {"remainders": list(range((a[5]*r+c)%stride, 6, stride)),
                    "each": Q(stride,6)}
                for r in range(12)}})
    assert sum(r["relative_simplex_volume"] for r in rows) == full_vol
    assert sum(r["original_direct_pi_mass"] for r in rows) == full_mass
    return {"pi_relative_simplex_volume": full_vol,
            "pi_q_volume": full_vol/factorial(5),
            "pi_original_direct_weight": full_mass,
            "distinct_exit_normals": len(desc),
            "nonempty_strict_cells": sum(r["relative_simplex_volume"] > 0 for r in rows),
            "rows": rows}


if __name__ == "__main__":
    print(json.dumps(run(), default=str, indent=2))

