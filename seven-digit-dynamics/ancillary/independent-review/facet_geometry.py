"""Lemma A.3 remaining assertions + Corollary 6, own derivation.

(1) vertices of the canonical sign regions, computed from two descriptions:
    (a) C(K) in (x,y,z) with the canonical rotation inequalities and the sign of x-3z;
    (b) homogeneous limits of the source inequalities (eq. source) + ordered coordinates + sign;
(2) the residue (fast-column) tables from my own edge graph (edge.py) and the facet update;
(3) M_Q and M_P by exact integration of rho over my own polytope decomposition in xi-space
    (vertex enumeration -> facet fans -> tetrahedra -> sympy exact integration), no use of the
    authors' tetrahedra or barycentric formula; Monte Carlo cross-check.
"""
from fractions import Fraction as Fr
from itertools import combinations
import math, random, sys
import sympy as sp

from edge import build_graph, analyse

x, y, z = sp.symbols("x y z")


def lin(expr):
    """affine expr in x,y,z -> (cx, cy, cz, c0) Fractions"""
    e = sp.expand(expr)
    return tuple(Fr(str(e.coeff(s))) for s in (x, y, z)) + (Fr(str(e.subs({x: 0, y: 0, z: 0}))),)


def solve3(rows):
    # rows: (cx,cy,cz,c0) meaning cx x + cy y + cz z + c0 = 0
    M = [[r[0], r[1], r[2], -r[3]] for r in rows]
    n = 3
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [M[r][k] - f * M[col][k] for k in range(4)]
    return tuple(M[i][3] / M[i][i] for i in range(3))


def ev(g, p):
    return g[0] * p[0] + g[1] * p[1] + g[2] * p[2] + g[3]


def vertices(ineqs):
    V = set()
    for tri in combinations(ineqs, 3):
        p = solve3(tri)
        if p is None:
            continue
        if all(ev(g, p) >= 0 for g in ineqs):
            V.add(p)
    return V


# ---------------- (1) canonical regions -----------------
a = 2 * x - 1
b = 2 * y - a
c = 2 * z - b
ordered_abc = [1 - a, a - b, b - c, c]
Kineq = [2 * a - 1, a + b - 1, a - b - c, b + 2 * c - 2 * a + 1, 2 - a - 2 * b - c]
canon = [x + z - 2 * y, x - y - z]
d = x - 3 * z
desc_a = {s: [lin(e) for e in ordered_abc + Kineq + canon + [s * d]] for s in (1, -1)}
# homogeneous source: u<=B, w>=0, m>=0, S-m>=0, E-m>=0, J,K1,L>=0 with B=1, plus ordered coords
src = [1 - x, z, 2 * y - 1, (x + z - 1) - (2 * y - 1), (2 * x - 2 * z - 1) - (2 * y - 1),
       1 - x - y + 2 * z, 1 - 2 * x + 2 * y - z, 1 - y - z, x - y, y - z]
desc_b = {s: [lin(e) for e in src + [s * d]] for s in (1, -1)}
named = {"A": (Fr(3, 4), Fr(1, 2), Fr(1, 4)), "U": (Fr(5, 6), Fr(1, 2), Fr(1, 6)), "G": (Fr(6, 7), Fr(1, 2), Fr(2, 7)),
         "D": (Fr(9, 10), Fr(1, 2), Fr(1, 5)), "F": (Fr(5, 6), Fr(1, 2), Fr(1, 3)), "V": (Fr(1), Fr(2, 3), Fr(1, 3))}
inv = {v: k for k, v in named.items()}
geo_ok = True
for s in (1, -1):
    Va, Vb = vertices(desc_a[s]), vertices(desc_b[s])
    la = sorted(inv.get(p, str(p)) for p in Va)
    lb = sorted(inv.get(p, str(p)) for p in Vb)
    claim = sorted(["A", "U", "G", "D", "V"] if s == 1 else ["A", "G", "F", "V"])
    print("sign %+d: vertices from C(K) description %s ; from source description %s ; paper %s ; match=%s"
          % (s, la, lb, claim, la == lb == claim))
    geo_ok &= (la == lb == claim)
    dv = [p[0] - 3 * p[2] for p in Va]
    print("   d=x-3z over vertices: min %s max %s" % (min(dv), max(dv)))
    tf = [(9 * p[0] - 3 * p[2] - 6) / 8 if s == 1 else 3 * (p[0] + p[2] - 1) / 4 for p in Va]
    ts = [(5 * p[0] - 3 * p[2] - 3) / 4 if s == 1 else (2 * p[0] + 6 * p[2] - 3) / 4 for p in Va]
    print("   t_f over vertices %s ; t_s over vertices %s" % (sorted(set(tf)), sorted(set(ts))))

# ---------------- (2) residue tables -----------------
rate = {}
for beta in (0, 1):
    start, rule, fails = build_graph(beta)
    cycles, *_ = analyse(beta, start, rule)
    cyc_of = {}
    for key, (cyc, ell, D) in cycles.items():
        for p in cyc:
            cyc_of[p] = Fr(3 * ell, D)

    def r_of(node, rule=rule, cyc_of=cyc_of):
        seen = set()
        while node not in cyc_of:
            seen.add(node)
            node = rule[node][0]
        return cyc_of[node]

    rate[beta] = r_of


def mstep(mm):
    return abs(mm + 1) - 5


paper = {  # row -> (pos even, pos odd, neg even, neg odd)
    0: ("12", "02", "0123", "."), 1: ("03", "02", ".", "-"), 2: ("03", "-", "03", "."), 3: ("03", "02", ".", "023"),
    4: ("03", "02", "0123", "."), 5: ("12", "02", ".", "-"), 6: ("12", "-", "12", "."), 7: ("12", "02", ".", "012")}
mine = {}
counts = {}
for side in (1, -1):
    for beta in (0, 1):
        cols = [r for r in range(8) if r % 2 == beta]
        nfast = 0
        ntot = 0
        for row in range(8):
            if side == -1 and row % 2 != beta:
                mine[(row, side, beta)] = "."
                continue
            fastcols = ""
            for ci, mres in enumerate(cols):
                outcomes = set()
                for SR in (2000, 2408, 4000, 123456 * 8):  # several representatives
                    Sv = SR + ((row - SR) % 8)
                    for mR in (0, 8, 40, 800):
                        mv = mR + mres
                        if side == 1:
                            N = (Sv - 16) // 4
                            s_ = Sv - 4 * N
                        else:
                            N = (Sv - 12) // 4
                            s_ = Sv - 4 * N
                        mm = mv
                        for _ in range(N):
                            mm = mstep(mm)
                        if side == 1:
                            node = (mm, s_)
                        else:
                            e = s_
                            Sp = (e + mm) // 2 - 1 if mm >= -1 else (e - mm) // 2 - 2
                            node = (e + 2, Sp)
                        outcomes.add(rate[beta](node))
                assert len(outcomes) == 1, (side, beta, row, mres, outcomes)
                r = outcomes.pop()
                ntot += 1
                if r == Fr(3, 2):
                    fastcols += str(ci)
                    nfast += 1
                else:
                    assert r == 2
            mine[(row, side, beta)] = fastcols if fastcols else "-"
        counts[(side, beta)] = (nfast, ntot)
tab_ok = True
print("row | d>0 even | d>0 odd | d<0 even | d<0 odd   (mine vs paper)")
for row in range(8):
    got = (mine[(row, 1, 0)], mine[(row, 1, 1)], mine[(row, -1, 0)], mine[(row, -1, 1)])
    print("  %d | %s  %s" % (row, " ".join("%5s" % g for g in got), "OK" if got == paper[row] else "DIFF %s" % (paper[row],)))
    tab_ok &= got == paper[row]
for k, (nf, nt) in sorted(counts.items()):
    print("  side %+d beta=%d: fast %d / %d = %s" % (k[0], k[1], nf, nt, Fr(nf, nt)))

# ---------------- (3) masses -----------------
rho = 5040 * (1 - x) * (x - y) * (y - z) * z
s1, s2, s3 = sp.symbols("s1 s2 s3")
std_cache = {}


def integrate_tet(P0, P1, P2, P3, f=rho):
    M = sp.Matrix([[sp.Rational(str(P1[i] - P0[i])), sp.Rational(str(P2[i] - P0[i])), sp.Rational(str(P3[i] - P0[i]))] for i in range(3)])
    det = abs(M.det())
    if det == 0:
        return sp.Integer(0)
    sub = {x: sp.Rational(str(P0[0])) + M[0, 0] * s1 + M[0, 1] * s2 + M[0, 2] * s3,
           y: sp.Rational(str(P0[1])) + M[1, 0] * s1 + M[1, 1] * s2 + M[1, 2] * s3,
           z: sp.Rational(str(P0[2])) + M[2, 0] * s1 + M[2, 1] * s2 + M[2, 2] * s3}
    g = sp.Poly(sp.expand(f.subs(sub)), s1, s2, s3)
    tot = sp.Integer(0)
    for (k1, k2, k3), coef in g.terms():
        # integral of s1^k1 s2^k2 s3^k3 over standard simplex = k1!k2!k3!/(k1+k2+k3+3)!
        tot += coef * sp.factorial(k1) * sp.factorial(k2) * sp.factorial(k3) / sp.factorial(k1 + k2 + k3 + 3)
    return tot * det


def polytope_integral(ineqs, f=rho):
    V = list(vertices(ineqs))
    if len(V) < 4:
        return sp.Integer(0), 0
    # check full-dimensional
    v0 = V[0]
    tets = []
    seen_faces = set()
    for g in ineqs:
        on = [p for p in V if ev(g, p) == 0]
        if len(on) < 3 or v0 in on:
            continue
        fk = frozenset(on)
        if fk in seen_faces:  # the same facet can be cut out by several (duplicate) inequalities
            continue
        seen_faces.add(fk)
        # order cyclically: project onto plane using float coordinates
        cen = [sum(float(p[i]) for p in on) / len(on) for i in range(3)]
        nrm = [float(g[0]), float(g[1]), float(g[2])]
        # basis in plane
        t = [1.0, 0.0, 0.0] if abs(nrm[0]) < 0.9 * max(abs(n) for n in nrm) else [0.0, 1.0, 0.0]
        e1 = [t[1] * nrm[2] - t[2] * nrm[1], t[2] * nrm[0] - t[0] * nrm[2], t[0] * nrm[1] - t[1] * nrm[0]]
        e2 = [nrm[1] * e1[2] - nrm[2] * e1[1], nrm[2] * e1[0] - nrm[0] * e1[2], nrm[0] * e1[1] - nrm[1] * e1[0]]
        ang = lambda p: math.atan2(sum((float(p[i]) - cen[i]) * e2[i] for i in range(3)), sum((float(p[i]) - cen[i]) * e1[i] for i in range(3)))
        on = sorted(set(on), key=ang)
        for i in range(1, len(on) - 1):
            tets.append((v0, on[0], on[i], on[i + 1]))
    tot = sum((integrate_tet(*t, f=f) for t in tets), sp.Integer(0))
    return sp.nsimplify(tot), len(tets)


# sanity: full ordered simplex D: int rho = 1, volume 1/6
Dineq = [lin(e) for e in (1 - x, x - y, y - z, z)]
I1, _ = polytope_integral(Dineq)
Ivol, _ = polytope_integral(Dineq, f=sp.Integer(1))
print("sanity: int_D rho =", I1, " vol(D) =", Ivol)

# M_Q: Q = C(K cap ordered simplex), described in xi coordinates
Qineq = [lin(e) for e in ordered_abc + Kineq]
MQ, nQ = polytope_integral(Qineq)
# M_P: P = {xi in D : sort(|2 xi - 1|) in K}; split D by signs of 2xi-1 and the radius order
MP = sp.Integer(0)
tot_pieces = sp.Integer(0)
from itertools import permutations, product
npieces = 0
for signs in product((1, -1), repeat=3):
    rad = [sg * (2 * t - 1) for sg, t in zip(signs, (x, y, z))]
    for perm in permutations(range(3)):
        A_, B_, C_ = rad[perm[0]], rad[perm[1]], rad[perm[2]]
        piece = [1 - x, x - y, y - z, z] + rad + [A_ - B_, B_ - C_]
        pin = [lin(e) for e in piece]
        full = [lin(e) for e in piece + [2 * A_ - 1, A_ + B_ - 1, A_ - B_ - C_, B_ + 2 * C_ - 2 * A_ + 1, 2 - A_ - 2 * B_ - C_]]
        Ipc, _ = polytope_integral(pin)
        tot_pieces += Ipc
        Ik, nt = polytope_integral(full)
        if nt:
            npieces += 1
        MP += Ik
print("sanity: sum over the 48 sign/order pieces of D of int rho =", sp.nsimplify(tot_pieces))
print("M_Q = %s  (paper 1339759/43200000)  match=%s   [%d tetrahedra]" % (MQ, MQ == sp.Rational(1339759, 43200000), nQ))
print("M_P = %s  (paper 937297/5400000)   match=%s   [%d nonempty pieces]" % (MP, sp.nsimplify(MP) == sp.Rational(937297, 5400000), npieces))

# Monte Carlo cross-check
rnd = random.Random(7)
NMC = 2000000
hq = hp = 0
sw = 0.0
for _ in range(NMC):
    xi = sorted((rnd.random(), rnd.random(), rnd.random()), reverse=True)
    X, Y, Z = xi
    wgt = 5040 * (1 - X) * (X - Y) * (Y - Z) * Z / 6  # density w.r.t. uniform on D (vol 1/6)
    aa = 2 * X - 1
    bb = 2 * Y - aa
    cc = 2 * Z - bb

    def inK(a_, b_, c_):
        return (1 >= a_ >= b_ >= c_ >= 0 and 2 * a_ >= 1 and a_ + b_ >= 1 and a_ - b_ - c_ >= 0
                and b_ + 2 * c_ - 2 * a_ + 1 >= 0 and 2 - a_ - 2 * b_ - c_ >= 0)
    if inK(aa, bb, cc):
        hq += wgt
    r = sorted((abs(2 * X - 1), abs(2 * Y - 1), abs(2 * Z - 1)), reverse=True)
    if inK(*r):
        hp += wgt
print("Monte Carlo (%d pts): M_Q ~ %.5f (exact %.5f), M_P ~ %.5f (exact %.5f)" % (NMC, hq / NMC, float(MQ), hp / NMC, float(MP)))
print("GEOMETRY OK" if geo_ok and tab_ok else "GEOMETRY/TABLE MISMATCH")
