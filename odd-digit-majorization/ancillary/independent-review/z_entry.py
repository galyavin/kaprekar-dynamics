"""Independent z3 encodings of Section 3 / Appendix A (real entry in nine digits).

L1  C^{-1} f C = H for 1>=a>=b>=c>=d>=0 (f by a 9-sorting network).
L1b first image: f(x) = C(sort|2x-1|) for 1>=x1>=..>=x4>=0.
L2  on the core: f(Cr) = R(Cr), and H(r) is again in the core (invariance).
L3  E consists of three-cycles (each phase maps to the next).
L4  the three exceptional points: fixed point and 2-cycle, each with a zero gap.
(i)-(iv) as stated in Appendix A, negated conclusion -> must be UNSAT.
"""
import sys, time, z3
from zmodel import *

QS = z3.Q


def base_r(s, name="r"):
    r = [z3.Real(f"{name}{i}") for i in range(4)]
    a, b, c, d = r
    s.add(1 > a, a >= b, b >= c, c >= d, d >= 0)
    return r


def iterates(r, n):
    out = [r]
    for _ in range(n):
        out.append(H_sym(out[-1]))
    return out


def report(name, s, expect=z3.unsat, t0=None):
    t = time.time()
    res = s.check()
    dt = time.time() - t
    ok = (res == expect)
    print(f"{name:55s} {str(res):8s} {'OK' if ok else '**FAIL**'}  {dt:7.1f}s", flush=True)
    if res == z3.sat:
        m = s.model()
        print("   model:", m)
    return ok


allok = True
# ---- L1b
s = z3.Solver()
x = [z3.Real(f"x{i}") for i in range(4)]
s.add(1 >= x[0], x[0] >= x[1], x[1] >= x[2], x[2] >= x[3], x[3] >= 0)
rad = sort_desc([zabs(2 * t - 1) for t in x])
fx = f_sym(x)
s.add(z3.Or([fx[i] != C_sym(rad)[i] for i in range(4)]))
allok &= report("L1b first image f(x)=C(sort|2x-1|)", s)

# ---- L1
s = z3.Solver()
r = [z3.Real(f"r{i}") for i in range(4)]
a, b, c, d = r
s.add(1 >= a, a >= b, b >= c, c >= d, d >= 0)
fy = f_sym(C_sym(r))
Hr = H_sym(r)
CH = C_sym(Hr)
s.add(z3.Or([fy[i] != CH[i] for i in range(4)]))
allok &= report("L1  f(C r) = C(H r)", s)

# ---- L2 (core: update is R and core invariant)
s = z3.Solver()
r = base_r(s)
s.add(in_core(r))
y = C_sym(r)
fy = f_sym(y)
Ry = R_sym(y)
s.add(z3.Or([fy[i] != Ry[i] for i in range(4)]))
allok &= report("L2a f(y)=R(y) on C(core)", s)
s = z3.Solver()
r = base_r(s)
s.add(in_core(r))
s.add(z3.Not(in_core(H_sym(r))))
allok &= report("L2b core invariant under H", s)
s = z3.Solver()
r = base_r(s)
s.add(in_core(r))
r3 = iterates(r, 3)[3]
s.add(z3.Or([r3[i] != r[i] for i in range(4)]))
allok &= report("L2c H^3 = id on core", s)
# is the core nonempty with interior, and does it avoid a=1?
s = z3.Solver()
r = [z3.Real(f"r{i}") for i in range(4)]
s.add(r[0] >= r[1], r[1] >= r[2], r[2] >= r[3], r[3] >= 0, r[0] <= 1)
s.add(in_core(r), r[0] >= 1)
report("L2d [info] core meets a=1 (boundary; excluded by a_4<1)", s, expect=z3.sat)

# ---- L3 E is a union of 3-cycles
s = z3.Solver()
A_, dd = z3.Reals("A dd")
s.add(A_ >= QS(5, 7), A_ <= QS(3, 4), dd > 0, dd <= 1 - A_)
ph1 = [A_, 2 * A_ - 1, 2 * A_ - 1, dd]
ph2 = [A_, 2 - 2 * A_ - dd, 3 * A_ - 2, 3 - 4 * A_]
ph3 = [A_, A_, 1 - A_ + dd, 1 - A_ - dd]
h1, h2, h3 = H_sym(ph1), H_sym(ph2), H_sym(ph3)
s.add(z3.Or(z3.Or([h1[i] != ph2[i] for i in range(4)]), z3.Or([h2[i] != ph3[i] for i in range(4)]),
            z3.Or([h3[i] != ph1[i] for i in range(4)])))
allok &= report("L3  E phases cycle 1->2->3->1", s)
# (phase predicates of App. A: dd=2-2a-b in (0,1-a] iff 1-a<=b<2-2a; dd'=1-a-d in (0,1-a] iff 0<=d<1-a: algebraic)

# ---- L4 exceptional points
from fractions import Fraction as Fr
import km
P0 = (Fr(1, 3),) * 4
P1 = (Fr(3, 5), Fr(3, 5), Fr(3, 5), Fr(1, 5))
P2 = (Fr(3, 5), Fr(1, 5), Fr(1, 5), Fr(1, 5))
print("L4  H(P0)==P0:", km.Hmap(P0) == P0, " H(P1)==P2:", km.Hmap(P1) == P2, " H(P2)==P1:", km.Hmap(P2) == P1,
      " zero gaps:", [min(km.Dvec(p)) for p in (P0, P1, P2)])
allok &= km.Hmap(P0) == P0 and km.Hmap(P1) == P2 and km.Hmap(P2) == P1

# ---- (i)
s = z3.Solver()
r = base_r(s)
it = iterates(r, 4)
for i in range(1, 5):
    s.add(it[i][0] == it[0][0])
s.add(it[4][0] < 1)
sD = [sort_desc(D_sym(t)) for t in it]
for i in range(1, 5):
    for j in range(5):
        s.add(sD[i][j] == sD[0][j])
r4 = it[4]
s.add(z3.Not(in_core(r4)), z3.Not(in_E(r4)))
for pt in (P0, P1, P2):
    s.add(z3.Not(eqv(r4, [z3.RealVal(str(t)) for t in pt])))
allok &= report("(i)  limit-cycle classification", s)

# ---- (ii) and (iii)
def hyp_ii(s, r, eta):
    a, b, c, d = r
    s.add(eta > 0, eta <= QS(1, 10000), a >= QS(5, 7) - eta, a <= QS(3, 4), b - c >= 0, b - c <= eta,
          zabs(c - 2 * a + 1) <= eta, d >= 100 * eta, d <= 1 - a + eta)


s = z3.Solver()
r = base_r(s)
eta = z3.Real("eta")
hyp_ii(s, r, eta)
it = iterates(r, 6)
s.add(it[6][0] != it[3][0])
allok &= report("(ii)  a_6 = a_3", s)

s = z3.Solver()
r = base_r(s)
eta = z3.Real("eta")
hyp_ii(s, r, eta)
s.add(r[0] >= QS(5, 7))
it = iterates(r, 6)
s.add(it[6][0] == r[0])
nu = lambda t: t[1] - t[2] + zabs(t[2] - 2 * t[0] + 1)
s.add(nu(it[6]) < nu(it[3]))
allok &= report("(iii) nu(r_6) >= nu(r_3)", s)

# ---- (iv)
s = z3.Solver()
r = base_r(s)
eta = z3.Real("eta")
s.add(eta > 0, W_sym(r) <= eta)
for t in D_sym(r):
    s.add(t >= 100 * eta)
it = iterates(r, 9)
s.add(z3.Not(in_core(it[9])))
s.add(z3.Or([it[9][i] - it[6][i] != it[6][i] - it[3][i] for i in range(4)]))
allok &= report("(iv)  r_9 - r_6 = r_6 - r_3 outside the core", s)

print("ALL ENTRY CHECKS OK" if allok else "SOME ENTRY CHECK FAILED")
