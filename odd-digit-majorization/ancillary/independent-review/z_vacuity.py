"""Non-vacuity and sensitivity checks for the independent z3 encodings.

(a) every hypothesis set used in an UNSAT query is satisfiable on its own;
(b) mutated statements (a constant tightened) become SAT, showing that the
    encodings can detect a violation and that the stated constants matter.
"""
import time, z3
from zmodel import *
from km import CORE

QS = z3.Q
B = z3.Int("B")


def rep(name, s, expect):
    t = time.time()
    res = s.check()
    ok = res == expect
    print(f"{name:70s} {str(res):6s} {'OK' if ok else '**UNEXPECTED**'} {time.time()-t:6.1f}s", flush=True)
    return ok


def base_r(s):
    r = [z3.Real(f"r{i}") for i in range(4)]
    s.add(1 > r[0], r[0] >= r[1], r[1] >= r[2], r[2] >= r[3], r[3] >= 0)
    return r


def its(r, n):
    out = [r]
    for _ in range(n):
        out.append(H_sym(out[-1]))
    return out


ok = True
# (i) hypotheses satisfiable, and with E's lower end 5/7 raised the conclusion fails
s = z3.Solver(); r = base_r(s); it = its(r, 4)
for i in range(1, 5):
    s.add(it[i][0] == it[0][0])
sD = [sort_desc(D_sym(t)) for t in it]
for i in range(1, 5):
    for j in range(5):
        s.add(sD[i][j] == sD[0][j])
s.add(it[4][0] < 1)
ok &= rep("(i) hypotheses satisfiable", s, z3.sat)
s.push(); s.add(z3.Not(in_core(it[4]))); ok &= rep("(i) hypotheses + outside core satisfiable (E points exist)", s, z3.sat); s.pop()


def in_E_mut(r, lo):
    a, b, c, d = r
    rng = z3.And(a >= lo, a <= QS(3, 4))
    p1 = z3.And(rng, b == 2 * a - 1, c == 2 * a - 1, d > 0, d <= 1 - a)
    p2 = z3.And(rng, c == 3 * a - 2, d == 3 - 4 * a, b >= 1 - a, b < 2 - 2 * a)
    p3 = z3.And(rng, b == a, c + d == 2 - 2 * a, d >= 0, d < 1 - a)
    return z3.Or(p1, p2, p3)


s.add(z3.Not(in_core(it[4])), z3.Not(in_E_mut(it[4], QS(72, 100))))
for pt in ((QS(1, 3),) * 4, (QS(3, 5), QS(3, 5), QS(3, 5), QS(1, 5)), (QS(3, 5), QS(1, 5), QS(1, 5), QS(1, 5))):
    s.add(z3.Not(z3.And([it[4][i] == pt[i] for i in range(4)])))
ok &= rep("(i) mutated: E range starts at 0.72 instead of 5/7 -> SAT", s, z3.sat)

# (ii) hypotheses satisfiable; with eta bound 1e-2 instead of 1e-4 and d>=eta, conclusion can fail?
def hyp(s, r, eta, etamax=QS(1, 10000), dmul=100):
    a, b, c, d = r
    s.add(eta > 0, eta <= etamax, a >= QS(5, 7) - eta, a <= QS(3, 4), b - c >= 0, b - c <= eta,
          zabs(c - 2 * a + 1) <= eta, d >= dmul * eta, d <= 1 - a + eta)


s = z3.Solver(); r = base_r(s); eta = z3.Real("eta"); hyp(s, r, eta)
ok &= rep("(ii) hypotheses satisfiable", s, z3.sat)
it = its(r, 6)
s.push(); s.add(it[6][0] == it[3][0], it[3][0] > r[0]); ok &= rep("(ii) hypotheses allow a_3 > a_0 (statement not trivial)", s, z3.sat); s.pop()
s = z3.Solver(); r = base_r(s); eta = z3.Real("eta"); hyp(s, r, eta, etamax=QS(1, 10), dmul=0); it = its(r, 6)
s.add(it[6][0] != it[3][0])
ok &= rep("(ii) mutated: eta<=1/10, d>=0 -> SAT", s, z3.sat)
# (iii) hypotheses satisfiable and nu strictly increases sometimes
s = z3.Solver(); r = base_r(s); eta = z3.Real("eta"); hyp(s, r, eta); s.add(r[0] >= QS(5, 7)); it = its(r, 6); s.add(it[6][0] == r[0])
nu = lambda t: t[1] - t[2] + zabs(t[2] - 2 * t[0] + 1)
ok &= rep("(iii) hypotheses satisfiable", s, z3.sat)
s.push(); s.add(nu(it[6]) > nu(it[3])); ok &= rep("(iii) nu can strictly increase (non-trivial)", s, z3.sat); s.pop()
# (iv) hypotheses satisfiable incl. r_9 outside core; mutated (no min D bound) SAT
s = z3.Solver(); r = base_r(s); eta = z3.Real("eta")
s.add(eta > 0, W_sym(r) <= eta)
for t in D_sym(r):
    s.add(t >= 100 * eta)
it = its(r, 9)
s.add(z3.Not(in_core(it[9])))
ok &= rep("(iv) hypotheses satisfiable (with r_9 outside core)", s, z3.sat)
s = z3.Solver(); r = base_r(s); eta = z3.Real("eta")
s.add(eta > 0, W_sym(r) <= eta)
for t in D_sym(r):
    s.add(t >= 2 * eta)
it = its(r, 9)
s.add(z3.Not(in_core(it[9])))
s.add(z3.Or([it[9][i] - it[6][i] != it[6][i] - it[3][i] for i in range(4)]))
ok &= rep("(iv) mutated: min D >= 2 eta instead of 100 eta -> SAT", s, z3.sat)

# ---------------- integer side
def qv(t):
    return [z3.Int(f"{t}{i}") for i in range(4)]


for i in range(1, 5):
    s = z3.Solver(); q = qv("q"); s.add(B >= 2, valid_q(q, B), q[3] >= 1)
    s.add(S_sym(gaps_sym(Fsym(q, B), B), i) > S_sym(gaps_sym(q, B), i))
    ok &= rep(f"Z2 mutated: S_{i} with +0 -> SAT (the +1 is attained)", s, z3.sat)

s = z3.Solver(); q0 = qv("q"); s.add(B >= 2, valid_q(q0, B), q0[3] >= 1)
qs = [q0]
for t in range(4):
    qs.append(Fsym(qs[-1], B))
for t in range(1, 4):
    s.add(qs[t][3] >= 1)
ok &= rep("Z3 hypotheses (four positive states) satisfiable", s, z3.sat)
s.add(qs[4][3] == 0)
ok &= rep("Z3 mutated: B>=2 instead of B>=40 -> SAT (bound matters)", s, z3.sat)
if s.check() == z3.sat:
    m = s.model(); print("     witness B =", m[B], [m.eval(t) for t in q0])

# cover: removing the |4u-3B|<=8 row at rank 3 must give SAT for some j
def ab(x):
    return z3.If(x >= 0, x, -x)


def cover3_minus5(q, B):
    u, v, w, z = q
    row1 = z3.And(v == w, z3.Or(ab(u - v) <= 2, ab(u + v - B) <= 2, ab(2 * v - B) <= 3))
    row2 = z3.And(u == v, z3.Or(ab(u + w - B) <= 2, ab(2 * u - B) <= 3))
    row3 = z3.And(ab(2 * u - B) <= 3, ab(2 * v - B) <= 3)
    row4 = z3.And(B - u <= 3, z3.Or(ab(2 * v - B) <= 3, v == w, ab(v + w - B) <= 3))
    return z3.Or(row1, row2, row3, row4)


found = False
for j in (1, 2, 3, 4):
    s = z3.Solver(); q0 = qv("q")
    s.add(B >= 40, valid_q(q0, B), q0[3] == 0, q0[2] > 0, z3.Not(cover3_minus5(q0, B)))
    qs = [q0]
    for t in range(j):
        qs.append(Fsym(qs[-1], B))
    for t in range(1, j):
        s.add(qs[t][3] >= 1)
    s.add(qs[j][3] == 0)
    if s.check() == z3.sat:
        m = s.model(); print(f"   cover row 5 needed: j={j} witness B={m[B]} q0={[m.eval(t) for t in q0]}")
        found = True
        break
print(f"{'Z4 mutated: rank-3 cover without row 5 -> SAT for some j':70s} {'OK' if found else '**UNEXPECTED**'}")
ok &= found
# returning zero states exist at rank 2 and 3 for every j? (non-vacuity of each query)
for rank in (2, 3):
    for j in (1, 2, 3, 4):
        s = z3.Solver(); q0 = qv("q")
        s.add(B >= 40, valid_q(q0, B), q0[3] == 0, q0[0] > 0)
        s.add(q0[1] > 0, q0[2] == 0) if rank == 2 else s.add(q0[2] > 0)
        qs = [q0]
        for t in range(j):
            qs.append(Fsym(qs[-1], B))
        for t in range(1, j):
            s.add(qs[t][3] >= 1)
        s.add(qs[j][3] == 0)
        res = s.check()
        print(f"   Z4 hypotheses rank {rank}, j={j}: {res}")

# bulk translation: margin 0 -> SAT (margin matters) ; margin satisfiable
s = z3.Solver(); q = qv("q"); s.add(B >= 2, valid_q(q, B))


def core_int(q, B):
    u, v, w, z = q
    a = 2 * u - B; b = 2 * v - 2 * u + B; c = 2 * w - 2 * v + 2 * u - B; d = 2 * z - 2 * w + 2 * v - 2 * u + B
    return [row[0] * a + row[1] * b + row[2] * c + row[3] * d + row[4] * B for row in CORE]


for t in core_int(q, B):
    s.add(t >= 16)
ok &= rep("Z5 hypotheses (slacks>=16) satisfiable", s, z3.sat)
s = z3.Solver(); q = qv("q"); s.add(B >= 2, valid_q(q, B))
for t in core_int(q, B):
    s.add(t >= 0)
q3 = q
for _ in range(3):
    q3 = Fsym(q3, B)
s.add(z3.Or([q3[i] != q[i] + (-3, -2, -1, -1)[i] for i in range(4)]))
ok &= rep("Z5 mutated: slack margin 0 -> SAT (faces stop the bulk drift)", s, z3.sat)
print("VACUITY/SENSITIVITY OK" if ok else "SOME UNEXPECTED")
