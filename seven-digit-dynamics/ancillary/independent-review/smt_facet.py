"""Lemma A.3 first two assertions (facet return, one-step exit), own unbounded-integer encoding.

Uses the general step (case split on w>0, v>0) with the 7-input network so that positivity
of intermediate states is not assumed.  UNSAT = implication holds for all integers, no upper
bound on B.  Sanity queries (expected SAT) probe that the guards are not vacuous and that the
claimed update is not trivially forced."""
import sys, time
import z3
from smt_common import *
import model

u, v, w, B = z3.Ints("u v w B")
q0 = (u, v, w)
m = 2 * v - B
S = u + w - B
E = 2 * u - 2 * w - B
d = u - 3 * w
J = B - u - v + 2 * w
K1 = B - 2 * u + 2 * v - w
L = B - v - w

q1 = step_general(q0, B)
q2 = step_general(q1, B)
q3 = step_general(q2, B)
m3 = 2 * q3[1] - B

guards = [valid(q0, B), B >= 100000, u <= B - 1000, w >= 1000, m >= -6,
          S >= 16, S - m >= 16, E >= 16, E - m >= 16, J >= 8, K1 >= 8, L >= 8]
results = []


def run(name, A, expect):
    s = z3.Solver()
    s.set("timeout", 600000)
    s.add(A)
    t = time.time()
    r = s.check()
    line = "%-60s %-6s expected %-6s %.1fs" % (name, r, expect, time.time() - t)
    if r == z3.sat:
        mm = s.model()
        BB = mm[B].as_long()
        qq = (mm[u].as_long(), mm[v].as_long(), mm[w].as_long())
        x = qq
        for _ in range(3):
            x = model.F(x, BB)
        line += "  model B=%d q=%s -> F^3=%s" % (BB, qq, x)
    ok = str(r) == expect
    results.append(ok)
    print(line + ("" if ok else "  <-- MISMATCH"))
    sys.stdout.flush()


run("A3 return: guards & not (u-3, w-1, |m+1|-5)", guards + [z3.Or(q3[0] != u - 3, q3[2] != w - 1, m3 != zabs(m + 1) - 5)], "unsat")
run("sanity: guards satisfiable with m=-6", guards + [m == -6], "sat")
run("sanity: guards satisfiable with m>=6", guards + [m >= 6], "sat")
run("sanity: weaker J>=0 instead of J>=8 breaks the rule",
    [g for g in guards if not z3.eq(g, J >= 8)] + [J >= 0, z3.Or(q3[0] != u - 3, q3[2] != w - 1, m3 != zabs(m + 1) - 5)], "sat")

# exit step (Lemma A.3, second assertion)
exitres = []
for e in (12, 13, 14, 15):
    for mN in (-5, -4, -3, -2, -1):
        if (mN - e) % 2:
            continue
        Sp = (e + mN) // 2 - 1 if mN >= -1 else (e - mN) // 2 - 2
        assert (e + mN) % 2 == 0
        hyp = [valid(q0, B), B >= 100000, -d >= 8000, B + 6 * d >= 24000, E == e, m == mN]
        u1, v1, w1 = step_general(q0, B)
        mp = 2 * v1 - B
        Spp = u1 + w1 - B
        concl = z3.And(u1 == u - 1, mp == e + 2, Spp == Sp,
                       zabs(mp) <= 20, zabs(Spp) <= 20,
                       4 * u1 >= 3 * B + 4000, 6 * u1 <= 5 * B - 6000)
        run("A3 exit e=%d m_N=%d -> m'=%d S'=%d" % (e, mN, e + 2, Sp), hyp + [z3.Not(concl)], "unsat")
        run("  sanity: exit hypotheses satisfiable e=%d m=%d" % (e, mN), hyp, "sat")
print("ALL AS EXPECTED" if all(results) else "SOME MISMATCH")
