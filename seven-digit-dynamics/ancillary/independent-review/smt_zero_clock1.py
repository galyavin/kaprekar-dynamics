"""Lemma A.4 (first access, zero-origin returns) and Lemma A.5 one-step bound, own encoding.

Every 'claim' query asserts hypotheses + negated conclusion; UNSAT means the implication holds.
'sanity' queries are deliberately false strengthenings and are expected SAT; their models
are replayed with the plain Python model (model.F)."""
import sys, time
import z3
from smt_common import *
import model

results = []


def run(name, assertions, expect, replay=None, timeout_ms=600000):
    s = z3.Solver()
    s.set("timeout", timeout_ms)
    for a in assertions:
        s.add(a)
    t = time.time()
    r = s.check()
    dt = time.time() - t
    line = "%-48s %-7s expected %-6s %6.1fs" % (name, r, expect, dt)
    ok = (str(r) == expect)
    if r == z3.sat and replay is not None:
        m = s.model()
        rep = replay(m)
        line += "  replay: " + rep
    print(line + ("" if ok else "   <-- MISMATCH"))
    sys.stdout.flush()
    results.append((name, str(r), expect, ok))
    return r


validate_networks()
u, v, w, B, b = z3.Ints("u v w B b")
q = (u, v, w)

# ---------- Lemma A.4 (i): three positive states -> positive fourth (B >= 18) ----------
q1 = step_pos(q, B)
q2 = step_pos(q1, B)
q3 = step_pos(q2, B)
hyp_i = [valid(q, B), w >= 1, q1[2] >= 1, q2[2] >= 1]
run("A4(i) q,Fq,F2q positive & F3q zero, B>=18", hyp_i + [B >= 18, q3[2] == 0], "unsat")


def replay_i(m):
    BB, qq = m[B].as_long(), (m[u].as_long(), m[v].as_long(), m[w].as_long())
    x = [qq]
    for _ in range(3):
        x.append(model.F(x[-1], BB))
    return "B=%d orbit %s" % (BB, x)


run("sanity: A4(i) without B>=18 (B>=2)", hyp_i + [B >= 2, q3[2] == 0], "sat", replay_i)

# ---------- Lemma A.4 (ii): zero start (u,v,0), v>0, per parity, b>=9 ----------
for beta in (0, 1):
    z0 = (u, v, z3.IntVal(0))
    base = [B == 2 * b + beta, b >= 9, B > u, u >= v, v > 0]
    f1 = step_general(z0, B)
    f2 = step_general(f1, B)
    f3 = step_general(f2, B)
    if beta == 0:
        inlist = z3.And(u == b, v == b)
    else:
        inlist = z3.Or(z3.And(u == b, v == b), z3.And(u == b + 1, v == b), z3.And(u == b + 1, v == b + 1))
    run("A4(ii)a beta=%d: F(u,v,0) not positive" % beta, base + [f1[2] <= 0], "unsat")
    run("A4(ii)b1 beta=%d: F2 zero but source not listed" % beta, base + [f2[2] == 0, z3.Not(inlist)], "unsat")
    run("A4(ii)b2 beta=%d: listed source but F2 positive" % beta, base + [inlist, f2[2] != 0], "unsat")
    run("A4(ii)b3 beta=%d: F2 zero but F2 != (b+1,1,0)" % beta,
        base + [f2[2] == 0, z3.Or(f2[0] != b + 1, f2[1] != 1)], "unsat")
    run("A4(ii)c beta=%d: F2 positive but F3 not" % beta, base + [f2[2] > 0, f3[2] <= 0], "unsat")
    # the exceptional image (b+1,1,0) is not itself listed for b>=9
    run("A4(ii)x beta=%d: (b+1,1,0) in list" % beta,
        [B == 2 * b + beta, b >= 9, z3.substitute(inlist, (u, b + 1), (v, z3.IntVal(1)))], "unsat")

# ---------- Lemma A.4 (iii): first image of (u,0,0), all B>=2 ----------
zz = (u, z3.IntVal(0), z3.IntVal(0))
f1 = step_general(zz, B)
exp = (zmax(B - u, u - 1), zmin(B - u, u - 1), 0)
run("A4(iii) F(u,0,0) formula, B>=2", [B >= 2, B > u, u > 0, z3.Or(*[f1[i] != exp[i] for i in range(3)])], "unsat")

# ---------- Lemma A.5 one-step: V(Fq) <= V(q)+8 for positive q, all B>=2 ----------
for net6 in (True, False):
    f = step_pos(q, B, net6=net6)
    run("A5 one-step V(Fq)>V(q)+8 (%s)" % ("net6" if net6 else "oet7"),
        [B >= 2, valid(q, B), w >= 1, Vz(f, B) > Vz(q, B) + 8], "unsat")


def replay_V(m):
    BB, qq = m[B].as_long(), (m[u].as_long(), m[v].as_long(), m[w].as_long())
    return "B=%d q=%s V=%d V(Fq)=%d" % (BB, qq, model.V(qq, BB), model.V(model.F(qq, BB), BB))


f = step_pos(q, B)
run("sanity: V(Fq)>V(q)+7 (false strengthening)", [B >= 2, valid(q, B), w >= 1, Vz(f, B) > Vz(q, B) + 7], "sat", replay_V)
# cross-check the two networks agree on the positive step
f6 = step_pos(q, B, True)
f7 = step_pos(q, B, False)
run("net6 vs oet7 positive step disagree", [B >= 2, valid(q, B), w >= 1, z3.Or(*[f6[i] != f7[i] for i in range(3)])], "unsat")
# B-1 is the maximum digit for any valid positive state
run("B-1 not maximal digit", [valid(q, B), w >= 1, z3.Or(u > B - 1, v > B - 1, w - 1 > B - 1, B - w - 1 > B - 1, B - v - 1 > B - 1, B - u > B - 1)], "unsat")

print("ALL AS EXPECTED" if all(r[3] for r in results) else "SOME MISMATCH")
