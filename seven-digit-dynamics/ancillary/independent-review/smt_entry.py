"""Lemma A.1 (eight real inequalities), own QF_LRA encoding.

H(a,b,c) = sort_desc(a, |a+b-1|, |b+c-1|), with max/min by nested If and the median as
max(min(x,y), min(max(x,y), z)) (the authors use sum-minus-extremes).
Each claim query is hypotheses + negated conclusion (UNSAT = implication holds).
Sanity queries are deliberately false variants (expected SAT); their models are replayed
in exact rational arithmetic.
Writes SMT-LIB2 copies to smt2/entry_*.smt2 for an independent Yices run."""
import os, sys, time
from fractions import Fraction as Fr
import z3

os.makedirs("smt2", exist_ok=True)


def zmax(x, y):
    return z3.If(x >= y, x, y)


def zmin(x, y):
    return z3.If(x >= y, y, x)


def zabs(x):
    return z3.If(x >= 0, x, -x)


def pos(x):
    return z3.If(x >= 0, x, 0)


def H(q):
    a, b, c = q
    x, y, z = a, zabs(a + b - 1), zabs(b + c - 1)
    hi = zmax(zmax(x, y), z)
    lo = zmin(zmin(x, y), z)
    md = zmax(zmin(x, y), zmin(zmax(x, y), z))
    return (hi, md, lo)


def W(q):
    a, b, c = q
    return pos(1 - a - b) + pos(b + c - a) + pos(2 * a - 1 - b - 2 * c) + pos(2 * b + c + a - 2)


# exact replay
def Hf(q):
    a, b, c = q
    return tuple(sorted((a, abs(a + b - 1), abs(b + c - 1)), reverse=True))


def Wf(q):
    a, b, c = q
    p = lambda t: t if t > 0 else 0
    return p(1 - a - b) + p(b + c - a) + p(2 * a - 1 - b - 2 * c) + p(2 * b + c + a - 2)


a, b, c = z3.Reals("a b c")
q = [(a, b, c)]
for j in range(6):
    q.append(H(q[-1]))
A = [qq[0] for qq in q]
dom = [0 <= c, c <= b, b <= a, a < 1]
half, f35 = z3.RealVal("1/2"), z3.RealVal("3/5")
s = a + b + c - 1
inZ = z3.And(half <= a, a <= f35, z3.Or(z3.And(b == 2 - 3 * a, c == 2 * a - 1), z3.And(b == a, c == 2 * a - 1)))

claims = {
    "i-1": [a < half, A[1] < half, 2 * A[1] + q[1][1] < 1],
    "i-2": [a < half, 2 * a + b >= 1,
            z3.Or(q[1][0] != a + pos(-s), q[1][1] != a - pos(s), q[1][2] != 1 - a - b)],
    "ii": [a >= half, W(q[3]) > W(q[0])],
    "iii-1": [a > f35, A[3] != A[2]],
    "iii-2": [a > f35, A[6] == a, W(q[6]) > zmax(0, W(q[0]) - zmin(5 * a - 3, 1 - a))],
    "iv": [a > half, A[6] < 1, W(q[6]) == W(q[0]), W(q[0]) > 0, z3.Not(inZ)],
    "v": [a > half, W(q[0]) < zmin(2 * a - 1, 1 - a) / 2, W(q[6]) != 0],
    "vi": [half < a, a <= f35, z3.Not(z3.And(b == 2 - 3 * a, c == 2 * a - 1)),
           zmax(zabs(b - 2 + 3 * a), zabs(c - 2 * a + 1)) < (2 * a - 1) / 100]
          + [z3.And(A[j] <= f35, W(q[j]) >= 2 * A[j] - 1) for j in range(5)],
}
sanity = {
    # strict three-step decrease is false
    "S-ii-strict": [a >= half, W(q[3]) >= W(q[0]), W(q[0]) > 0],
    # (iv) fails without the segment exclusion (Z points exist)
    "S-iv-Z": [a > half, A[6] < 1, W(q[6]) == W(q[0]), W(q[0]) > 0],
    # (vi) with only j<=1 allowed should fail
    "S-vi-short": [half < a, a <= f35, z3.Not(z3.And(b == 2 - 3 * a, c == 2 * a - 1)),
                   zmax(zabs(b - 2 + 3 * a), zabs(c - 2 * a + 1)) < (2 * a - 1) / 100]
                  + [z3.And(A[j] <= f35, W(q[j]) >= 2 * A[j] - 1) for j in range(2)],
    # (v) with a larger threshold (factor 2 instead of 1/2) should fail
    "S-v-wide": [a > half, W(q[0]) < 2 * zmin(2 * a - 1, 1 - a), W(q[6]) != 0],
}


def frac(v):
    v = v.as_fraction() if hasattr(v, "as_fraction") else None
    return v


ok = True
for name, hyp in list(claims.items()) + list(sanity.items()):
    sol = z3.Solver()
    sol.add(dom + hyp)
    with open("smt2/entry_%s.smt2" % name, "w") as fh:
        fh.write("(set-logic QF_LRA)\n" + sol.to_smt2())
    t = time.time()
    r = sol.check()
    exp = "sat" if name.startswith("S-") else "unsat"
    line = "A1 %-12s %-6s expected %-6s %.2fs" % (name, r, exp, time.time() - t)
    if r == z3.sat:
        m = sol.model()
        pt = tuple(m.eval(x, model_completion=True).as_fraction() for x in (a, b, c))
        orb = [pt]
        for _ in range(6):
            orb.append(Hf(orb[-1]))
        line += "  model (a,b,c)=(%s) W0=%s W3=%s W6=%s a6=%s" % (
            ",".join(str(x) for x in pt), Wf(orb[0]), Wf(orb[3]), Wf(orb[6]), orb[6][0])
    if str(r) != exp:
        ok = False
        line += "   <-- MISMATCH"
    print(line)
    sys.stdout.flush()
print("ALL AS EXPECTED" if ok else "SOME MISMATCH")
