"""Lemma A.5 depth-nine claim, own encoding.

Claim: for every integer B>=2 and positive q with q,F q,...,F^9 q all positive (w>=1) and
V(q)>200, V(F^9 q) < V(q).  Negation: ... and V_9 >= V_0.

Usage: smt_clock9.py MODE [args]
  MODE=plain               single query, no auxiliary constraints
  MODE=aux                 add V_j <= V_{j-1}+8 (one-step lemma, proved separately) and
                           V_j >= V_0 - 8(9-j) (consequence of V_9>=V_0 and the one-step lemma)
  MODE=part K              aux + source sign pattern K in 0..15 (4g_i >= B or < B)
  MODE=yices K             emit SMT-LIB2 for pattern K (for the Yices binary)
"""
import sys, time
import z3
from smt_common import *

DEPTH = 9


def build(mode, k=None):
    u, v, w, B = z3.Ints("u v w B")
    qs = [(u, v, w)]
    for j in range(DEPTH):
        qs.append(step_pos(qs[-1], B))
    Vs = [Vz(qq, B) for qq in qs]
    A = [B >= 2, valid(qs[0], B)] + [qq[2] >= 1 for qq in qs]
    A += [Vs[0] > 200, Vs[DEPTH] >= Vs[0]]
    if mode in ("aux", "part", "yices", "yhint"):
        Vv = z3.Ints(" ".join("V%d" % j for j in range(DEPTH + 1)))
        A += [Vv[j] == Vs[j] for j in range(DEPTH + 1)]
        for j in range(1, DEPTH + 1):
            A += [Vv[j] <= Vv[j - 1] + 8, Vv[j] >= Vv[0] - 8 * (DEPTH - j)]
    if mode == "yhint":
        # implied facts (each image of a sorted network is an ordered state below B); keeps the query equivalent
        for qq in qs[1:]:
            A += [qq[0] <= B - 1, qq[0] >= qq[1], qq[1] >= qq[2]]
    if k is not None:
        g = (B - u, u - v, v - w, w)
        for i in range(4):
            A.append(4 * g[i] >= B if (k >> i) & 1 else 4 * g[i] < B)
    return A


if __name__ == "__main__":
    mode = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] not in ("", "-") else None
    A = build(mode, k)
    if mode in ("yices", "yhint"):
        s = z3.Solver()
        for a in A:
            s.add(a)
        txt = s.to_smt2().replace("(check-sat)", "(set-logic QF_LIA)\n(check-sat)", 1)
        # put set-logic first
        txt = "(set-logic QF_LIA)\n" + txt.replace("(set-logic QF_LIA)\n", "")
        sys.stdout.write(txt)
        sys.exit(0)
    s = z3.Solver()
    tmo = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    s.set("timeout", tmo * 1000)
    for a in A:
        s.add(a)
    t = time.time()
    r = s.check()
    print("clock9 mode=%s part=%s -> %s (%.1fs)" % (mode, k, r, time.time() - t))
    if r == z3.sat:
        print(s.model())
