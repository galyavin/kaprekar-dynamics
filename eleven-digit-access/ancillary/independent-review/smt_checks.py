"""Independent SMT checks, decided by cvc5 1.4.2 (Python API) and Yices 2.7.0 (CLI).
cvc5 was not used by the authors.  All queries are QF_LRA (real-linear) unless noted.
  S1  eq.(paired-return) == Psi o Psi under the two retained guards        (expect unsat)
  S2  one-unit transfer identity z+ = Tz + clip(d1,-1,1)(0,-1,1,0,0,0)        (expect unsat)
  S3  Lemma 12 first claim with allowance 6 and 12: min ghat(Psi X) >= min ghat(X) - 6 (expect unsat)
  S4  Lemma 9(2) ray: c>=d, |xi|=1, |F_S^10 xi| <= 2 (expect unsat) + bisection of the sharp factor
  S5  Lemma 9(4) isolated cycle with T_I built from the radius map rho (own derivation, not the
      signed words): |x|=1, b>=c, a<=0 and -b+c+2d+h>=0 at x and T_I x       (expect unsat)
Install cvc5 1.4.2 in the Python environment.  Usage: python smt_checks.py OUTDIR
"""
import os
import subprocess
import sys
from fractions import Fraction as Fr

import cvc5


class Q:
    def __init__(self, logic="QF_LRA"):
        self.lines = [f"(set-logic {logic})"]
        self.n = 0

    def var(self, name=None, sort="Real"):
        if name is None:
            self.n += 1
            name = f"t{self.n}"
        self.lines.append(f"(declare-fun {name} () {sort})")
        return name

    def let(self, expr):
        v = self.var()
        self.lines.append(f"(assert (= {v} {expr}))")
        return v

    def a(self, expr):
        self.lines.append(f"(assert {expr})")

    def text(self):
        return "\n".join(self.lines + ["(check-sat)", "(exit)"]) + "\n"


def num(x):
    x = Fr(x)
    if x.denominator == 1:
        return str(x.numerator) if x >= 0 else f"(- {-x.numerator})"
    n, d = x.numerator, x.denominator
    return f"(/ {n} {d})" if n >= 0 else f"(- (/ {-n} {d}))"


def mx(q, a, b): return q.let(f"(ite (>= {a} {b}) {a} {b})")
def mn(q, a, b): return q.let(f"(ite (<= {a} {b}) {a} {b})")
def ab(q, a): return q.let(f"(ite (>= {a} 0) {a} (- {a}))")
def ps(q, a): return q.let(f"(ite (> {a} 0) {a} 0)")
def add(q, *xs): return q.let("(+ " + " ".join(xs) + ")") if len(xs) > 1 else xs[0]
def sub(q, a, b): return q.let(f"(- {a} {b})")
def mul(q, c, a): return q.let(f"(* {num(c)} {a})")
def neg(q, a): return q.let(f"(- {a})")


def psi(q, X):
    c = [num(-6)]
    for i in range(1, 6):
        a_ = sub(q, X[i - 1], num(6)) if i == 5 else X[i - 1]
        b_ = sub(q, neg(q, X[6 - i - 1]), num(6)) if i < 5 else neg(q, X[6 - i - 1])
        c += [mx(q, a_, b_), mn(q, a_, b_)]
    return [sub(q, c[k - 1], c[11 - k]) for k in range(1, 6)]


def ghat(q, X):
    return [neg(q, X[0])] + [sub(q, X[i], X[i + 1]) for i in range(4)] + [X[4]]


def S1():
    q = Q()
    l, w, p, s, r = (q.var(n) for n in ("l", "w", "p", "s", "r"))
    X = [neg(q, l), add(q, mul(q, -2, l), mul(q, 6, p)), add(q, mul(q, -3, l), mul(q, 6, s)),
         add(q, mul(q, 6, w), mul(q, -4, l)), add(q, l, mul(q, 6, r))]
    Y = psi(q, psi(q, X))
    l2 = neg(q, Y[0])
    p2 = mul(q, Fr(1, 6), add(q, Y[1], mul(q, 2, l2)))
    s2 = mul(q, Fr(1, 6), add(q, Y[2], mul(q, 3, l2)))
    w2 = mul(q, Fr(1, 6), add(q, Y[3], mul(q, 4, l2)))
    r2 = mul(q, Fr(1, 6), sub(q, Y[4], l2))
    Cf = ps(q, add(q, neg(q, l), mul(q, 2, s), num(1)))
    Af = mx(q, num(-1), neg(q, r))
    J = sub(q, ps(q, sub(q, neg(q, r), num(1))), ps(q, sub(q, r, num(1))))
    Kf = ps(q, add(q, r, num(1)))
    Hf = ps(q, sub(q, sub(q, sub(q, w, J), p), num(1)))
    t = sub(q, sub(q, s, p), Cf)
    Df = mx(q, sub(q, Af, num(1)), t)
    Ef = sub(q, mx(q, Af, sub(q, t, num(1))), mx(q, sub(q, neg(q, t), num(1)), neg(q, Af)))
    Ff = add(q, mn(q, Af, sub(q, t, num(1))), J, num(1))
    fl = sub(q, l, mul(q, 6, Df))
    fw = sub(q, sub(q, add(q, w, J), p), mul(q, 4, Df))
    fp = sub(q, Ef, mul(q, 2, Df))
    fs = add(q, sub(q, Ff, mul(q, 3, Df)), Hf)
    fr = sub(q, add(q, sub(q, neg(q, p), Kf), Df), Hf)
    q.a(f"(<= (+ (- {w} {l}) {p}) (- 1))")
    q.a(f"(>= (- (- {l} (* 2 {w})) (* 2 {Kf})) (- 1))")
    q.a(f"(or (distinct {l2} {fl}) (distinct {w2} {fw}) (distinct {p2} {fp}) (distinct {s2} {fs}) (distinct {r2} {fr}))")
    return q.text()


def S2():
    q = Q()
    X = [q.var(f"X{i}") for i in range(1, 6)]
    u = [num(0)] + [mul(q, Fr(1, 6), x) for x in X] + [num(0)]
    d = [add(q, u[i], u[6 - i]) for i in range(4)]
    z = [mul(q, Fr(1, 6), g) for g in ghat(q, X)]
    T = [None] * 6
    for i in range(3):
        m = sub(q, ab(q, add(q, d[i], num(1))), ab(q, add(q, d[i + 1], num(1))))
        zz = add(q, z[i], z[5 - i])
        T[2 * i] = mul(q, Fr(1, 2), add(q, zz, m))
        T[2 * i + 1] = mul(q, Fr(1, 2), sub(q, zz, m))
    m0 = mn(q, z[5], add(q, z[0], num(1)))
    T[0] = m0
    T[1] = sub(q, add(q, z[0], z[5]), m0)
    cl = mx(q, num(-1), mn(q, num(1), d[1]))
    pred = [T[0], sub(q, T[1], cl), add(q, T[2], cl), T[3], T[4], T[5]]
    zp = [mul(q, Fr(1, 6), g) for g in ghat(q, psi(q, X))]
    q.a("(or " + " ".join(f"(distinct {a} {b})" for a, b in zip(zp, pred)) + ")")
    return q.text()


def S3(allow):
    q = Q()
    X = [q.var(f"X{i}") for i in range(1, 6)]
    m = q.var("m")
    for g in ghat(q, X):
        q.a(f"(>= {g} {m})")
    q.a("(or " + " ".join(f"(< {g} (- {m} {allow}))" for g in ghat(q, psi(q, X))) + ")")
    return q.text()


def unit_norm(q, xs):
    for x in xs:
        q.a(f"(<= {x} 1)")
        q.a(f"(>= {x} (- 1))")
    q.a("(or " + " ".join(f"(= {x} 1) (= {x} (- 1))" for x in xs) + ")")


def S4(K):
    q = Q()
    v, c, d, z = (q.var(n) for n in ("v", "c", "d", "z"))
    q.a(f"(>= {c} {d})")
    unit_norm(q, [v, c, d, z])
    xi = [v, c, d, z]
    for _ in range(10):
        v_, c_, d_, z_ = xi
        ncd = sub(q, neg(q, c_), d_)
        xi = [sub(q, neg(q, d_), z_), mx(q, v_, ncd), mn(q, v_, ncd), sub(q, neg(q, v_), c_)]
    for x in xi:
        q.a(f"(<= {x} {num(K)})")
        q.a(f"(>= {x} {num(-K)})")
    return q.text()


def rho(q, r):
    raw = [r[0]] + [ab(q, sub(q, add(q, r[i], r[i + 1]), num(1))) for i in range(4)]
    a = list(raw)
    for i in range(4):            # bubble network, descending
        for j in range(4 - i):
            hi, lo = mx(q, a[j], a[j + 1]), mn(q, a[j], a[j + 1])
            a[j], a[j + 1] = hi, lo
    return a


def S5(eps):
    q = Q()
    x = [q.var(n) for n in ("a", "b", "c", "d", "h")]
    rA = [Fr(4, 5), Fr(3, 5), Fr(3, 5), Fr(2, 5), Fr(1, 5)]
    unit_norm(q, x)
    q.a(f"(>= {x[1]} {x[2]})")
    r = [add(q, num(rA[i]), mul(q, eps, x[i])) for i in range(5)]
    for _ in range(3):
        r = rho(q, r)
    xp = [mul(q, 1 / Fr(eps), sub(q, r[i], num(rA[i]))) for i in range(5)]
    for (a_, b_, c_, d_, h_) in (x, xp):
        q.a(f"(<= {a_} 0)")
        q.a(f"(>= (+ (- {c_} {b_}) (* 2 {d_}) {h_}) 0)")
    return q.text()


def S5_return_samples(eps):
    """Also print T_I at the unit vectors (for the report)."""
    pass


def run_cvc5(txt):
    tm = cvc5.TermManager()
    s = cvc5.Solver(tm)
    sm = cvc5.SymbolManager(tm)
    ip = cvc5.InputParser(s, sm)
    ip.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6, txt, "q")
    res = None
    while True:
        cmd = ip.nextCommand()
        if cmd.isNull():
            break
        out = cmd.invoke(s, sm)
        if out.strip() in ("sat", "unsat", "unknown"):
            res = out.strip()
    return res


def run_yices(path):
    r = subprocess.run(["yices-smt2", path], capture_output=True, text=True, timeout=3000)
    return r.stdout.strip().splitlines()[0] if r.stdout.strip() else r.stderr.strip()[:200]


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    jobs = [("S1_paired_return", S1(), "unsat"), ("S2_transfer_identity", S2(), "unsat"),
            ("S3_allow6", S3(6), "unsat"), ("S3_allow12", S3(12), "unsat"), ("S3_allow5_control", S3(5), "sat"),
            ("S4_ray_K2", S4(2), "unsat"),
            ("S5_isolated_eps1e-4", S5(Fr(1, 10**4)), "unsat"), ("S5_isolated_eps1e-6", S5(Fr(1, 10**6)), "unsat")]
    lines = []
    for name, txt, exp in jobs:
        path = os.path.join(outdir, name + ".smt2")
        open(path, "w").write(txt)
        rc = run_cvc5(txt)
        ry = run_yices(path)
        lines.append(f"{name}: cvc5={rc} yices={ry} expected={exp} {'OK' if rc == ry == exp else 'CHECK'}")
        print(lines[-1], flush=True)
    # bisection for the sharp ten-step ray factor
    lo, hi = Fr(2), Fr(6)
    for _ in range(40):
        mid = (lo + hi) / 2
        if run_cvc5(S4(mid)) == "unsat":
            lo = mid
        else:
            hi = mid
    lines.append(f"S4 bisection: min |F_S^10 xi|_inf over c>=d,|xi|=1 lies in [{float(lo):.9f}, {float(hi):.9f}]; "
                 f"295/77={295/77:.9f}; at 295/77: cvc5={run_cvc5(S4(Fr(295, 77)))} (sat expected if attained), "
                 f"just below: cvc5={run_cvc5(S4(Fr(295, 77) - Fr(1, 10**9)))}")
    print(lines[-1], flush=True)
    open(os.path.join(outdir, "smt_results.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
