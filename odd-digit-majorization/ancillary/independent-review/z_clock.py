"""Independent z3 encodings for the integer statements of Appendix C and the
bulk motion of Section 4.1 (nine digits).  Transition: borrow-formula digits
(proved equal to column subtraction by Z1), Batcher 9-sorting network,
symmetric differences.  All variables are unbounded integers; no upper bound
on B is imposed anywhere.

Usage: z_clock.py [names...]   (default: all)
"""
import sys, time, z3
from zmodel import *

SEL = set(sys.argv[1:])


def want(n):
    return not SEL or n in SEL


def report(name, s, expect=z3.unsat):
    t = time.time()
    res = s.check()
    dt = time.time() - t
    ok = res == expect
    print(f"{name:62s} {str(res):8s} {'OK' if ok else '**FAIL**'} {dt:8.1f}s", flush=True)
    if res == z3.sat:
        print("   model:", s.model())
    return ok


def qvars(tag):
    return [z3.Int(f"{tag}{i}") for i in range(4)]


allok = True
B = z3.Int("B")

# Z1: column subtraction == borrow formula, all ranks
if want("Z1"):
    s = z3.Solver()
    q = qvars("q")
    s.add(B >= 2, valid_q(q, B))
    cs, br = col_sub_digits(q, B)
    bd = borrow_digits(q, B)
    s.add(z3.Or(br != 0, z3.Or([cs[i] != bd[i] for i in range(9)])))
    allok &= report("Z1 column subtraction == borrow digits (all ranks, B>=2)", s)

# Z2: S_i(g(Fq)) <= S_i(g(q)) + 1, q4 > 0, B >= 2
if want("Z2"):
    for i in range(1, 5):
        s = z3.Solver()
        q = qvars("q")
        s.add(B >= 2, valid_q(q, B), q[3] >= 1)
        q1 = Fsym(q, B)
        s.add(S_sym(gaps_sym(q1, B), i) > S_sym(gaps_sym(q, B), i) + 1)
        allok &= report(f"Z2 S_{i}(g(Fq)) <= S_{i}(g(q))+1  (q4>0, B>=2)", s)

# Z3: first access: four consecutive positive states force a positive fifth (B>=40)
def first_access(s, Bmin, parity=None):
    q0 = qvars("q")
    s.add(B >= Bmin, valid_q(q0, B), q0[3] >= 1)
    if parity is not None:
        s.add(B % 2 == parity)
    qs = [q0]
    for t in range(4):
        qs.append(Fsym(qs[-1], B))
    for t in range(1, 4):
        s.add(qs[t][3] >= 1)
    s.add(qs[4][3] == 0)


if want("Z3"):
    for par in (0, 1):
        s = z3.Solver()
        first_access(s, 40, par)
        allok &= report(f"Z3 q0..q3 positive => q4 positive  (B>=40, B%2={par})", s)
    s = z3.Solver()
    first_access(s, 33, 1)
    allok &= report("Z3b same for odd B>=33 (paper: a>=16)", s)
    # even B, stronger statement in the text: positive->positive->zero impossible
    s = z3.Solver()
    q0 = qvars("q")
    s.add(B >= 40, B % 2 == 0, valid_q(q0, B), q0[3] >= 1)
    q1 = Fsym(q0, B)
    q2 = Fsym(q1, B)
    s.add(q1[3] >= 1, q2[3] == 0)
    allok &= report("Z3c even B>=40: positive->positive->zero impossible", s)


# Z4: cover of returning zero states
def ab(x):
    return z3.If(x >= 0, x, -x)


def cover2(q, B):
    u, v, w, z = q
    return z3.Or(ab(u - v) <= 3, ab(u + v - B) <= 3, ab(2 * u - B) <= 4, ab(2 * v - B) <= 4, B - u <= 3)


def cover3(q, B):
    u, v, w, z = q
    row1 = z3.And(v == w, z3.Or(ab(u - v) <= 2, ab(u + v - B) <= 2, ab(2 * v - B) <= 3))
    row2 = z3.And(u == v, z3.Or(ab(u + w - B) <= 2, ab(2 * u - B) <= 3))
    row3 = z3.And(ab(2 * u - B) <= 3, ab(2 * v - B) <= 3)
    row4 = z3.And(B - u <= 3, z3.Or(ab(2 * v - B) <= 3, v == w, ab(v + w - B) <= 3))
    row5 = z3.And(ab(4 * u - 3 * B) <= 8, z3.Or([ab(4 * y - j * B) <= 12 for y in (v, w) for j in (1, 2, 3)]))
    return z3.Or(row1, row2, row3, row4, row5)


if want("Z4"):
    for rank in (2, 3):
        for j in (1, 2, 3, 4):
            for par in (0, 1):
                s = z3.Solver()
                q0 = qvars("q")
                s.add(B >= 40, B % 2 == par, valid_q(q0, B), q0[3] == 0, q0[0] > 0)
                if rank == 2:
                    s.add(q0[1] > 0, q0[2] == 0)
                    s.add(z3.Not(cover2(q0, B)))
                else:
                    s.add(q0[2] > 0)
                    s.add(z3.Not(cover3(q0, B)))
                qs = [q0]
                for t in range(j):
                    qs.append(Fsym(qs[-1], B))
                for t in range(1, j):
                    s.add(qs[t][3] >= 1)
                s.add(qs[j][3] == 0)
                allok &= report(f"Z4 rank {rank} zero state, first return at j={j}, B%2={par}, outside cover", s)

# Z5: bulk three-step translation by (-3,-2,-1,-1) on the core with integer margin K
def core_int_slacks(q, B):
    u, v, w, z = q
    a = 2 * u - B
    b = 2 * v - 2 * u + B
    c = 2 * w - 2 * v + 2 * u - B
    d = 2 * z - 2 * w + 2 * v - 2 * u + B
    from km import CORE
    return [row[0] * a + row[1] * b + row[2] * c + row[3] * d + row[4] * B for row in CORE]


if want("Z5"):
    for K in (40, 24, 16):
        s = z3.Solver()
        q = qvars("q")
        s.add(B >= 2, valid_q(q, B))
        for t in core_int_slacks(q, B):
            s.add(t >= K)
        q3 = q
        for _ in range(3):
            q3 = Fsym(q3, B)
        dlt = (-3, -2, -1, -1)
        s.add(z3.Or([q3[i] != q[i] + dlt[i] for i in range(4)]))
        allok &= report(f"Z5 F^3 q = q+(-3,-2,-1,-1) when all 12 integer core slacks >= {K}", s)

print("ALL OK" if allok else "SOME FAILED")
