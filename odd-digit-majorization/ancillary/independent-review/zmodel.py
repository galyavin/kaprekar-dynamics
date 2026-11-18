"""z3 encodings of the transition model, written independently.

Integer model: digits by column subtraction (lemma Z1 proves it equals the
borrow-formula digit list for every rank), then a Batcher sorting network of
max/min comparators on all nine digits, then symmetric differences.
Real model: the leading map f via a 9-element sorting network, and the
paper's radius map H via a 4-element sorting network (lemma L1 proves
C^{-1} f C = H).
"""
import z3
from km import batcher_network, check_network, CORE

NET = {n: batcher_network(n) for n in (4, 5, 8, 9)}
for n, net in NET.items():
    assert check_network(net, n)


def zmax(a, b):
    return z3.If(a >= b, a, b)


def zmin(a, b):
    return z3.If(a >= b, b, a)


def sort_desc(vals):
    v = list(vals)
    for i, j in NET[len(v)]:
        a, b = v[i], v[j]
        v[i], v[j] = zmax(a, b), zmin(a, b)
    return v


def zabs(x):
    return z3.If(x >= 0, x, -x)


# ------------------------------------------------------------ integer model
def col_sub_digits(q, B):
    """Column subtraction desc-asc for the sorted realization
    (q1,q2,q3,q4,q4,0,0,0,0); returns digits most-significant first."""
    desc = [q[0], q[1], q[2], q[3], q[3], 0, 0, 0, 0]
    asc = desc[::-1]
    out = [None] * 9
    borrow = z3.IntVal(0)
    for j in range(8, -1, -1):
        t = desc[j] - asc[j] - borrow
        out[j] = z3.If(t < 0, t + B, t)
        borrow = z3.If(t < 0, z3.IntVal(1), z3.IntVal(0))
    return out, borrow


def borrow_digits(q, B):
    """Borrow-formula digits, all ranks (rank r: lead q1..q_{r-1},q_r-1,
    9-2r copies of B-1, then B-1-q_r..B-1-q_2, B-q_1)."""
    u, v, w, z = q
    r4 = [u, v, w, z - 1, B - 1, B - 1 - z, B - 1 - w, B - 1 - v, B - u]
    r3 = [u, v, w - 1, B - 1, B - 1, B - 1, B - 1 - w, B - 1 - v, B - u]
    r2 = [u, v - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - 1 - v, B - u]
    r1 = [u - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - u]
    r0 = [z3.IntVal(0)] * 9
    out = []
    for i in range(9):
        out.append(z3.If(z > 0, r4[i], z3.If(w > 0, r3[i], z3.If(v > 0, r2[i],
                   z3.If(u > 0, r1[i], r0[i])))))
    return out


def Fsym(q, B):
    d = sort_desc(borrow_digits(q, B))
    return [d[i] - d[8 - i] for i in range(4)]


def gaps_sym(q, B):
    u, v, w, z = q
    return [B - u, u - v, v - w, w - z, z]


def S_sym(g, j):
    s = sort_desc(g)
    return sum(s[:j])


def valid_q(q, B):
    u, v, w, z = q
    return z3.And(B - 1 >= u, u >= v, v >= w, w >= z, z >= 0)


# ------------------------------------------------------------ real model
def f_sym(x):
    one = z3.RealVal(1)
    vals = [one] + list(x) + [one - t for t in x]
    d = sort_desc(vals)
    return [d[i] - d[8 - i] for i in range(4)]


def C_sym(r):
    a, b, c, d = r
    return [(1 + a) / 2, (a + b) / 2, (b + c) / 2, (c + d) / 2]


def Cinv_sym(y):
    u, v, w, z = y
    a = 2 * u - 1
    b = 2 * v - a
    c = 2 * w - b
    d = 2 * z - c
    return [a, b, c, d]


def H_sym(r):
    a, b, c, d = r
    return sort_desc([a, zabs(a + b - 1), zabs(b + c - 1), zabs(c + d - 1)])


def D_sym(r):
    a, b, c, d = r
    return [1 - a, 1 - b, a - c, b - d, c + d]


def core_rows(r):
    return [row[0] * r[0] + row[1] * r[1] + row[2] * r[2] + row[3] * r[3] + row[4] for row in CORE]


def in_core(r):
    return z3.And([s >= 0 for s in core_rows(r)])


def R_sym(y):
    u, v, w, z = y
    return [u, u - z, v - z, v - w]


def in_E(r):
    a, b, c, d = r
    rng = z3.And(a >= z3.Q(5, 7), a <= z3.Q(3, 4))
    p1 = z3.And(rng, b == 2 * a - 1, c == 2 * a - 1, d > 0, d <= 1 - a)
    p2 = z3.And(rng, c == 3 * a - 2, d == 3 - 4 * a, b >= 1 - a, b < 2 - 2 * a)
    p3 = z3.And(rng, b == a, c + d == 2 - 2 * a, d >= 0, d < 1 - a)
    return z3.Or(p1, p2, p3)


def eqv(r, pt):
    return z3.And([r[i] == pt[i] for i in range(4)])


def W_sym(r):
    m = core_rows(r)
    mn = m[0]
    for t in m[1:]:
        mn = zmin(mn, t)
    return zmax(z3.RealVal(0), -mn)


def run(solver_name, s, timeout_ms=None):
    if timeout_ms:
        s.set("timeout", timeout_ms)
    r = s.check()
    return r
