"""Independent Z3 encodings of the top-level computer-assisted claims of
"A common arrival-time law for five-digit Kaprekar dynamics" (P2).

Derived independently from the paper's text only, before
reading any project encoding.  The step map is encoded *differently* from the
paper's 26 guarded branches: the five borrowed digits of K_B(u,v) are put
through a bubble-sort comparator network (fresh integer variables, If-then-else
min/max), and the next pair is (s1-s5, s2-s4).  The comparator network is
checked to be a sorting network by the 0-1 principle, and the symbolic step is
cross-checked against direct numerical subtraction at concrete points.

Every claim is posed as "hypotheses AND NOT conclusion"; the expected answer is
UNSAT.  No upper bound on B is imposed anywhere.

Run:  kaprekar-investigation/.venv-research/bin/python \
        docs/drafts/p2-five-digit/review/z3_independent.py
"""
import itertools
import random
import sys
import time

import z3

COUNTER = itertools.count()


def fresh(prefix):
    return z3.Int(f"{prefix}_{next(COUNTER)}")


# ----------------------------------------------------------------------------
# Sorting network (descending) -- bubble sort, 10 comparators on 5 wires.
NETWORK = [(j, j + 1) for i in range(4) for j in range(4 - i)]


def network_is_sorting():
    for bits in itertools.product((0, 1), repeat=5):
        w = list(bits)
        for i, j in NETWORK:
            if w[i] < w[j]:
                w[i], w[j] = w[j], w[i]
        if any(w[k] < w[k + 1] for k in range(4)):
            return False
    return True


def step(s, B, u, v, tag):
    """Return (u', v') as fresh Int variables constrained in solver s."""
    pos = v > 0
    digits = [z3.If(pos, u, u - 1), z3.If(pos, v - 1, B - 1), B - 1,
              z3.If(pos, B - v - 1, B - 1), B - u]
    w = []
    for k, d in enumerate(digits):
        x = fresh(f"{tag}_d{k}")
        s.add(x == d)
        w.append(x)
    for i, j in NETWORK:
        hi, lo = fresh(f"{tag}_hi"), fresh(f"{tag}_lo")
        s.add(hi == z3.If(w[i] >= w[j], w[i], w[j]))
        s.add(lo == z3.If(w[i] >= w[j], w[j], w[i]))
        w[i], w[j] = hi, lo
    u2, v2 = fresh(f"{tag}_u"), fresh(f"{tag}_v")
    s.add(u2 == w[0] - w[4], v2 == w[1] - w[3])
    return u2, v2


def py_step_numeric(B, u, v):
    """Direct numerical subtraction from the definition (reference)."""
    N = (B**4 - 1) * u + (B**3 - B) * v       # = desc - asc for sorted pair
    d = []
    for _ in range(5):
        d.append(N % B)
        N //= B
    d.sort(reverse=True)
    return d[0] - d[4], d[1] - d[3]


def dom(B, u, v):
    return z3.And(1 <= u, u <= B - 1, 0 <= v, v <= u)


def H(B, u, v):
    a, b = u + v - B, 2 * u - v - B
    return z3.And(2 * v <= B - 1, 2 * (u - v) <= B - 3, a >= -2, b >= -2)


def J(B, u, v):
    return z3.And(2 * u >= B - 1, 2 * v <= B + 1, 2 * (u - v) <= B + 2)


def T2(B, m, u, v):
    # B = 2m in the caller when this is used; else empty
    return z3.Or(z3.And(u == m, v == m - 1), z3.And(u == m + 1, v == 0))


def margins(B, u, v):
    return u + v - B, 2 * u - v - B


RESULTS = []


def check(name, s, expect=z3.unsat):
    t0 = time.time()
    r = s.check()
    dt = time.time() - t0
    ok = (r == expect)
    extra = ""
    if r == z3.sat:
        m = s.model()
        extra = " model: " + ", ".join(
            f"{d.name()}={m[d]}" for d in m.decls() if "_" not in d.name())
    RESULTS.append((name, str(r), ok, dt))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {r} ({dt:.2f}s){extra}",
          flush=True)
    return r, s


def base_solver(lo=30):
    s = z3.Solver()
    B = z3.Int("B")
    s.add(B >= lo)
    return s, B


# ----------------------------------------------------------------------------
def sanity_step():
    """Cross-check the symbolic step against numerical subtraction."""
    rng = random.Random(20261117)
    bad = 0
    pts = []
    for _ in range(300):
        B = rng.randrange(2, 5000)
        u = rng.randrange(1, B)
        v = rng.choice([0, u, rng.randrange(0, u + 1)])
        pts.append((B, u, v))
    for B in (2, 3, 4, 30, 31):
        for u in range(1, B):
            for v in range(0, u + 1):
                pts.append((B, u, v))
    for (Bv, uv, vv) in pts:
        s = z3.Solver()
        B, u, v = z3.Ints("B u v")
        s.add(B == Bv, u == uv, v == vv)
        u2, v2 = step(s, B, u, v, "t")
        assert s.check() == z3.sat
        m = s.model()
        got = (m.eval(u2).as_long(), m.eval(v2).as_long())
        if got != py_step_numeric(Bv, uv, vv):
            bad += 1
    print(f"[{'PASS' if bad == 0 else 'FAIL'}] symbolic step == numeric "
          f"subtraction on {len(pts)} concrete points (mismatches {bad})")
    RESULTS.append(("step sanity", f"{len(pts)} pts", bad == 0, 0.0))


def borrowed_digit_identity():
    """(17): digits of K_B(u,v) -- polynomial identity + digit ranges."""
    s = z3.Solver()
    B, u, v = z3.Ints("B u v")
    s.add(B >= 2, dom(B, u, v))
    K = (B**4 - 1) * u + (B**3 - B) * v
    pos_val = u * B**4 + (v - 1) * B**3 + (B - 1) * B**2 + (B - v - 1) * B + (B - u)
    zero_val = (u - 1) * B**4 + (B - 1) * B**3 + (B - 1) * B**2 + (B - 1) * B + (B - u)
    s.add(z3.Or(z3.And(v > 0, pos_val != K), z3.And(v == 0, zero_val != K)))
    check("(17) borrowed-digit value identity", s)
    s = z3.Solver()
    s.add(B >= 2, dom(B, u, v), v > 0)
    ds = [u, v - 1, B - 1, B - v - 1, B - u]
    s.add(z3.Or([z3.Or(d < 0, d > B - 1) for d in ds]))
    check("(17) borrowed digits lie in [0,B-1] (v>0)", s)
    s = z3.Solver()
    s.add(B >= 2, dom(B, u, v), v == 0)
    ds = [u - 1, B - 1, B - 1, B - 1, B - u]
    s.add(z3.Or([z3.Or(d < 0, d > B - 1) for d in ds]))
    check("(17) borrowed digits lie in [0,B-1] (v=0)", s)


def claim_core_rule():
    """(20) on H_B and forward invariance of H_B."""
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v))
    u2, v2 = step(s, B, u, v, "c")
    a, b = margins(B, u, v)
    a2, b2 = margins(B, u2, v2)
    ra = z3.If(a >= 1, b, z3.If(a >= -1, b, b + 1))
    rb = z3.If(a >= 1, a - 3, z3.If(a >= -1, -2 * a, 3))
    s.add(z3.Not(z3.And(dom(B, u2, v2), H(B, u2, v2), a2 == ra, b2 == rb)))
    check("A1 core: f_B(H_B) in H_B and margins obey (20)", s)
    # also: in (u,v) coordinates as stated in the text
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v))
    u2, v2 = step(s, B, u, v, "c")
    a, b = margins(B, u, v)
    eu = z3.If(a >= 1, u - 1, B - v)
    ev = z3.If(a >= 1, u - v + 1, z3.If(a >= -1, 2 * u - B, u - v - 1))
    s.add(z3.Or(u2 != eu, v2 != ev))
    check("A1' core: (u,v)-form of the three rules", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v), z3.Or(2 * u < B, v <= 0))
    check("A1'' core implies 2u>=B and v>0", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    a, b = margins(B, u, v)
    s.add(dom(B, u, v), H(B, u, v), z3.Or(2 * a > B - 3, 2 * b > B - 3))
    check("A1''' on H_B, max(a,b) <= (B-3)/2", s)


def claim_strip_prep():
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), v > 0)
    u2, v2 = step(s, B, u, v, "p")
    s.add(z3.Not(z3.And(2 * u2 >= B - 1, 2 * (u2 - v2) <= B + 2,
                        dom(B, u2, v2))))
    check("A2 v>0 => 2u'>=B-1, 2(u'-v')<=B+2, valid pair", s)
    s, B = base_solver(2)
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v))
    u2, v2 = step(s, B, u, v, "p")
    s.add(z3.Not(dom(B, u2, v2)))
    check("A2' every pair maps to a valid pair (all B>=2)", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), v == 0)
    u2, v2 = step(s, B, u, v, "p")
    s.add(v2 != z3.If(u - 1 <= B - u, u - 1, B - u))
    check("A2'' v=0 => v' = min(u-1, B-u)", s)


def claim_high_region():
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), 2 * v >= B + 2, u > v)
    u2, v2 = step(s, B, u, v, "h")
    s.add(z3.Or(u2 != u - 1, v2 != u + v - B + 1))
    check("A3 high off-diagonal branch (u-1, u+v-B+1)", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), 2 * v >= B + 2, u == v)
    u2, v2 = step(s, B, u, v, "h")
    s.add(z3.Or(u2 != u, v2 != 2 * u - B))
    check("A3' high diagonal branch (u, 2u-B)", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), 2 * v >= B + 2, u > v)
    u2, v2 = step(s, B, u, v, "h")
    s.add(u2 == v2, z3.Not(z3.And(u == B - 1, v == B - 2)))
    check("A3'' high off-diag step makes a diagonal only from (B-1,B-2)", s)
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), 2 * v >= B + 2)
    u2, v2 = step(s, B, u, v, "h")
    s.add(z3.Not(z3.And(2 * u2 >= B - 1, 2 * (u2 - v2) <= B + 2, v2 > 0)))
    check("A3''' high steps keep 2u>=B-1, 2(u-v)<=B+2, v>0", s)
    # Theorem 5 uses the threshold v > H = floor((B+1)/2), i.e. 2v >= B+1
    s, B = base_solver()
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), 2 * v >= B + 1, u > v)
    u2, v2 = step(s, B, u, v, "h")
    s.add(z3.Or(u2 != u - 1, v2 != u + v - B + 1))
    check("A3'''' branch (u-1,u+v-B+1) already valid for 2v>=B+1, u>v", s)


def orbit(s, B, u0, v0, n, tag):
    qs = [(u0, v0)]
    for j in range(n):
        qs.append(step(s, B, qs[-1][0], qs[-1][1], f"{tag}{j}"))
    return qs


def claim_six_step(nsteps=6, expect=z3.unsat):
    s, B = base_solver()
    m, p = z3.Ints("m p")
    s.add(B == 2 * m + p, z3.Or(p == 0, p == 1))
    u0, v0 = z3.Ints("u0 v0")
    s.add(dom(B, u0, v0), J(B, u0, v0))
    qs = orbit(s, B, u0, v0, nsteps, "e")
    for (uj, vj) in qs:
        s.add(z3.Not(H(B, uj, vj)))
        s.add(z3.Not(z3.And(p == 0, T2(B, m, uj, vj))))
    return check(f"A4 J_B reaches H_B u T2 within {nsteps} steps", s, expect)


def rule_py(a, b):
    if a >= 1:
        return b, a - 3
    if a >= -1:
        return b, -2 * a
    assert a == -2
    return b + 1, 3


def claim_two_step():
    a, b = z3.Ints("a b")

    def rule(a, b):
        return (z3.If(a >= 1, b, z3.If(a >= -1, b, b + 1)),
                z3.If(a >= 1, a - 3, z3.If(a >= -1, -2 * a, 3)))
    s = z3.Solver()
    s.add(a >= -2, b >= -2)
    M = z3.If(a >= b, a, b)
    mn = z3.If(a >= b, b, a)
    s.add(M >= 10)
    a1, b1 = rule(a, b)
    a2, b2 = rule(a1, b1)
    M2 = z3.If(a2 >= b2, a2, b2)
    s.add(M2 != M - 3 + z3.If(mn == -2, 1, 0))
    check("A5 two-step descent identity (21)", s)
    # minimum stays >= -2 under the base-free rule
    s = z3.Solver()
    s.add(a >= -2, b >= -2)
    a1, b1 = rule(a, b)
    s.add(z3.Or(a1 < -2, b1 < -2))
    check("A5' base-free rule keeps both margins >= -2", s)
    # 'at most once': after an exceptional block the block-start minimum is
    # 3 or 0, and this is inductive while M >= 10.
    s = z3.Solver()
    s.add(a >= -2, b >= -2, M >= 10, mn == -2)
    s.add(z3.Not(z3.If(a2 >= b2, b2, a2) == 3))
    check("A5'' exceptional block outputs minimum 3", s)
    s = z3.Solver()
    s.add(a >= -2, b >= -2, M >= 10, z3.Or(mn == 3, mn == 0))
    mn2 = z3.If(a2 >= b2, b2, a2)
    s.add(z3.Not(z3.Or(mn2 == 0, mn2 == 3)))
    check("A5''' block-start minimum in {0,3} is inductive for M>=10", s)
    # intermediate (odd-time) state of a block with M>=10 has max > 3
    s = z3.Solver()
    s.add(a >= -2, b >= -2, M >= 10)
    s.add(z3.If(a1 >= b1, a1, b1) <= 3)
    check("A5'''' odd-time state inside an M>=10 block is not cyclic (max>3)", s)


def claim_rare_basin():
    s, B = base_solver()
    m = z3.Int("m")
    s.add(B == 2 * m)
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v))
    R = lambda x, y: z3.Or(z3.And(x == m, y == m - 1), z3.And(x == m + 1, y == 0),
                           z3.And(x == m, y == 0))
    s.add(z3.Not(R(u, v)))
    u2, v2 = step(s, B, u, v, "r")
    s.add(R(u2, v2))
    check("A6 even B: no pair outside R_B maps into R_B", s)
    # forward closure
    for (nm, x, y, ex, ey) in [("r1->r2", "m", "m-1", "m+1", "0"),
                               ("r2->r1", "m+1", "0", "m", "m-1"),
                               ("r3->r1", "m", "0", "m", "m-1")]:
        s, B = base_solver()
        m = z3.Int("m")
        s.add(B == 2 * m)
        env = {"m": m}
        X, Y = eval(x, {}, env), eval(y, {}, env)
        EX, EY = eval(ex, {}, env), eval(ey, {}, env)
        u2, v2 = step(s, B, X, Y, "f")
        s.add(z3.Or(u2 != EX, v2 != EY))
        check(f"A6' even B: {nm}", s)
    s, B = base_solver()
    m = z3.Int("m")
    s.add(B == 2 * m)
    s.add(z3.Or(H(B, m, m - 1), H(B, m + 1, 0), H(B, m, 0)))
    check("A6'' R_B is disjoint from H_B", s)


def claim_residue():
    s, B = base_solver()
    k = z3.Int("k")
    s.add(B == 3 * k + 1)
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v))
    u2, v2 = step(s, B, u, v, "q")
    s.add((u + v) % 3 == 0, (u2 + v2) % 3 != 0)
    check("A7 B=1 mod 3, H_B: u+v=0 => u'+v'=0 (mod 3)", s)
    s, B = base_solver()
    k = z3.Int("k")
    s.add(B == 3 * k + 1)
    u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v))
    u2, v2 = step(s, B, u, v, "q")
    s.add((u + v) % 3 != 0, (u2 + v2) % 3 == 0)
    check("A7' B=1 mod 3, H_B: u'+v'=0 => u+v=0 (mod 3)", s)


def claim_entry_margin(const=40):
    allok = True
    for j in range(7):
        s, B = base_solver()
        u0, v0 = z3.Ints("u0 v0")
        s.add(dom(B, u0, v0), J(B, u0, v0), 4 * u0 >= 3 * B)
        qs = orbit(s, B, u0, v0, j, f"m{j}_")
        for i in range(j):
            s.add(z3.Not(H(B, *qs[i])))
        s.add(H(B, *qs[j]))
        aj, bj = margins(B, *qs[j])
        s.add(z3.Or(2 * aj > 2 * u0 - B + const, 2 * bj > 2 * u0 - B + const))
        r, _ = check(f"A8 (36) first core entry at j={j}: "
                     f"2max(a,b) <= 2u0-B+{const}", s)
        allok &= (r == z3.unsat)
    return allok


def claim_overshoot_states():
    # (43) and exit state in H_B, with 16 <= R <= B/32
    s, B = base_solver()
    k, R = z3.Ints("k R")
    s.add(B == 2 * k, R >= 16, 32 * R <= B)
    u2, v2 = step(s, B, 2 * k - R, k, "o")
    s.add(z3.Not(z3.And(u2 == 2 * k - R - 1, v2 == k - R + 1, H(B, u2, v2))))
    check("A9 (43) even: (2k-R,k) -> (2k-R-1,k-R+1) in H_B", s)
    s, B = base_solver()
    s.add(B == 2 * k + 1, R >= 16, 32 * R <= B)
    u2, v2 = step(s, B, 2 * k + 1 - R, k + 1, "o")
    s.add(z3.Not(z3.And(u2 == 2 * k - R, v2 == k + 2 - R, H(B, u2, v2))))
    check("A9' (43) odd: (2k+1-R,k+1) -> (2k-R,k+2-R) in H_B", s)
    s, B = base_solver()
    nu, Hh = z3.Ints("nu Hh")
    s.add(R >= 16, 32 * R <= B, nu >= 1, nu <= R - 3,
          z3.Or(B == 2 * Hh, B == 2 * Hh - 1))
    s.add(z3.Not(z3.And(dom(B, B - R, Hh - nu), H(B, B - R, Hh - nu))))
    check("A9'' exit (B-R, H-nu) in H_B for 1<=nu<=R-3", s)


def box_enumeration():
    from fractions import Fraction as Fr
    cycles = set()
    maxh = 0
    lo, hi = None, None
    box_ok = True
    inter_max = -99
    for a0 in range(-2, 10):
        for b0 in range(-2, 10):
            # two-step map stays in box; intermediate coordinates
            a1, b1 = rule_py(a0, b0)
            a2, b2 = rule_py(a1, b1)
            inter_max = max(inter_max, a1, b1)
            if not (-2 <= a2 <= 9 and -2 <= b2 <= 9):
                box_ok = False
            seen = {}
            path = []
            st = (a0, b0)
            while st not in seen:
                seen[st] = len(path)
                path.append(st)
                st = rule_py(*st)
            h = seen[st]                      # first index on the cycle
            cyc = path[h:]
            k = cyc.index(min(cyc))
            cycles.add(tuple(cyc[k:] + cyc[:k]))
            maxh = max(maxh, h)
            d = Fr(h) - Fr(2, 3) * max(a0, b0)
            lo = d if lo is None else min(lo, d)
            hi = d if hi is None else max(hi, d)
    print(f"box: two-step stays in [-2,9]^2: {box_ok}; max intermediate "
          f"coordinate {inter_max}; max transient {maxh}; "
          f"h-2M/3 in [{lo}, {hi}]")
    for c in sorted(cycles, key=len):
        print("   cycle", c, "len", len(c), "a+b mod 3 =",
              {(x + y) % 3 for x, y in c})
    expected = {((0, 0),),
                ((-2, 0), (1, 3), (3, -2)),
                ((-1, -1), (-1, 2), (2, 2), (2, -1)),
                ((-1, 0), (0, 2), (2, 0), (0, -1))}
    ok = (box_ok and inter_max <= 10 and maxh == 10 and lo == -3
          and hi == Fr(16, 3) and cycles == expected)
    print(f"[{'PASS' if ok else 'FAIL'}] box enumeration matches (22) and "
          f"the four listed cycles")
    RESULTS.append(("box enumeration (22)", "enum", ok, 0.0))
    # residue criterion (31) on cycles: T4 class a=b=2 mod 3
    t4 = expected and ((-1, -1), (-1, 2), (2, 2), (2, -1))
    t3 = ((-2, 0), (1, 3), (3, -2))
    ok2 = all(a % 3 == 2 and b % 3 == 2 for a, b in t4) and \
        all(not (a % 3 == 2 and b % 3 == 2) for a, b in t3)
    print(f"[{'PASS' if ok2 else 'FAIL'}] T4 cycle inside (2,2) residue class,"
          f" T3 cycle outside")
    RESULTS.append(("cycle residue classes", "enum", ok2, 0.0))
    # cycle states lie in H_B for every B>=30 with matching residue
    s = z3.Solver()
    B = z3.Int("B")
    s.add(B >= 30)
    conds = []
    for c in expected:
        for (a, b) in c:
            # u=(a+b+2B)/3, v=(2a-b+B)/3; in H_B iff 4a-2b<=B-3, 4b-2a<=B-9
            conds.append(z3.Or(4 * a - 2 * b > B - 3, 4 * b - 2 * a > B - 9))
    s.add(z3.Or(conds))
    check("A10 every listed cycle state lies in H_B for all B>=30", s)


def leading_matrices():
    # (13) linear parts and determinant / residue permutation
    mats = {"x<1/2": ((0, -1), (-1, -1)), "x>=1/2,x+y<1": ((0, -1), (1, -1)),
            "x+y>=1,y<1/2": ((1, 0), (1, -1)), "y>=1/2": ((1, 0), (1, 1))}
    ok = True
    for nm, ((p, q), (r, t)) in mats.items():
        det = p * t - q * r
        imgs = {((p * x + q * y) % 3, (r * x + t * y) % 3)
                for x in range(3) for y in range(3)}
        ok &= abs(det) == 1 and len(imgs) == 9
        print(f"   leading cell {nm}: matrix {((p, q), (r, t))} det {det} "
              f"permutes (Z/3)^2: {len(imgs) == 9}")
    print(f"[{'PASS' if ok else 'FAIL'}] leading matrices unimodular")
    RESULTS.append(("leading matrices", "enum", ok, 0.0))


def tightness_probes():
    print("--- informational tightness probes (SAT expected) ---")
    r, s = claim_six_step(5, expect=z3.sat)
    if r == z3.sat:
        m = s.model()
        print("    5-step counterexample: B =", m[z3.Int("B")],
              "q0 =", (m[z3.Int("u0")], m[z3.Int("v0")]))
    # smallest constant that still works in (36) -- scan downwards
    for c in (40, 30, 20, 16, 12, 10, 8, 6, 4, 2, 0):
        ok = True
        for j in range(7):
            s, B = base_solver()
            u0, v0 = z3.Ints("u0 v0")
            s.add(dom(B, u0, v0), J(B, u0, v0), 4 * u0 >= 3 * B)
            qs = orbit(s, B, u0, v0, j, f"t{j}_")
            for i in range(j):
                s.add(z3.Not(H(B, *qs[i])))
            s.add(H(B, *qs[j]))
            aj, bj = margins(B, *qs[j])
            s.add(z3.Or(2 * aj > 2 * u0 - B + c, 2 * bj > 2 * u0 - B + c))
            if s.check() != z3.unsat:
                ok = False
                break
        print(f"    (36) with +{c//2} in place of +20: "
              f"{'holds' if ok else 'fails'}")
        if not ok:
            break


if __name__ == "__main__":
    print("z3", z3.get_version_string(), "python", sys.version.split()[0])
    assert network_is_sorting()
    print("[PASS] bubble network sorts all 0-1 inputs (0-1 principle)")
    sanity_step()
    borrowed_digit_identity()
    claim_core_rule()
    claim_strip_prep()
    claim_high_region()
    claim_six_step()
    claim_two_step()
    claim_rare_basin()
    claim_residue()
    claim_entry_margin()
    claim_overshoot_states()
    box_enumeration()
    leading_matrices()
    nfail = sum(1 for r in RESULTS if not r[2])
    print(f"=== {len(RESULTS)} checks, {nfail} failures ===")
    if "--probe" in sys.argv:
        tightness_probes()
    if "--probe" not in sys.argv:
        sys.exit(1 if nfail else 0)


# ----------------------------------------------------------------------------
# Non-vacuity and mutation probes: hypotheses alone must be SAT, and weakened
# or perturbed versions of the claims must be SAT (i.e. the encodings bite).
def nonvacuity():
    print("--- non-vacuity / mutation probes (each line states expectation) ---")
    out = []

    def probe(name, s, expect):
        r = s.check()
        ok = r == expect
        extra = ""
        if r == z3.sat:
            m = s.model()
            extra = " e.g. " + ", ".join(f"{d.name()}={m[d]}" for d in m.decls()
                                          if "_" not in d.name())
        print(f"   [{'ok' if ok else 'UNEXPECTED'}] {name}: {r}{extra}")
        out.append(ok)

    s, B = base_solver(); u, v = z3.Ints("u v")
    s.add(dom(B, u, v), H(B, u, v)); probe("H_B nonempty", s, z3.sat)
    s, B = base_solver(); m, p = z3.Ints("m p"); u, v = z3.Ints("u0 v0")
    s.add(B == 2 * m + p, z3.Or(p == 0, p == 1), dom(B, u, v), J(B, u, v),
          z3.Not(H(B, u, v)), z3.Not(z3.And(p == 0, T2(B, m, u, v))))
    probe("J_B minus (H_B u T2) nonempty", s, z3.sat)
    # six-step claim is tight: 5 steps do not suffice
    s, B = base_solver(); m, p = z3.Ints("m p"); u0, v0 = z3.Ints("u0 v0")
    s.add(B == 2 * m + p, z3.Or(p == 0, p == 1), dom(B, u0, v0), J(B, u0, v0))
    qs = orbit(s, B, u0, v0, 5, "z")
    for (uj, vj) in qs:
        s.add(z3.Not(H(B, uj, vj)), z3.Not(z3.And(p == 0, T2(B, m, uj, vj))))
    probe("5 steps do NOT suffice for J_B -> H_B u T2", s, z3.sat)
    # without T2 the six-step claim fails for even B
    s, B = base_solver(); u0, v0 = z3.Ints("u0 v0")
    s.add(dom(B, u0, v0), J(B, u0, v0))
    qs = orbit(s, B, u0, v0, 6, "y")
    for (uj, vj) in qs:
        s.add(z3.Not(H(B, uj, vj)))
    probe("dropping T2 from the target breaks six-step entry", s, z3.sat)
    # rare basin: removing r3 from R_B must be SAT (r3 maps in)
    s, B = base_solver(); m = z3.Int("m"); u, v = z3.Ints("u v")
    s.add(B == 2 * m, dom(B, u, v))
    R2 = lambda x, y: z3.Or(z3.And(x == m, y == m - 1), z3.And(x == m + 1, y == 0))
    s.add(z3.Not(R2(u, v)))
    u2, v2 = step(s, B, u, v, "w"); s.add(R2(u2, v2))
    probe("2-cycle alone is not predecessor-closed (r3 enters)", s, z3.sat)
    # each first-entry time j in (36) is realisable
    for j in range(7):
        s, B = base_solver(); u0, v0 = z3.Ints("u0 v0")
        s.add(dom(B, u0, v0), J(B, u0, v0), 4 * u0 >= 3 * B)
        qs = orbit(s, B, u0, v0, j, f"v{j}_")
        for i in range(j):
            s.add(z3.Not(H(B, *qs[i])))
        s.add(H(B, *qs[j]))
        probe(f"(36) hypotheses with first entry j={j} satisfiable", s, None)
    # residue claim hypotheses satisfiable
    s, B = base_solver(); k = z3.Int("k"); u, v = z3.Ints("u v")
    s.add(B == 3 * k + 1, dom(B, u, v), H(B, u, v), (u + v) % 3 == 0)
    probe("B=1 mod 3, H_B, u+v=0 mod 3 nonempty", s, z3.sat)
    # smallest bases: where do the claims start to hold?
    for lo in (2, 8, 10, 16, 20, 24, 26, 28, 29):
        s, B = base_solver(lo); u, v = z3.Ints("u v")
        s.add(dom(B, u, v), H(B, u, v))
        u2, v2 = step(s, B, u, v, "s")
        a, b = margins(B, u, v); a2, b2 = margins(B, u2, v2)
        ra = z3.If(a >= 1, b, z3.If(a >= -1, b, b + 1))
        rb = z3.If(a >= 1, a - 3, z3.If(a >= -1, -2 * a, 3))
        s.add(z3.Not(z3.And(dom(B, u2, v2), H(B, u2, v2), a2 == ra, b2 == rb)))
        r = s.check()
        mdl = s.model() if r == z3.sat else None
        print(f"   [info] core rule/invariance with B>={lo}: {r}"
              + (f" (B={mdl[B]})" if mdl else ""))
    return all(x for x in out if x is not None)


if __name__ == "__main__" and "--probe" in sys.argv:
    nonvacuity()
