"""Lemma A.2 / Proposition 4 (edge clock), own derivation.

State at an edge block start, base B = 2b+beta, offsets (m,S) with m = beta mod 2:
    u free,  v = (B+m)/2 = b + (beta+m)/2,  w = B - u + S.
All digits are affine forms in (u, b, 1).  For each node (m,S) the three-step block is
computed symbolically: at the representative point (B=100000+beta, u=80000) the seven digits
are sorted by value; the resulting order is then certified over the whole strip
    R = {B >= 100000, B = beta mod 2, 3B/4+128 <= u <= 5B/6-128}
by checking every adjacent comparison (non-strict suffices: sorting is by value) and every
intermediate positivity w>=1 on the real polyhedron R (two vertices + two recession rays).
Only non-strict comparisons are needed because sorted VALUES are what enter the differences.
"""
from fractions import Fraction as Fr
from collections import defaultdict
import sys

REP_B0 = 100000
REP_U = 80000


class Aff:
    __slots__ = ("cu", "cb", "c0")

    def __init__(self, cu, cb, c0):
        self.cu, self.cb, self.c0 = Fr(cu), Fr(cb), Fr(c0)

    def __add__(s, o):
        o = o if isinstance(o, Aff) else Aff(0, 0, o)
        return Aff(s.cu + o.cu, s.cb + o.cb, s.c0 + o.c0)

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, Aff) else Aff(0, 0, o)
        return Aff(s.cu - o.cu, s.cb - o.cb, s.c0 - o.c0)

    def __rsub__(s, o):
        return Aff(0, 0, o) - s

    def __mul__(s, k):
        return Aff(s.cu * k, s.cb * k, s.c0 * k)

    __rmul__ = __mul__

    def at(s, u, b):
        return s.cu * u + s.cb * b + s.c0

    def key(s):
        return (s.cu, s.cb, s.c0)

    def __repr__(s):
        return "%s*u+%s*b+%s" % (s.cu, s.cb, s.c0)


def region_vertices_rays(beta, lo=128, hi=128, Bmin=REP_B0):
    """Real polyhedron in (u,b): 2b+beta>=Bmin, u>=3B/4+lo, u<=5B/6-hi."""
    b0 = Fr(Bmin - beta, 2)
    B0 = 2 * b0 + beta
    verts = [(Fr(3, 4) * B0 + lo, b0), (Fr(5, 6) * B0 - hi, b0)]
    # rays: b += 1 (B += 2): u += 3/2 or 5/3
    rays = [(Fr(3, 2), Fr(1)), (Fr(5, 3), Fr(1))]
    return verts, rays


def nonneg_on(g, verts, rays):
    if any(g.at(u, b) < 0 for (u, b) in verts):
        return False
    if any(g.cu * du + g.cb * db < 0 for (du, db) in rays):
        return False
    return True


def borrow(u, v, w, B):
    return [u, v, w - 1, B - 1, B - w - 1, B - v - 1, B - u]


def block(beta, m, S, lo=128, hi=128, Bmin=REP_B0):
    """Return (m', S', du_final, [du1, du2], nguards, failures) for one three-step block."""
    B = Aff(0, 2, beta)
    u = Aff(1, 0, 0)
    assert (m - beta) % 2 == 0
    v = Aff(0, 1, Fr(beta + m, 2))
    w = B - u + S
    q = (u, v, w)
    verts, rays = region_vertices_rays(beta, lo, hi, Bmin)
    repb = Fr(REP_B0 + beta - beta, 2)  # b at representative B=100000+beta
    repb = Fr(REP_B0, 2)
    ur = REP_U
    guards = []
    us = []
    for step in range(3):
        uu, vv, ww = q
        guards.append(ww - 1)  # positive state: borrow formula valid
        guards.append(uu - vv)
        guards.append(vv - ww)
        guards.append((B - 1) - uu)
        d = borrow(uu, vv, ww, B)
        vals = [x.at(ur, repb) for x in d]
        order = sorted(range(7), key=lambda i: -vals[i])
        for i in range(6):
            guards.append(d[order[i]] - d[order[i + 1]])
        s = [d[i] for i in order]
        q = (s[0] - s[6], s[1] - s[5], s[2] - s[4])
        us.append(q[0])
    uu, vv, ww = q
    fails = [g for g in guards if not nonneg_on(g, verts, rays)]
    m2 = 2 * vv - B
    S2 = uu + ww - B
    du = uu - u
    for x in (m2, S2, du, us[0] - u, us[1] - u):
        if x.cu != 0 or x.cb != 0:
            fails.append(("nonconstant", x))
    return (int(m2.c0), int(S2.c0), int(du.c0), [int((us[0] - u).c0), int((us[1] - u).c0)], len(guards), fails)


def build_graph(beta, R=20, **kw):
    start = [(m, S) for m in range(-R, R + 1) for S in range(-R, R + 1) if (m - beta) % 2 == 0]
    rule = {}
    todo = list(start)
    nguards_distinct = set()
    allfails = []
    while todo:
        node = todo.pop()
        if node in rule:
            continue
        m2, S2, du, partial, ng, fails = block(beta, *node, **kw)
        if fails:
            allfails.append((node, fails))
        rule[node] = ((m2, S2), du, partial)
        todo.append((m2, S2))
    return start, rule, allfails


def analyse(beta, start, rule):
    # cycles
    cycles = {}
    node_cycle = {}
    for n0 in rule:
        path = []
        seen = {}
        n = n0
        while n not in seen:
            seen[n] = len(path)
            path.append(n)
            n = rule[n][0]
        cyc = path[seen[n]:]
        key = min(cyc)
        if key not in cycles:
            # rotate to start at key
            i = cyc.index(key)
            cyc = cyc[i:] + cyc[:i]
            ell = len(cyc)
            D = -sum(rule[c][1] for c in cyc)
            cycles[key] = (cyc, ell, D)
    # transient data from the starting square
    maxk = 0
    maxdisp = 0
    maxpartial = 0
    entries = []
    for n0 in start:
        n = n0
        k = 0
        disp = 0
        oncyc = None
        visited = []
        while True:
            # is n on a cycle?
            found = None
            for key, (cyc, ell, D) in cycles.items():
                if n in cyc:
                    found = key
            if found is not None:
                break
            (n2, du, partial) = rule[n]
            maxpartial = max(maxpartial, max(abs(p) for p in partial), abs(du))
            disp += du
            maxdisp = max(maxdisp, abs(disp))
            n = n2
            k += 1
        e = disp
        cyc, ell, D = cycles[found]
        # one further circuit from entry phase n
        d2 = disp
        x = n
        for _ in range(ell):
            (x2, du, partial) = rule[x]
            maxpartial = max(maxpartial, max(abs(p) for p in partial), abs(du))
            d2 += du
            maxdisp = max(maxdisp, abs(d2))
            x = x2
        assert x == n
        maxk = max(maxk, k)
        entries.append((n0, k, e, found, n))
    return cycles, maxk, maxdisp, maxpartial, entries


if __name__ == "__main__":
    out = {}
    for beta in (0, 1):
        start, rule, fails = build_graph(beta)
        cycles, maxk, maxdisp, maxpartial, entries = analyse(beta, start, rule)
        print("beta=%d initial nodes=%d reachable nodes=%d guard failures=%d" % (beta, len(start), len(rule), len(fails)))
        for f in fails[:5]:
            print("   FAIL", f)
        for key, (cyc, ell, D) in sorted(cycles.items()):
            print("   cycle", cyc, "ell=%d D=%d r=3l/D=%s" % (ell, D, Fr(3 * ell, D)))
        print("   max transient blocks=%d  max |cum. displacement| (transient+1 circuit)=%d  max partial/blk |du|=%d" % (maxk, maxdisp, maxpartial))
        out[beta] = (start, rule, cycles, entries)
