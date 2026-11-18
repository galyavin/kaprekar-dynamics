"""Fresh-base tests of Proposition 4 (edge clock) and Proposition 5 (facet transit).

T is computed by word-level digit subtraction (./traj diffT); h = T - 1 (word transient).
Rates come from my own edge graph (edge.py) and my own facet bookkeeping.
Edge: the proof gives h - r(u-3B/4) = 3k + r e + h_end - r(u_end-3B/4), which Lemma A.2
asserts lies in [-43/2, 343/8]; we test that sharper interval, not only |.|<=200.
Facet: |T - B t_r(q/B)| <= 216 with r from the residue rule; also report the distance to the
other clock (must be large, so the rate assignment is observable)."""
import random, subprocess, sys
from fractions import Fraction as Fr
from edge import build_graph, analyse

rates = {}
for beta in (0, 1):
    start, rule, fails = build_graph(beta)
    cycles, *_ = analyse(beta, start, rule)
    cyc_of = {p: Fr(3 * ell, D) for key, (cyc, ell, D) in cycles.items() for p in cyc}
    rates[beta] = (rule, cyc_of)


def edge_rate(beta, node):
    rule, cyc_of = rates[beta]
    while node not in cyc_of:
        node = rule[node][0]
    return cyc_of[node]


def T_of(B, q):
    out = subprocess.run(["./traj", "diffT", str(B), str(q[0]), str(q[1]), str(q[2])], capture_output=True, text=True).stdout.split()
    return int(out[0])


def facet_rate(B, q):
    u, v, w = q
    beta = B % 2
    m = 2 * v - B
    S = u + w - B
    E = 2 * u - 2 * w - B
    d = u - 3 * w
    mm = m
    if d > 0:
        N = (S - 16) // 4
        s_ = S - 4 * N
        for _ in range(N):
            mm = abs(mm + 1) - 5
        return edge_rate(beta, (mm, s_)), 1
    N = (E - 12) // 4
    e = E - 4 * N
    for _ in range(N):
        mm = abs(mm + 1) - 5
    Sp = (e + mm) // 2 - 1 if mm >= -1 else (e - mm) // 2 - 2
    return edge_rate(beta, (e + 2, Sp)), -1


def in_source(B, q):
    u, v, w = q
    if (B + 0) % 1:
        return None
    m = 2 * v - B
    S = u + w - B
    E = 2 * u - 2 * w - B
    d = u - 3 * w
    J = B - u - v + 2 * w
    K1 = B - 2 * u + 2 * v - w
    L = B - v - w
    base = (B >= 100000 and u <= B - 1000 and w >= 1000 and m >= 0 and min(S - m, E - m, J, K1, L) >= 1000
            and B > u >= v >= w)
    if not base:
        return None
    if d >= 8000 and 3 * d <= B - 24000:
        return 1
    if -d >= 8000 and B + 6 * d >= 24000:
        return -1
    return None


if __name__ == "__main__":
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    nE = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    nF = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    # ---- Proposition 4 ----
    lo, hi, worst = None, None, 0
    cases = []
    for i in range(nE):
        B = rnd.randrange(100000, 400000)
        beta = B % 2
        m = rnd.choice([t for t in range(-20, 21) if (t - beta) % 2 == 0])
        S = rnd.randrange(-20, 21)
        ulo = -((-3 * B) // 4) + 1000
        uhi = (5 * B) // 6 - 1000
        if i % 10 == 0:
            u = ulo
        elif i % 10 == 1:
            u = uhi
        else:
            u = rnd.randrange(ulo, uhi + 1)
        cases.append((B, m, S, u))
    for corner in ((-20, -20), (-20, 20), (20, -20), (20, 20), (-19, -20), (19, 20)):
        for B in (100000, 100001, 100002, 100003):
            m, S = corner
            if (m - B) % 2:
                continue
            for u in (-((-3 * B) // 4) + 1000, (5 * B) // 6 - 1000):
                cases.append((B, m, S, u))
    for (B, m, S, u) in cases:
        q = (u, (B + m) // 2, B - u + S)
        h = T_of(B, q) - 1
        r = edge_rate(B % 2, (m, S))
        disc = h - r * (Fr(u) - Fr(3 * B, 4))
        lo = disc if lo is None or disc < lo else lo
        hi = disc if hi is None or disc > hi else hi
    print("Prop 4: %d cases (fresh B in [1e5,4e5) + corners/extremes); h - r(u-3B/4) in [%s, %s] = [%.3f, %.3f]; Lemma A.2 interval [-43/2, 343/8]; inside=%s"
          % (len(cases), lo, hi, float(lo), float(hi), lo >= Fr(-43, 2) and hi <= Fr(343, 8)))
    sys.stdout.flush()
    # ---- Proposition 5 ----
    done = {1: 0, -1: 0}
    maxdev = {1: 0, -1: 0}
    mingap = {1: 10**9, -1: 10**9}
    fastcount = {}
    tries = 0
    while min(done.values()) < nF // 2 and tries < 10**6:
        tries += 1
        B = rnd.randrange(100000, 300000)
        x = rnd.uniform(0.75, 1.0)
        y = rnd.uniform(0.5, 2 / 3)
        z = rnd.uniform(1 / 6, 1 / 3)
        q = (int(x * B), int(y * B), int(z * B))
        side = in_source(B, q)
        if side is None or done[side] >= nF // 2:
            continue
        r, side2 = facet_rate(B, q)
        assert side2 == side
        u, v, w = q
        if side == 1:
            tf = Fr(9 * u - 3 * w - 6 * B, 8)
            ts = Fr(5 * u - 3 * w - 3 * B, 4)
        else:
            tf = Fr(3 * (u + w - B), 4)
            ts = Fr(2 * u + 6 * w - 3 * B, 4)
        T = T_of(B, q)
        pred, other = (tf, ts) if r == Fr(3, 2) else (ts, tf)
        dev = abs(T - pred)
        maxdev[side] = max(maxdev[side], dev)
        mingap[side] = min(mingap[side], abs(T - other))
        key = (side, B % 2)
        fastcount.setdefault(key, [0, 0])
        fastcount[key][0] += (r == Fr(3, 2))
        fastcount[key][1] += 1
        done[side] += 1
    for side in (1, -1):
        print("Prop 5 side %+d: %d fresh source states; max |T - B t_r| = %.2f (claim <= 216); min distance to the other clock = %.1f"
              % (side, done[side], float(maxdev[side]), float(mingap[side])))
    for k, (a, n) in sorted(fastcount.items()):
        print("   side %+d beta=%d: fast-rate fraction %d/%d = %.3f" % (k[0], k[1], a, n, a / n))
