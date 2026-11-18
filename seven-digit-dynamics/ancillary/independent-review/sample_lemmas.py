"""Stratified random checks of Lemmas A.4(i) and A.5 at large fresh bases (supplements the
exhaustive bases 18..300 and the unbounded solver proofs).
Strata: near the attractor (3B/4,B/2,B/4) with V just above 200; uniform states; boundary
states (w=1, v=w, u=v, u=B-1, gaps equal to B/4 +- small)."""
import random, sys
import model

rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
stats = {"one": 0, "nine": 0, "ai": 0, "fail": 0}
worst9 = -10**9


def check(q, B):
    global worst9
    u, v, w = q
    if not (B > u >= v >= w >= 1):
        return
    F = model.F
    q1 = F(q, B)
    stats["one"] += 1
    if model.V(q1, B) > model.V(q, B) + 8:
        stats["fail"] += 1
        print("ONE-STEP FAIL", B, q)
    qs = [q]
    for _ in range(9):
        qs.append(F(qs[-1], B))
    if all(x[2] >= 1 for x in qs[:3]):
        stats["ai"] += 1
        if qs[3][2] == 0 and B >= 18:
            stats["fail"] += 1
            print("A4(i) FAIL", B, q)
    if all(x[2] >= 1 for x in qs) and model.V(q, B) > 200:
        stats["nine"] += 1
        dv = model.V(qs[9], B) - model.V(q, B)
        worst9 = max(worst9, dv)
        if dv >= 0:
            stats["fail"] += 1
            print("NINE-STEP FAIL", B, q)


for i in range(N):
    B = rnd.choice([rnd.randrange(400, 5000), rnd.randrange(10**5, 10**6), rnd.randrange(10**8, 10**9)])
    kind = i % 4
    if kind == 0:  # near attractor, small potential
        r = rnd.randrange(30, 150)
        q = (3 * B // 4 + rnd.randrange(-r, r + 1), B // 2 + rnd.randrange(-r, r + 1), B // 4 + rnd.randrange(-r, r + 1))
    elif kind == 1:  # uniform
        a = sorted((rnd.randrange(1, B), rnd.randrange(1, B), rnd.randrange(1, B)), reverse=True)
        q = tuple(a)
    elif kind == 2:  # boundaries
        u = rnd.randrange(2, B)
        v = rnd.randrange(1, u + 1)
        w = rnd.randrange(1, v + 1)
        c = rnd.randrange(5)
        if c == 0:
            w = 1
        elif c == 1:
            w = v
        elif c == 2:
            v = u
            w = min(w, v)
        elif c == 3:
            u = B - 1
        else:
            v, w = u, u
        q = (u, v, w)
    else:  # gaps near B/4 (where |4g-B| changes sign)
        g = [B // 4 + rnd.randrange(-60, 61) for _ in range(3)]
        w = g[2]
        v = w + g[1]
        u = v + g[0]
        q = (u, v, w)
    check(q, B)
print("sampled: one-step %d, A4(i) premises %d, nine-step premises %d, failures %d, max V9-V0 = %d"
      % (stats["one"], stats["ai"], stats["nine"], stats["fail"], worst9))
