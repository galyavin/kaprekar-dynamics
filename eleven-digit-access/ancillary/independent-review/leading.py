"""Independent checks of the leading map K in exact rationals (Sections 2-4, App. A):
 - the five all-positive cells reaching e>=s_* (ray, triangle+2 edges, isolated orbit)
 - triangle/edge energy infima and supremum (5/384, 1/48, 1/48; 1/21, unattained)
 - the two-cycle family eq.(two), its phase swap and the leading adjacent slacks
 - classification of lattice periodic points produced by lattice.c
 - Lemma 2 majorization and Lemma 3 entry on random rational orbits
Usage: python3 leading.py [lattice_output_file ...]
"""
import random
import sys
from fractions import Fraction as Fr

S_STAR = Fr(1, 30)
SIXTH = Fr(1, 6)


def K(x):
    s = sorted([Fr(1)] + list(x) + [1 - t for t in x], reverse=True)
    return tuple(s[i] - s[10 - i] for i in range(5))


def g_of_x(x):
    return (1 - x[0], x[0] - x[1], x[1] - x[2], x[2] - x[3], x[3] - x[4], x[4])


def x_of_g(g):
    x1 = 1 - g[0]
    x2 = x1 - g[1]
    x3 = x2 - g[2]
    x4 = x3 - g[3]
    x5 = x4 - g[4]
    assert x5 == g[5]
    return (x1, x2, x3, x4, x5)


def energy_g(g):
    return sum((t - SIXTH) ** 2 for t in g)


def energy(x):
    return energy_g(g_of_x(x))


def Kn(x, n):
    for _ in range(n):
        x = K(x)
    return x


def gv(*a, d):
    return tuple(Fr(t, d) for t in a)


V1 = gv(2, 4, 2, 0, 3, 3, d=14)
V2 = gv(1, 2, 1, 1, 1, 2, d=8)
V3 = gv(1, 2, 1, 1, 2, 1, d=8)
ISO = gv(1, 1, 2, 3, 1, 2, d=10)


def comb(lams, vs):
    return tuple(sum(l * v[i] for l, v in zip(lams, vs)) for i in range(6))


def report(lines, s):
    lines.append(s)
    print(s, flush=True)


def check_cells(lines, rnd):
    # ray: gaps 1/6 + f(-1,1,1,-1,-1,1)
    ok = True
    for f in [Fr(1, 13), Fr(1, 12), Fr(1, 10), Fr(1, 7), Fr(1, 6) - Fr(1, 1000), Fr(3, 40), Fr(1, 20), Fr(1, 100)]:
        g = tuple(SIXTH + f * c for c in (-1, 1, 1, -1, -1, 1))
        x = x_of_g(g)
        ok &= K(x) == x and energy_g(g) == 6 * f * f
    report(lines, f"ray: fixed by K at 8 rational f (incl. f<1/sqrt180); e=6f^2: {ok}; e>=s_* iff f>=1/sqrt(180)")
    # triangle: period three on vertices, interior, edges; K affine on the closed triangle
    xs = [x_of_g(v) for v in (V1, V2, V3)]
    per3 = all(Kn(x, 3) == x and K(x) != x for x in xs)
    aff = True
    per3_pts = True
    for _ in range(300):
        l = [Fr(rnd.randint(0, 50)) for _ in range(3)]
        if sum(l) == 0:
            continue
        l = [t / sum(l) for t in l]
        g = comb(l, (V1, V2, V3))
        x = x_of_g(g)
        for k in (1, 2, 3):
            imgs = [Kn(v, k) for v in xs]
            lhs = Kn(x, k)
            rhs = tuple(sum(li * im[i] for li, im in zip(l, imgs)) for i in range(5))
            aff &= lhs == rhs
        per3_pts &= Kn(x, 3) == x
    report(lines, f"triangle: vertices period 3: {per3}; K, K^2, K^3 affine on closed triangle at 300 random rational points: {aff}; "
                  f"K^3=id there: {per3_pts}")
    report(lines, "  triangle phase vertex images (gaps): " + "; ".join(
        str([tuple(str(t) for t in g_of_x(Kn(x, k))) for x in xs]) for k in (1, 2)))
    # energy ranges: closed-triangle min and where attained; vertex energies
    eV = [energy_g(v) for v in (V1, V2, V3)]

    def seg_min(P, Q):  # min of convex quadratic on segment [P,Q]
        d = [Q[i] - P[i] for i in range(6)]
        a = sum(t * t for t in d)
        b = 2 * sum((P[i] - SIXTH) * d[i] for i in range(6))
        t = -b / (2 * a)
        t = min(max(t, Fr(0)), Fr(1))
        pt = tuple(P[i] + t * d[i] for i in range(6))
        return energy_g(pt), t

    m12, t12 = seg_min(V1, V2)
    m13, t13 = seg_min(V1, V3)
    m23, t23 = seg_min(V2, V3)
    # interior stationary point of e on the triangle's affine plane
    # minimize |P + s u + t w - c|^2
    c = (SIXTH,) * 6
    u = [V2[i] - V1[i] for i in range(6)]
    w = [V3[i] - V1[i] for i in range(6)]
    r0 = [V1[i] - c[i] for i in range(6)]
    uu = sum(a * a for a in u); ww = sum(a * a for a in w); uw = sum(a * b for a, b in zip(u, w))
    ur = sum(a * b for a, b in zip(u, r0)); wr = sum(a * b for a, b in zip(w, r0))
    det = uu * ww - uw * uw
    s_ = (-ur * ww + wr * uw) / det
    t_ = (-wr * uu + ur * uw) / det
    inside = s_ > 0 and t_ > 0 and s_ + t_ < 1
    report(lines, f"triangle energies: e(V1)={eV[0]} e(V2)={eV[1]} e(V3)={eV[2]}; min on edge V1V2={m12} at t={t12}; "
                  f"min on edge V1V3={m13} at t={t13}; min on edge V2V3={m23} at t={t23}; plane stationary point inside: {inside} (s={s_}, t={t_})")
    report(lines, f"  => interior inf={min(m12, m13, m23) if not inside else 'interior min'} (unattained if on excluded edge V2V3), "
                  f"edge infima {m12},{m13} (at V2/V3 endpoints, excluded), common sup {max(eV)} (at V1, excluded)")
    # isolated orbit
    xi = x_of_g(ISO)
    report(lines, f"isolated orbit: period3={Kn(xi, 3) == xi and K(xi) != xi}, e={energy_g(ISO)} (=s_* : {energy_g(ISO) == S_STAR}); "
                  f"phases={[tuple(str(t) for t in g_of_x(Kn(xi, k))) for k in range(3)]}")
    # sharp work orbit
    xs_ = gv(12, 8, 6, 6, 3, d=14)
    orb = [Kn(xs_, k) for k in range(4)]
    report(lines, f"sharp orbit (12,8,6,6,3)/14 -> {[tuple(int(t * 14) for t in o) for o in orb]}; period 3: {orb[3] == orb[0]}; "
                  f"energies {[str(energy(o)) for o in orb[:3]]}")
    return (V1, V2, V3), xs


def check_two_cycle(lines, rnd):
    ok_swap = ok_slack = True
    nonper_outside = 0
    outside_ball_tried = 0
    tried = 0
    for _ in range(1500):
        a = Fr(rnd.randint(0, 1200), 36000)
        b = Fr(rnd.randint(0, 3000), 36000)
        if b < a:
            a, b = b, a
        A = SIXTH + a; D = SIXTH + b; C = 1 - 4 * A - D
        p1 = (A, A, A, D, C, A); p2 = (A, A, C, A, D, A)
        if min(C, A, D) <= 0:
            continue
        inside = 3 * a + 2 * b <= SIXTH
        x1, x2 = x_of_g(p1), x_of_g(p2)
        if inside:
            tried += 1
            ok_swap &= K(x1) == x2 and K(x2) == x1
            # leading levels (x_i, 1 - x_{6-i}); slack = min(level i) - max(level i+1)
            for x, exp in ((x1, (C, 1 - 5 * A, 1 - 5 * A, C)), (x2, (A, 1 - 3 * A - 2 * D, 1 - 3 * A - 2 * D, A))):
                lev = [(x[i], 1 - x[4 - i]) for i in range(5)]
                sl = tuple(min(lev[i]) - max(lev[i + 1]) for i in range(4))
                ok_slack &= sl == exp
        elif energy_g(p1) < S_STAR:
            outside_ball_tried += 1
            if not (K(K(x1)) == x1):
                nonper_outside += 1
    report(lines, f"two-cycle eq.(two): {tried} rational (a,b) with b>=a>=0, 3a+2b<=1/6: K swaps phases: {ok_swap}; "
                  f"leading adjacent slacks as displayed: {ok_slack}; in-ball samples with 3a+2b>1/6: {outside_ball_tried}, of which not 2-periodic: {nonper_outside}")


def in_tri(g, vs):
    # barycentric coordinates of g in conv(vs) (6-dim gaps, sum 1); exact
    P, Q, R = vs
    u = [Q[i] - P[i] for i in range(6)]
    w = [R[i] - P[i] for i in range(6)]
    r = [g[i] - P[i] for i in range(6)]
    # solve using normal equations, then verify exactly
    uu = sum(a * a for a in u); ww = sum(a * a for a in w); uw = sum(a * b for a, b in zip(u, w))
    ur = sum(a * b for a, b in zip(u, r)); wr = sum(a * b for a, b in zip(w, r))
    det = uu * ww - uw * uw
    s = (ur * ww - wr * uw) / det
    t = (wr * uu - ur * uw) / det
    if any(P[i] + s * u[i] + t * w[i] != g[i] for i in range(6)):
        return None
    return (1 - s - t, s, t)


def classify(files, lines, tri_vs, tri_xs):
    tri_phases = []
    for k in range(3):
        tri_phases.append(tuple(g_of_x(Kn(x, k)) for x in tri_xs))
    xi = x_of_g(ISO)
    iso_ph = {g_of_x(Kn(xi, k)) for k in range(3)}
    cntI = cntO = 0
    unexplained = []
    kinds = {}
    for fn in files:
        for line in open(fn):
            parts = line.split()
            if not parts or parts[0] not in ("I", "O"):
                continue
            N = int(parts[1])
            q = [Fr(int(t), N) for t in parts[2:7]]
            g = g_of_x(q)
            if parts[0] == "I":
                cntI += 1
                kind = None
                if g[0] == g[3] == g[4] and g[1] == g[2] == g[5] and g[1] > g[0]:
                    kind = "ray"
                elif g in iso_ph:
                    kind = "isolated"
                else:
                    for k, ph in enumerate(tri_phases):
                        lam = in_tri(g, ph)
                        if lam is None:
                            continue
                        l1, l2, l3 = lam
                        if min(lam) < 0:
                            continue
                        if l1 > 0 and l2 > 0 and l3 > 0:
                            kind = f"triangle-int"
                        elif l1 > 0 and l2 > 0 and l3 == 0:
                            kind = "triangle-edge12"
                        elif l1 > 0 and l3 > 0 and l2 == 0:
                            kind = "triangle-edge13"
                        if kind:
                            break
                if kind is None:
                    unexplained.append(("I", N, tuple(str(t) for t in g), str(energy_g(g))))
                else:
                    kinds[kind] = kinds.get(kind, 0) + 1
            else:
                cntO += 1
                a = g[0] - SIXTH
                ok = False
                if all(t == SIXTH for t in g):
                    ok = True
                for ph in ((0, 1, 2, 3, 4, 5),):
                    pass
                A = g[0]
                # canonical (A,A,A,D,C,A) or other (A,A,C,A,D,A)
                for pat in ((0, 1, 2, 5), (0, 1, 3, 5)):
                    if all(g[i] == A for i in pat):
                        rest = [i for i in range(6) if i not in pat]
                        if pat == (0, 1, 2, 5):
                            D, C = g[3], g[4]
                        else:
                            C, D = g[2], g[4]
                        b = D - SIXTH
                        if b >= a >= 0 and 3 * a + 2 * b <= SIXTH and C == 1 - 4 * A - D:
                            ok = True
                if not ok:
                    unexplained.append(("O", N, tuple(str(t) for t in g), str(energy_g(g))))
    report(lines, f"lattice classification: tagI(interior, e>=s_*)={cntI} by kind {kinds}; tagO(g0>=1/6, e<s_*)={cntO}; "
                  f"unexplained={len(unexplained)}")
    for u in unexplained[:30]:
        report(lines, f"  UNEXPLAINED {u}")


def random_orbits(lines, rnd, n=4000, Nden=1000003):
    """Lemma 2 (majorization, monotone energy) and Lemma 3 (entry) on random rational orbits."""
    fails = []
    maj_bad = 0
    steps_hist = {}
    for _ in range(n):
        q = sorted((rnd.randint(0, Nden) for _ in range(5)), reverse=True)
        x = tuple(Fr(t, Nden) for t in q)
        seen = {}
        j = 0
        while True:
            e = energy(x)
            if e < S_STAR:
                steps_hist[j] = steps_hist.get(j, 0) + 1
                break
            if x in seen:
                fails.append((tuple(str(t) for t in g_of_x(x)), str(e)))
                break
            seen[x] = j
            y = K(x)
            gs = sorted(g_of_x(x), reverse=True); gk = sorted(g_of_x(y), reverse=True)
            ps = pk = 0
            for a_, b_ in zip(gs, gk):
                ps += a_; pk += b_
                if pk > ps:
                    maj_bad += 1
                    break
            x = y
            j += 1
    report(lines, f"random rational orbits (denominator {Nden}, n={n}): majorization violations={maj_bad}; "
                  f"orbits never entering e<s_*: {len(fails)}; entry-step histogram={dict(sorted(steps_hist.items()))}")
    for f in fails[:10]:
        report(lines, f"  NON-ENTERING {f}")


def main():
    rnd = random.Random(20261118)
    lines = []
    tri_vs, tri_xs = check_cells(lines, rnd)
    check_two_cycle(lines, rnd)
    if len(sys.argv) > 1:
        classify(sys.argv[1:], lines, tri_vs, tri_xs)
    random_orbits(lines, rnd)


if __name__ == "__main__":
    main()
