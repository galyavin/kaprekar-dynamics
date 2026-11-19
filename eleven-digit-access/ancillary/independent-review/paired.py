"""Independent tests of the two-cycle/edge identities (Sections 4.2, 6.2; Appendix D).
Own encodings from the paper's definitions:
  Psi   : centered paired map, eq.(centered)
  P_B   : paired rule (definition in Section 4.2)
  F_B   : positional subtraction of the canonical word (no use of eq.(digits))
Tests
  (A) Psi(6q-Bv) == 6 P_B(q) - Bv for positive-innermost q
  (B) eq.(paired-return) == Psi o Psi under the two retained guards
  (C) displayed first output of the return under the guards
  (D) the two gap-increment formulas
  (E) eq.(edge) F_B = P_B + Delta e_4 near the other phase on the open edge 3a+2b=1/6
      and at its lower endpoint; F_B = P_B near the canonical phase
  (F) actual canonical return near the edge = paired return with Delta added to w'
  (G) one-unit transfer identity z+ = Tz + clip(d1,-1,1)(0,-1,1,0,0,0)
  (H) pairing criterion of Section 6.2: all gaps >= m, 7m >= B+1 => F_B = P_B
"""
import random
import sys
from fractions import Fraction as Fr

V = (5, 4, 3, 2, 1)


def kap_word(w, B):
    D = sorted(w, reverse=True)
    out = [0] * 11
    borrow = 0
    for k in range(10, -1, -1):
        t = D[k] - D[10 - k] - borrow
        if t < 0:
            t += B
            borrow = 1
        else:
            borrow = 0
        out[k] = t
    return out


def F(q, B):
    D = sorted(kap_word(list(q) + [0] * 6, B), reverse=True)
    return tuple(D[i] - D[10 - i] for i in range(5))


def P(q, B):
    L = {i: q[i - 1] - (i == 5) for i in range(1, 6)}
    R = {i: B - q[i - 1] - (i > 1) for i in range(1, 6)}
    s = [B - 1]
    for i in range(1, 6):
        a, b = L[i], R[6 - i]
        s += [max(a, b), min(a, b)]
    return tuple(s[i] - s[10 - i] for i in range(5))


def Psi(X):
    c = [-6]
    for i in range(1, 6):
        a = X[i - 1] - 6 * (i == 5)
        b = -X[6 - i - 1] - 6 * (i < 5)
        c += [max(a, b), min(a, b)]
    return tuple(c[k - 1] - c[11 - k] for k in range(1, 6))


def center(q, B):
    return tuple(6 * q[i] - B * V[i] for i in range(5))


def gaps(q, B):
    return (B - q[0], q[0] - q[1], q[1] - q[2], q[2] - q[3], q[3] - q[4], q[4])


def ghat(X):
    return (-X[0], X[0] - X[1], X[1] - X[2], X[2] - X[3], X[3] - X[4], X[4])


def pos(t):
    return t if t > 0 else 0


def ret_coords(X):
    l = -X[0]
    p = Fr(X[1] + 2 * l, 6)
    s = Fr(X[2] + 3 * l, 6)
    w = Fr(X[3] + 4 * l, 6)
    r = Fr(X[4] - l, 6)
    return l, w, p, s, r


def X_of(l, w, p, s, r):
    return (-l, -2 * l + 6 * p, -3 * l + 6 * s, 6 * w - 4 * l, l + 6 * r)


def ret_aux(l, w, p, s, r):
    Cf = pos(-l + 2 * s + 1)
    Af = max(-1, -r)
    J = pos(-r - 1) - pos(r - 1)
    Kf = pos(r + 1)
    Hf = pos(w - J - p - 1)
    t = s - p - Cf
    Df = max(Af - 1, t)
    Ef = max(Af, t - 1) - max(-t - 1, -Af)
    Ff = min(Af, t - 1) + J + 1
    return Cf, Af, J, Kf, Hf, t, Df, Ef, Ff


def guards(l, w, p, s, r):
    Cf, Af, J, Kf, Hf, t, Df, Ef, Ff = ret_aux(l, w, p, s, r)
    return (w - l + p <= -1) and (l - 2 * w - 2 * Kf >= -1)


def ret_formula(l, w, p, s, r):
    Cf, Af, J, Kf, Hf, t, Df, Ef, Ff = ret_aux(l, w, p, s, r)
    return (l - 6 * Df, w + J - p - 4 * Df, Ef - 2 * Df, Ff - 3 * Df + Hf, -p - Kf + Df - Hf)


def first_out(l, w, p, s, r):
    Cf, Af, J, Kf, Hf, t, Df, Ef, Ff = ret_aux(l, w, p, s, r)
    return (-l + 6 * Af, -2 * l + 6 * J, 3 * l - 6 * w - 6 * Kf, 2 * l + 6 * p - 6 * w, l + 6 * (p - s + Cf))


def incr(l, w, p, s, r):
    Cf, Af, J, Kf, Hf, t, Df, Ef, Ff = ret_aux(l, w, p, s, r)
    return r + J + Kf + Hf, -Cf - pos(t - Af - 1) + Hf


def gap45(X):  # q4-q5 and q3-q4 up to the common B/6 shift, in integer-gap units
    g = ghat(X)
    return Fr(g[4], 6), Fr(g[3], 6)


def test_A(rnd, out):
    bad = 0; n = 0
    for _ in range(200000):
        B = rnd.choice([rnd.randint(7, 60), rnd.randint(61, 10**4), rnd.randint(10**5, 10**9)])
        q = sorted((rnd.randint(1, B - 1) for _ in range(5)), reverse=True)
        n += 1
        X = center(q, B)
        Y = center(P(q, B), B)
        if Psi(X) != Y:
            bad += 1
    out.append(f"(A) Psi vs P_B in centered coordinates: {n} random positive-innermost states, mismatches={bad}")


def test_BCD(rnd, out):
    nB = badB = badC = badD = 0
    fails = []
    for box in (6, 20, 60, 400):
        for _ in range(60000):
            l = rnd.randint(-box, 3 * box)
            w = rnd.randint(-3 * box, box)
            p = rnd.randint(-box, box)
            s = rnd.randint(-box, box)
            r = rnd.randint(-box, box)
            if not guards(l, w, p, s, r):
                continue
            nB += 1
            X = X_of(l, w, p, s, r)
            Y1 = Psi(X)
            Y2 = Psi(Y1)
            if tuple(ret_coords(Y2)) != ret_formula(l, w, p, s, r):
                badB += 1
                if len(fails) < 8:
                    fails.append(((l, w, p, s, r), ret_coords(Y2), ret_formula(l, w, p, s, r)))
            if Y1 != first_out(l, w, p, s, r):
                badC += 1
            a0, b0 = gap45(X)
            a2, b2 = gap45(Y2)
            if (a2 - a0, b2 - b0) != incr(l, w, p, s, r):
                badD += 1
    out.append(f"(B) paired return eq.(paired-return) vs Psi^2 under the two guards: {nB} integer states "
               f"(boxes 6..400), mismatches={badB}")
    for f in fails:
        out.append(f"     B-fail state={f[0]} Psi^2={f[1]} formula={f[2]}")
    out.append(f"(C) displayed first output vs Psi: mismatches={badC}")
    out.append(f"(D) increment formulas Delta(q4-q5), Delta(q3-q4) vs Psi^2: mismatches={badD}")
    return fails


def edge_points():
    pts = []
    for a in (Fr(0), Fr(1, 60), Fr(1, 90), Fr(1, 200), Fr(1, 31), Fr(7, 300)):
        b = (Fr(1, 6) - 3 * a) / 2
        A = Fr(1, 6) + a; D = Fr(1, 6) + b; C = 1 - 4 * A - D
        canon = (A, A, A, D, C, A)
        other = (A, A, C, A, D, A)
        pts.append((a, b, canon, other))
    return pts


def x_of_g(g):
    x = [1 - g[0]]
    for i in range(1, 5):
        x.append(x[-1] - g[i])
    return tuple(x)


def test_E_F(rnd, out):
    badE = badEc = badF = nE = nEc = nF = nDelta = 0
    nFdelta = 0
    for a, b, canon, other in edge_points():
        xc, xo = x_of_g(canon), x_of_g(other)
        for B in (10**6 + rnd.randint(0, 59) for _ in range(12)):
            for _ in range(2500):
                R0 = rnd.choice([3, 20, 200, B // 2000])
                qo = [int(B * xo[i]) + rnd.randint(-R0, R0) for i in range(5)]
                qo = tuple(sorted(qo, reverse=True))
                if qo[4] > 0 and qo[0] < B:
                    nE += 1
                    Delta = pos(qo[2] + qo[3] - B + 1)
                    if Delta:
                        nDelta += 1
                    pred = tuple(P(qo, B)[i] + (Delta if i == 3 else 0) for i in range(5))
                    if F(qo, B) != pred:
                        badE += 1
                qc = [int(B * xc[i]) + rnd.randint(-R0, R0) for i in range(5)]
                qc = tuple(sorted(qc, reverse=True))
                if qc[4] > 0 and qc[0] < B:
                    nEc += 1
                    if F(qc, B) != P(qc, B):
                        badEc += 1
                    # (F) canonical return
                    X = center(qc, B)
                    l, w, p, s, r = ret_coords(X)
                    Cf, Af, J, Kf, Hf, t, Df, Ef, Ff = ret_aux(l, w, p, s, r)
                    if guards(l, w, p, s, r) and Hf == 0:
                        nF += 1
                        mid = F(qc, B)
                        Dm = pos(mid[2] + mid[3] - B + 1)
                        if Dm:
                            nFdelta += 1
                        end = F(mid, B)
                        lf, wf, pf, sf, rf = ret_formula(l, w, p, s, r)
                        pred = (lf, wf + Dm, pf, sf, rf)
                        if tuple(ret_coords(center(end, B))) != pred:
                            badF += 1
    out.append(f"(E) other phase near open edge 3a+2b=1/6 (6 edge points incl. a=0, offsets up to B/2000, B~1e6): "
               f"{nE} states ({nDelta} with Delta>0), mismatches with eq.(edge)={badE}; canonical phase F_B=P_B: "
               f"{nEc} states, mismatches={badEc}")
    out.append(f"(F) canonical returns with guards and H_f=0: {nF} ({nFdelta} with Delta>0); "
               f"actual 2-step = formula with Delta added to w' only: mismatches={badF}")


def test_G(rnd, out):
    bad = 0; n = 0
    for _ in range(200000):
        box = rnd.choice([3, 12, 50, 1000])
        beta = rnd.randint(0, 5)
        X = tuple(rnd.randint(-box, box) * 6 + i * beta for i in range(1, 6))
        u = [Fr(0)] + [Fr(t, 6) for t in X] + [Fr(0)]
        d = [u[i] + u[6 - i] for i in range(4)]   # d_0..d_3
        z = [Fr(t, 6) for t in ghat(X)]
        Tz = [None] * 6
        for i in range(3):
            m = (abs(d[i] + 1) - abs(d[i + 1] + 1))
            Tz[2 * i] = (z[i] + z[5 - i] + m) / 2
            Tz[2 * i + 1] = (z[i] + z[5 - i] - m) / 2
        mn = min(z[5], z[0] + 1)
        Tz[0], Tz[1] = mn, z[0] + z[5] - mn
        cl = max(-1, min(1, d[1]))
        pred = [Tz[0], Tz[1] - cl, Tz[2] + cl, Tz[3], Tz[4], Tz[5]]
        zp = [Fr(t, 6) for t in ghat(Psi(X))]
        n += 1
        if zp != pred:
            bad += 1
            if bad <= 3:
                out.append(f"     G-fail X={X} z+={[str(t) for t in zp]} pred={[str(t) for t in pred]}")
    out.append(f"(G) one-unit transfer identity on {n} physical-lattice X: mismatches={bad}")


def test_H(rnd, out):
    bad = 0; n = 0
    for _ in range(100000):
        B = rnd.randint(13, 10**7)
        m = -(-(B + 1) // 7)
        if 6 * m > B:
            continue
        extra = B - 6 * m
        cuts = sorted(rnd.randint(0, extra) for _ in range(5))
        h = [cuts[0]] + [cuts[i] - cuts[i - 1] for i in range(1, 5)] + [extra - cuts[4]]
        G = [m + t for t in h]
        rnd.shuffle(G)
        q = []
        acc = B
        for i in range(5):
            acc -= G[i]
            q.append(acc)
        q = tuple(q)
        assert gaps(q, B) == tuple(G)
        n += 1
        if F(q, B) != P(q, B):
            bad += 1
    out.append(f"(H) Section 6.2 pairing criterion (min gap m, 7m>=B+1 => F_B=P_B): {n} random states at threshold m, mismatches={bad}")


def main():
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    out = []
    test_A(rnd, out)
    test_BCD(rnd, out)
    test_E_F(rnd, out)
    test_G(rnd, out)
    test_H(rnd, out)
    print("\n".join(out))


if __name__ == "__main__":
    main()
