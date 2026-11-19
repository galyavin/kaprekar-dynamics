"""Independent checks of the radius-coordinate claims (Lemma 3 proof, Appendix C), exact rationals,
with the radius map derived from the definition of K:
   r(x) = sort_desc(|2x_i - 1|),  K(x) = C(r(x)),  r(Kx) = rho(r(x)),
   rho(r) = sort_desc(r1, |r1+r2-1|, |r2+r3-1|, |r3+r4-1|, |r4+r5-1|).
 (R0) r(Kx) == rho(r(x)) and K(x) == C(r(x)) on random x
 (R1) ray: normals xi=(a+b,c,d,h-a) of rho(r*+x)-r* equal F_S(xi) (Lemma 9(1))
 (R3) triangle: affine normals N(r)=(r3-r4, r1+r4+r5-1, 2r1-r3-r4-r5-1) vanish on the triangle's
      radius section, and N(rho^3 r) = (|L|, |v|+w_+-L_+, L-|w|), L=u+w-v, near interior and the
      two open edges (Lemma 9(3))
 (R4) isolated cycle: r_A=(4,3,3,2,1)/5 has period 3 under rho; energy expansion identity at r_A
"""
import random
from fractions import Fraction as Fr

S6 = Fr(1, 6)


def K(x):
    s = sorted([Fr(1)] + list(x) + [1 - t for t in x], reverse=True)
    return tuple(s[i] - s[10 - i] for i in range(5))


def rad(x):
    return tuple(sorted((abs(2 * t - 1) for t in x), reverse=True))


def C(r):
    return (Fr(1, 2) * (1 + r[0]), Fr(1, 2) * (r[0] + r[1]), Fr(1, 2) * (r[1] + r[2]), Fr(1, 2) * (r[2] + r[3]), Fr(1, 2) * (r[3] + r[4]))


def rho(r):
    return tuple(sorted([r[0]] + [abs(r[i] + r[i + 1] - 1) for i in range(4)], reverse=True))


def g_of_x(x):
    return (1 - x[0], x[0] - x[1], x[1] - x[2], x[2] - x[3], x[3] - x[4], x[4])


def x_of_g(g):
    x = [1 - g[0]]
    for i in range(1, 5):
        x.append(x[-1] - g[i])
    return tuple(x)


def energy_g(g):
    return sum((t - S6) ** 2 for t in g)


def main():
    rnd = random.Random(3)
    out = []
    # R0
    ok = True
    for _ in range(20000):
        N = rnd.randint(2, 10**6)
        x = tuple(sorted((Fr(rnd.randint(0, N), N) for _ in range(5)), reverse=True))
        ok &= K(x) == C(rad(x)) and rad(K(x)) == rho(rad(x))
    out.append(f"(R0) K = C o r and r o K = rho o r on 20000 random rational x: {ok}")
    # R1 ray
    ok = True
    n = 0
    for _ in range(20000):
        f = Fr(rnd.randint(1, 999), 6000)
        rs = (Fr(2, 3) + 2 * f, Fr(2, 3) - 2 * f, Fr(1, 3), Fr(1, 3), 2 * f)
        eps = Fr(1, 10**7)
        x = [eps * Fr(rnd.randint(-1000, 1000), 1000) for _ in range(5)]
        if x[2] < x[3]:
            x[2], x[3] = x[3], x[2]
        r = tuple(rs[i] + x[i] for i in range(5))
        a, b, c, d, h = x
        xi = (a + b, c, d, h - a)
        y = [t - u for t, u in zip(rho(r), rs)]
        xi2 = (y[0] + y[1], y[2], y[3], y[4] - y[0])
        v, c_, d_, z = xi
        FS = (-d_ - z, max(v, -c_ - d_), min(v, -c_ - d_), -v - c_)
        ok &= xi2 == FS
        n += 1
    out.append(f"(R1) ray one-step normals == F_S at {n} random small offsets, f in (0,1/6): {ok}")
    # R3 triangle
    V1 = tuple(Fr(t, 14) for t in (2, 4, 2, 0, 3, 3)); V2 = tuple(Fr(t, 8) for t in (1, 2, 1, 1, 1, 2)); V3 = tuple(Fr(t, 8) for t in (1, 2, 1, 1, 2, 1))

    def Nf(r):
        return (r[2] - r[3], r[0] + r[3] + r[4] - 1, 2 * r[0] - r[2] - r[3] - r[4] - 1)

    def Phi(u, v, w):
        L = u + w - v
        return (abs(L), abs(v) + max(w, 0) - max(L, 0), L - abs(w))

    res = {}
    for kind in ("interior", "edge12", "edge13"):
        okv = okp = True
        cnt = 0
        for _ in range(4000):
            if kind == "interior":
                l = [Fr(rnd.randint(1, 100)) for _ in range(3)]
            elif kind == "edge12":
                l = [Fr(rnd.randint(1, 100)), Fr(rnd.randint(1, 100)), Fr(0)]
            else:
                l = [Fr(rnd.randint(1, 100)), Fr(0), Fr(rnd.randint(1, 100))]
            l = [t / sum(l) for t in l]
            g = tuple(sum(li * V[i] for li, V in zip(l, (V1, V2, V3))) for i in range(6))
            x0 = x_of_g(g)
            # find the phase whose radius section lies on the stated plane
            phases = [x0, K(x0), K(K(x0))]
            r0 = None
            for ph in phases:
                if Nf(rad(ph)) == (0, 0, 0):
                    r0 = rad(ph)
                    break
            if r0 is None:
                okv = False
                continue
            assert rho(rho(rho(r0))) == r0
            eps = Fr(1, 10**8)
            x = [eps * Fr(rnd.randint(-1000, 1000), 1000) for _ in range(5)]
            r = [r0[i] + x[i] for i in range(5)]
            r = tuple(sorted(r, reverse=True))          # keep a valid sorted radius vector
            u, v, w = Nf(r)
            r3 = rho(rho(rho(r)))
            okp &= Nf(r3) == Phi(u, v, w)
            cnt += 1
        res[kind] = (okv, okp, cnt)
    out.append(f"(R3) triangle: section on plane / return identity / samples: {res}")
    # R4 isolated
    rA = (Fr(4, 5), Fr(3, 5), Fr(3, 5), Fr(2, 5), Fr(1, 5))
    secs = [rA, rho(rA), rho(rho(rA))]
    ok = rho(secs[2]) == rA
    okE = True
    sstar = Fr(1, 30)
    for _ in range(5000):
        x = [Fr(rnd.randint(-1000, 1000), 10**5) for _ in range(5)]
        a, b, c, d, h = x
        r = tuple(rA[i] + x[i] for i in range(5))
        g = g_of_x(C(r))
        lhs = energy_g(g) - sstar
        vec = (-a, -b, a - c, b - d, c - h, d + h)
        rhs = Fr(1, 10) * (-b + c + 2 * d + h) + Fr(1, 4) * sum(t * t for t in vec)
        okE &= lhs == rhs
    out.append(f"(R4) isolated: r_A period 3 under rho: {ok}; sections {[[str(t) for t in s] for s in secs]}; "
               f"energy expansion identity at 5000 offsets: {okE}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
