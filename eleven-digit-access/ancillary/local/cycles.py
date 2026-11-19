"""Exact verification of the bounded cycles at the upper Two endpoint, for each base class mod 5.

For B = 5m + r the cycle states are written q = m*P + e with P one of the two
leading phases (4,3,3,2,1), (4,3,2,1,1) of the upper Two endpoint and e a fixed integer residual.
Two independent checks:
  1. Full-integer subtraction on a realizing word for many m (big integers, no update formula).
  2. A symbolic affine proof for all m >= m0: the borrowed digit rows (coefficient of m,
     constant) sort in the same order for every m >= m0, and their symmetric differences
     reproduce the next affine state.  This proves the cycle for every B = 5m + r, m >= m0.
Usage: python3 -B verify_cycles.py
"""
from fractions import Fraction

PH = ((4, 3, 3, 2, 1), (4, 3, 2, 1, 1))


def full_step(B, q):
    word = (*q, 0, 0, 0, 0, 0, 0)          # sorted digits with differences q
    hi = lo = 0
    for a, b in zip(word, reversed(word)):
        hi = hi * B + a
        lo = lo * B + b
    n = hi - lo
    dig = []
    for _ in range(11):
        n, r = divmod(n, B)
        dig.append(r)
    assert n == 0
    dig.sort(reverse=True)
    return tuple(dig[i] - dig[10 - i] for i in range(5))


def gaps(B, q):
    return (B - q[0], q[0] - q[1], q[1] - q[2], q[2] - q[3], q[3] - q[4], q[4])


def orbit(B, q):
    seen = [q]
    while True:
        q = full_step(B, q)
        if q == seen[0]:
            return seen
        assert q not in seen and len(seen) < 100
        seen.append(q)


def affine(m, q):
    best = min(PH, key=lambda P: max(abs(x - m * p) for x, p in zip(q, P)))
    return best, tuple(x - m * p for x, p in zip(q, best))


def symbolic(r, cyc):
    """cyc: list of (P, e).  Returns the least m0 for which the affine sort is valid for all m>=m0."""
    m0 = 1
    for j, (P, e) in enumerate(cyc):
        # physical validity: B-1 >= q1 >= ... >= q5 >= 1 with B = 5m + r
        lin = [(5 - P[0], r - 1 - e[0])] + [(P[i] - P[i + 1], e[i] - e[i + 1]) for i in range(4)] \
              + [(P[4], e[4] - 1)]
        for a, b in lin:
            assert a > 0 or b >= 0
            if a > 0 and b < 0:
                m0 = max(m0, -(b // a))
        rows = [(5, r - 1)] + [(P[i], e[i] - (i == 4)) for i in range(5)] \
               + [(5 - P[i], r - e[i] - (i != 0)) for i in range(5)]
        rows.sort(reverse=True)               # lexicographic: coefficient of m, then constant
        for (a, b), (c, d) in zip(rows, rows[1:]):
            if a > c:                          # need (a-c) m + (b-d) >= 0
                m0 = max(m0, -((b - d) // (a - c)) if b < d else 1)
            else:
                assert a == c and b >= d
        nxt = tuple((rows[i][0] - rows[10 - i][0], rows[i][1] - rows[10 - i][1]) for i in range(5))
        P2, e2 = cyc[(j + 1) % len(cyc)]
        assert nxt == tuple(zip(P2, e2)), (j, nxt, P2, e2)
    return m0


def main():
    reps = {0: [(0, 0, -2, 0, -1)], 1: [(0, 1, 2, 4, -1), (0, 0, 1, 3, -1), (-1, -1, -1, 1, -1)],
            2: [(-1, 0, -1, 2, -1)]}
    for r, offs in reps.items():
        for off in offs:
            m = 200000
            B = 5 * m + r
            # localscan reports offsets relative to round(B p/5) at phase (4,3,2,1,1)
            base = tuple((2 * B * p + 5) // 10 for p in PH[1])
            q = tuple(b + o for b, o in zip(base, off))
            cyc_q = orbit(B, q)
            cyc = [affine(m, x) for x in cyc_q]
            m0 = symbolic(r, cyc)
            # full-integer replay at many m, including small ones above m0
            for mm in list(range(max(m0, 2), max(m0, 2) + 40)) + [10**3 + 7, 10**9 + 3, 10**15 + 1]:
                BB = 5 * mm + r
                qs = [tuple(mm * p + x for p, x in zip(P, e)) for P, e in cyc]
                for a, b in zip(qs, qs[1:] + qs[:1]):
                    assert full_step(BB, a) == b
                    g = gaps(BB, a)
                    assert a[4] > 0 and min(g) <= 5
            mins = [min(gaps(B, x)) for x in cyc_q]
            per = len(cyc)
            print(f"B=5m+{r}: period {per}, min gaps along cycle {mins}, positive innermost at every state,"
                  f" valid for all m >= {m0} (symbolic), full-integer replay OK")
            for P, e in cyc:
                print(f"    m*{P} + {e}")


if __name__ == '__main__':
    main()
