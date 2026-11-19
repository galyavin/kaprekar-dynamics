"""Independent checks (exact integers): exceptional cycles of Lemma 13, the
immediate-exit example of Section 7, written from the paper's definitions.
F_B is computed by positional subtraction of the canonical word."""
import sys


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
    assert borrow == 0
    return out


def diffs(w):
    D = sorted(w, reverse=True)
    return tuple(D[i] - D[10 - i] for i in range(5))


def F(q, B):
    return diffs(kap_word(list(q) + [0] * 6, B))


def gaps(q, B):
    return (B - q[0], q[0] - q[1], q[1] - q[2], q[2] - q[3], q[3] - q[4], q[4])


def MB(B):
    return -(-(B + 1) // 7) + 64


def orbit_period(q, B, maxlen=10000):
    seen = {}
    cur = tuple(q)
    for j in range(maxlen):
        if cur in seen:
            return seen[cur], j - seen[cur]
        seen[cur] = j
        cur = F(cur, B)
    return None, None


def check_four_cycle(m):
    B = 5 * m
    rows = [((4, 3, 3, 2, 1), (1, 0, -1, 0, 2)), ((4, 3, 2, 1, 1), (0, 0, -2, 0, -1)),
            ((4, 3, 3, 2, 1), (1, 0, 0, 1, 2)), ((4, 3, 2, 1, 1), (0, 0, -2, 1, -1))]
    st = [tuple(c * m + k for c, k in zip(a, b)) for a, b in rows]
    ok = all(F(st[i], B) == st[(i + 1) % 4] for i in range(4))
    mg = [min(gaps(s, B)) for s in st]
    pos = all(s[4] > 0 for s in st)
    return ok, max(mg), pos, st[0] == (4 * m + 1, 3 * m, 3 * m - 1, 2 * m, m + 2)


def check_family(r, const, period, m):
    B = 5 * m + r
    q0 = tuple(m * c + k for c, k in zip((4, 3, 2, 1, 1), const))
    if not (B > q0[0] >= q0[1] >= q0[2] >= q0[3] >= q0[4] >= 0):
        return False, None, None, None
    pre, per = orbit_period(q0, B, maxlen=4 * period + 10)
    if pre != 0:
        return False, pre, per, None
    cyc = [q0]
    for _ in range(per - 1):
        cyc.append(F(cyc[-1], B))
    mg = max(min(gaps(s, B)) for s in cyc)
    pos = all(s[4] > 0 for s in cyc)
    return per == period, per, mg, pos


def main():
    out = []
    # four-cycle, B=5m: threshold 20 and fresh m
    ms = list(range(20, 400)) + [1001, 4999, 12347, 200003, 10**6 + 3, 10**9 + 7]
    bad = [m for m in ms if check_four_cycle(m) != (True, 2, True, True)]
    out.append(f"four-cycle B=5m: tested {len(ms)} values of m in [20,399] U fresh large; failures={bad[:10]}")
    lowm = [m for m in range(1, 20) if check_four_cycle(m)[0]]
    out.append(f"  four-cycle also closes for m in {lowm} (below stated threshold 20)")
    fam = [(1, (1, 2, 2, 4, -1), 8, 6), (1, (1, 1, 1, 3, -1), 8, 6), (1, (0, 0, -1, 1, -1), 10, 4),
           (2, (1, 1, 0, 2, -1), 8, 4)]
    for r, c, per, m0 in fam:
        ms = list(range(m0, 600)) + [1001, 4999, 12347, 200003, 10**6 + 3, 10**9 + 7]
        res = [(m, check_family(r, c, per, m)) for m in ms]
        badm = [m for m, x in res if not (x[0] and x[2] is not None and x[2] <= 5 and x[3])]
        maxmg = max(x[2] for m, x in res if x[2] is not None)
        below = [m for m in range(1, m0) if check_family(r, c, per, m)[0]]
        out.append(f"B=5m+{r} const={c} period={per} from m={m0}: tested {len(ms)}; failures={badm[:10]}; "
                   f"max over phases of min gap={maxmg}; also cycles below threshold at m={below}")
    # sigma = infinity: all phases have min gap <= 5 < M_B once B large
    out.append(f"M_B at B=100: {MB(100)} (>5, so no listed cycle state is in S_B for any B>=2)")
    # immediate-exit example
    fails = []
    tested = 0
    for t in list(range(390, 3001)) + [10**4 + 1, 123457, 10**6 + 1, 10**9 + 9]:
        B = 7 * t
        M = t + 65
        A = B - 5 * M
        assert MB(B) == M
        q = (B - A, B - A - M, B - A - 2 * M, B - A - 3 * M, M)
        assert gaps(q, B) == (A, M, M, M, M, M)
        insS = min(gaps(q, B)) >= MB(B)
        qn = F(q, B)
        gn = gaps(qn, B)
        tested += 1
        if t >= 400 and not (insS and gn == (M, A + 1, M - 1, M, M, M) and min(gn) < MB(B)):
            fails.append(t)
        if t < 400 and t in (390, 395, 399):
            out.append(f"  t={t}: start in S_B={insS}, next gaps match={gn == (M, A + 1, M - 1, M, M, M)}")
    out.append(f"immediate exit B=7t: tested {tested} values t in [390,3000] U large; failures(t>=400)={fails[:10]}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
