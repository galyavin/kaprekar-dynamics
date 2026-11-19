"""Cross-check the lane summaries recorded in the six-residue verification record against
proof.md's table and the clock formula, using mixture.py's transcription.

Checks per beta: (1) lane set equals the table's; (2) recorded clock form
equals P + kappa A; (3) recorded arrival form equals P + A/(2 v_b) (neutral:
recorded as 0); (4) sum of count/index per lane equals the table mass, and
per-beta total is 1 (index histogram = the leaves' 1/I weights);
(5) leaf count equals the stated count; (6) offsets inside (-144, 144);
(7) my sampled offsets (if a samples file list is given) lie inside the
recorded per-lane offset interval.
Usage: tree_check.py tree.txt [sample files...]
"""
import sys
from fractions import Fraction as F

import mixture as M

STATED_LEAVES = {0: 361560, 1: 393589, 2: 851056, 3: 328355, 4: 1035718, 5: 318584}


def form_from(nums):
    """(a1..a5, const, den) -> affine form in y via x"""
    *a, const, den = nums
    return M.lin(*[(F(a[i], den), M.x_form(i + 1)) for i in range(5)], (F(const, den), (F(1), [F(0)] * 5)))


def same(f, g):
    return f[0] == g[0] and all(u == v for u, v in zip(f[1], g[1]))


def main():
    path = sys.argv[1]
    lanes = {}
    leaves = {}
    beta = None
    for line in open(path):
        t = line.split()
        if not t:
            continue
        if t[0] == 'beta':
            beta = int(t[1])
            leaves[beta] = int(t[t.index('leaves') + 1])
            lanes[beta] = []
        elif t[0] == 'lane':
            d, na, nb, wd, wn = map(int, t[1:6])
            va, vb = F(na, d), F(nb, d)
            w = None if wd == 0 else F(wn, wd)
            i = t.index('arrival')
            arr = list(map(int, t[i + 1:i + 8]))
            i = t.index('clock')
            clk = list(map(int, t[i + 1:i + 8]))
            i = t.index('offsets')
            omin, omax = F(t[i + 1]), F(t[i + 2])
            i = t.index('indices')
            hist = [tuple(map(int, s.split(':'))) for s in t[i + 1:]]
            lanes[beta].append(((va, vb, w), arr, clk, omin, omax, hist))
    allok = True
    rec = {}
    for b in range(6):
        table = {lane: (k, p) for lane, k, p in M.lanes(b)}
        got = {}
        total = F(0)
        nleaves = 0
        for lane, arr, clk, omin, omax, hist in lanes[b]:
            mass = sum(F(c, I) for I, c in hist)
            nleaves += sum(c for I, c in hist)
            total += mass
            k, p = table.get(lane, (None, None))
            ok_lane = k is not None
            ok_clock = ok_lane and same(form_from(clk), M.clock_form(k))
            va, vb, w = lane
            if va == 2 * vb:
                ok_arr = arr[:6] == [0] * 6
            else:
                ok_arr = same(form_from(arr), M.lin((F(1), M.P_FORM), (1 / (2 * vb), M.A_FORM)))
            ok_mass = ok_lane and mass == p
            ok_off = -144 < omin and omax < 144
            got[lane] = True
            rec[(b, lane)] = (omin, omax)
            ok = ok_lane and ok_clock and ok_arr and ok_mass and ok_off
            allok &= ok
            print(f"beta {b} lane {tuple(str(v) for v in lane)}: in table {ok_lane}, clock=P+kA {ok_clock}, "
                  f"arrival {ok_arr}, sum count/I = {mass} matches {ok_mass}, offsets [{float(omin):.3f},{float(omax):.3f}] "
                  f"in (-144,144) {ok_off}")
        miss = set(table) - set(got)
        ok = total == 1 and not miss and nleaves == leaves[b] == STATED_LEAVES[b]
        allok &= ok
        print(f"beta {b}: total weight {total}, table lanes missing from tree {miss or 'none'}, "
              f"histogram leaves {nleaves} = recorded {leaves[b]} = stated {STATED_LEAVES[b]}: {ok}")
    print(f"grand total leaves {sum(leaves.values())}")
    if len(sys.argv) > 2:
        off = M.samples(sys.argv[2:]) if False else None
        rows = M.parse(sys.argv[2:])
        outside = 0
        n = 0
        for r in rows:
            B, G, T = r[0], r[1:7], r[7]
            res = M.attribute(B, G, T)
            hits = [x for x in res if x[0] < M.OFFSET_BOUND]
            if len(hits) != 1:
                continue
            _, o, lane, _ = hits[0]
            lo, hi = rec[(B % 6, lane)]
            n += 1
            if not (lo <= o <= hi):
                outside += 1
                print(f"OFFSET OUTSIDE RECORDED RANGE: B={B} G={G} lane={lane} offset={o}")
        print(f"sampled offsets checked against recorded lane ranges: {n}, outside: {outside}")
        allok &= outside == 0
    print("ALL TREE SUMMARY CHECKS PASS" if allok else "SOME TREE CHECK FAILED")


if __name__ == "__main__":
    main()
