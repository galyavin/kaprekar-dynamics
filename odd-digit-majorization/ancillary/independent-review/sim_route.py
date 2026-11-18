"""Fresh-base test of the Section 4.4 table (route A0,B1,AB) by direct digit
subtraction.  Inputs: the 16 inverse branches x of the radius (18,14,6,3)/25
(x_i=(1+-r_i)/2), integer q = round(Bx) + uniform offset in [-R,R]^4.
Per sample the prediction uses that sample's own y=f(q/B) (real leading map):
times L+c(s-4/5), c in {7/2,9/2,50/9,37/6}; terminal Q((s+4)/6) for c=7/2 and
Q(4/5) otherwise.  Observed T = 1+mu (mu = first cyclic difference time, the
input string is non-cyclic since q is transient) and Z = successor of the
first cyclic difference.

usage: sim_route.py NSAMP RFRAC B1 [B2 ...]
"""
import sys, random, subprocess, itertools, json, math
from fractions import Fraction as Fr
from concurrent.futures import ThreadPoolExecutor
from km import f_real

HERE = __file__.rsplit("/", 1)[0]
r = [Fr(18, 25), Fr(14, 25), Fr(6, 25), Fr(3, 25)]
BR = []
for sg in itertools.product((1, -1), repeat=4):
    BR.append(sorted([(1 + s * ri) / 2 for s, ri in zip(sg, r)], reverse=True))
CS = [Fr(7, 2), Fr(9, 2), Fr(50, 9), Fr(37, 6)]
TABLE = {0: [Fr(11, 48), Fr(5, 48), Fr(137, 288), Fr(55, 288)],
         1: [Fr(29, 128), Fr(5, 48), Fr(47, 96), Fr(23, 128)],
         2: [Fr(29, 144), Fr(5, 48), Fr(607, 1152), Fr(193, 1152)],
         3: [Fr(107, 576), Fr(5, 48), Fr(69, 128), Fr(197, 1152)],
         4: [Fr(233, 1152), Fr(5, 48), Fr(73, 144), Fr(215, 1152)]}
for b, row in TABLE.items():
    assert sum(row) == 1


def Q(s):
    return (s, 2 * s - 1, 3 * s - 2, 1 - s)


def gen(B, n, R, seed):
    rnd = random.Random(seed)
    out = []
    while len(out) < n:
        bi = rnd.randrange(16)
        x = BR[bi]
        q = [int(round(float(xi) * B)) + rnd.randint(-R, R) for xi in x]
        if not (B - 1 >= q[0] > q[1] > q[2] > q[3] > 0):
            continue
        out.append((bi, q))
    return out


def run_chunk(B, items):
    inp = "".join(f"{B} {q[0]} {q[1]} {q[2]} {q[3]}\n" for _, q in items)
    res = subprocess.run([HERE + "/kap9", "orbits"], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    return [tuple(map(int, l.split())) for l in res if l.strip()]


def analyse(B, items, outs):
    beta = B % 5
    cnt = [0] * 4
    unmatched = 0
    maxdt = 0.0
    maxdz = 0.0
    mind2 = 1.0
    notpaired = 0
    perbranch = {}
    for (bi, q), o in zip(items, outs):
        mu, lam, z = o[0], o[1], o[2:6]
        allpos, allpaired = o[6], o[7]
        if not allpaired:
            notpaired += 1
        T = mu + 1
        x = [Fr(t, B) for t in q]
        y = f_real(x)
        u, v, w, zz = y
        L = Fr(7, 4) * u - v - zz / 4 - Fr(3, 4)
        s = -u / 4 + v / 2 - zz / 4 + Fr(3, 4)
        ts = [float(L + c * (s - Fr(4, 5))) for c in CS]
        tobs = T / B
        d = sorted((abs(tobs - t), i) for i, t in enumerate(ts))
        i = d[0][1]
        zpred = Q((s + 4) / 6) if i == 0 else Q(Fr(4, 5))
        dz = max(abs(z[j] / B - float(zpred[j])) for j in range(4))
        sep = min(abs(ts[a] - ts[b]) for a in range(4) for b in range(a + 1, 4))
        if d[0][0] > 0.25 * sep or dz > 1e-3:
            unmatched += 1
            continue
        maxdt = max(maxdt, d[0][0])
        maxdz = max(maxdz, dz)
        mind2 = min(mind2, sep)
        cnt[i] += 1
        perbranch.setdefault(bi, [0] * 4)[i] += 1
    n = sum(cnt)
    freq = [c / n for c in cnt] if n else [0] * 4
    pred = [float(p) for p in TABLE[beta]]
    se = [math.sqrt(p * (1 - p) / n) if n else 0 for p in pred]
    z = [(f - p) / e if e else 0 for f, p, e in zip(freq, pred, se)]
    return dict(B=B, beta=beta, n=n, unmatched=unmatched, notpaired=notpaired, counts=cnt, freq=freq, pred=pred,
                zscore=z, max_time_err=maxdt, max_terminal_err=maxdz, min_atom_sep=mind2, perbranch=perbranch)


def main():
    n = int(sys.argv[1])
    rfrac = float(sys.argv[2])
    bases = [int(b) for b in sys.argv[3:]]
    results = []
    for B in bases:
        R = max(1, int(rfrac * B))
        items = gen(B, n, R, seed=B)
        k = 8
        chunks = [items[i::k] for i in range(k)]
        with ThreadPoolExecutor(k) as ex:
            outs_c = list(ex.map(lambda c: run_chunk(B, c), chunks))
        items2, outs = [], []
        for c, o in zip(chunks, outs_c):
            items2 += c
            outs += o
        res = analyse(B, items2, outs)
        res["R"] = R
        results.append(res)
        print(f"B={B} beta={res['beta']} R={R} n={res['n']} unmatched={res['unmatched']} notpaired={res['notpaired']} "
              f"freq={[round(f,4) for f in res['freq']]} pred={[round(p,4) for p in res['pred']]} "
              f"z={[round(t,2) for t in res['zscore']]} maxdt={res['max_time_err']:.2e} maxdz={res['max_terminal_err']:.2e}",
              flush=True)
    json.dump(results, open(HERE + f"/sim_route_{bases[0]}.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
