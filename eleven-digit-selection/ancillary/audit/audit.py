"""Check the even exclusion leaf dump independently of its search and lattice code.

Run with Python 3.9+ in a private snapshot:
    python3 -B audit.py leaves.jsonl.gz /path/to/compute.lock
The digest of a gzip file can depend on its modification time. The report
identifies the copy reviewed; the mathematical checks use the decoded data.
"""
import argparse
import fcntl
import gzip
import json
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import factorial, lcm, prod


CENTER = (26293, 27557, 28481, 40447, 29865, 27357)
EPSILON = Q(1, 10**6)
CUTOFF = 271137
RAMP = Q(131, 5000)
EXPECTED = {
    Q(9451, 337500): (50204, Q(-16099, 33)),
    Q(295607, 10080000): (99024, Q(-36471, 56)),
    Q(10371, 350000): (90772, Q(-76877, 154)),
    Q(179987, 6000000): (121560, Q(-2645, 4)),
}


def basis_form(rows):
    """Canonical row Hermite basis of a full-rank integral row lattice."""
    a = [list(row) for row in rows]
    dimension = len(a[0])
    for col in range(dimension):
        pivot = next(i for i in range(col, len(a)) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        for i in range(col + 1, len(a)):
            while a[i][col]:
                q = a[col][col] // a[i][col]
                reduced = [x - q*y for x, y in zip(a[col], a[i])]
                a[col], a[i] = a[i], reduced
        if a[col][col] < 0:
            a[col] = [-v for v in a[col]]
        for i in range(col):
            q = a[i][col] // a[col][col]
            a[i] = [x - q*y for x, y in zip(a[i], a[col])]
    assert all(not any(row) for row in a[dimension:])
    return tuple(tuple(row) for row in a[:dimension])


def representative(p, h):
    """Reduce by the five triangular pivots, from the first coordinate on."""
    a, b, c, d, e = p
    i, a = divmod(a, h[0][0])
    b, c, d, e = b-i*h[0][1], c-i*h[0][2], d-i*h[0][3], e-i*h[0][4]
    i, b = divmod(b, h[1][1])
    c, d, e = c-i*h[1][2], d-i*h[1][3], e-i*h[1][4]
    i, c = divmod(c, h[2][2])
    d, e = d-i*h[2][3], e-i*h[2][4]
    i, d = divmod(d, h[3][3])
    return a, b, c, d, (e-i*h[3][4]) % h[4][4]


def divided(v):
    assert all(type(n) is int and n % 6 == 0 for n in v)
    return tuple(n // 6 for n in v)


@lru_cache(None)
def leaf_basis(raw):
    return basis_form(tuple(divided(row) for row in raw))


def centered(y):
    return tuple((i+1) * sum(y) - 6*sum(y[:i+1]) for i in range(5))


def dot(a, b):
    return sum(i*j for i, j in zip(a, b))


@lru_cache(None)
def corner_values(row):
    cs = [Q(c, sum(CENTER)) for c in CENTER[:5]]
    values = []
    for signs in product((-1, 1), repeat=5):
        y = tuple(c+s*EPSILON for c, s in zip(cs, signs))
        values.append(dot(row, centered((*y, 1-sum(y)))))
    return min(values), max(values)


@lru_cache(None)
def integer_clock(row, offset):
    scale = lcm(*(v.denominator for v in (*row, offset)))
    return scale, tuple(int(scale*r) for r in row), int(scale*offset)


@lru_cache(None)
def integral_increments(row, h):
    # The centered state X equals 6p in the partition coordinates.
    return all((6*dot(row, v)).denominator == 1 for v in h)


def centered_step(x):
    # Build offsets in the five permanently paired rank levels.
    entries = [-6]
    for i, a in enumerate(x):
        u, v = a - 6*(i == 4), -x[4-i] - 6*(i < 4)
        entries.extend((max(u, v), min(u, v)))
    return tuple(entries[i] - entries[-1-i] for i in range(5))


def whole_step(x, base):
    # G1 gives the next integer word directly from the five differences.
    q = [divmod(a + (5-i)*base, 6) for i, a in enumerate(x)]
    assert all(r == 0 for _, r in q)
    digits = []
    n = sum((base**(10-i) - base**i)*v for i, (v, _) in enumerate(q))
    for _ in range(11):
        n, d = divmod(n, base)
        digits.append(d)
    assert n == 0
    digits.sort()
    return tuple(6*(digits[-1-i]-digits[i]) - (5-i)*base for i in range(5))


def residual(x):
    return max(x[0], -x[-1], *(b-a for a, b in zip(x, x[1:])))


@lru_cache(None)
def trace(x, base):
    seen = {}
    while x not in seen:
        seen[x] = len(seen)
        assert len(seen) <= 2000
        y = centered_step(x)
        assert whole_step(x, base) == y
        x = y
    return seen[x], len(seen) - seen[x]


def check_odd():
    base = 180000
    seed = (*CENTER[:4], CENTER[4]+1, CENTER[5])
    shift = tuple(a-b for a, b in zip(centered(seed), centered(CENTER)))
    assert shift == (1, 2, 3, 4, -1)
    assert base % 3840 == 3360 and 3840 % 6 == 0
    # Translating every gap by k leaves X fixed on both arms.
    for k in (1, 317, 639):
        for g in (CENTER, seed):
            assert centered(tuple(n+k for n in g)) == centered(g)
    odd = tuple(map(Q, ('205/1152','715/1152','-65/192','-715/1152','1235/1152')))
    assert corner_values(odd) == (Q(2219399,86400000), Q(2222801,86400000))
    assert corner_values(odd)[1] < Q(13,500)
    # The parity reconstruction record, lines 183--199, in paired stage order.
    # Each entry is (five gap coefficients, base coefficient, moduli for
    # even/odd). Compute the order of the *image* in the product of cyclic
    # groups, independently of the successive domain-kernel calculation.
    congruences = (
        ((-66,-60,-24,-42,-96), 48, (36,36)),
        ((-78,-90,-24,-54,-150), 66, (72,72)),
        ((-48,-54,-36,-36,-114), 48, (72,144)),
        ((132,141,48,90,237), -108, (36,72)),
        ((-198,-195,-30,-120,-285), 138, (360,288)),
        ((-66,-65,-10,-40,-95), 46, (480,768)),
    )
    for row, b, _ in congruences:
        assert sum(row)+6*b == 0
    rows = [row for row,_,mods in congruences for _ in mods]
    moduli = [m for _,_,mods in congruences for m in mods]
    diagonal = [[m*int(i==j) for j in range(len(moduli))]
                for i,m in enumerate(moduli)]
    h = basis_form((*diagonal,*zip(*rows)))
    index = Q(prod(moduli),prod(h[i][i] for i in range(len(moduli))))
    assert index == 3317760
    # Independent tensor Simpson integral on the five free gap coordinates.
    a = EPSILON/2
    integral = 0
    for node in product((-1, 0, 1), repeat=5):
        y = [Q(n,base)+a*j for n,j in zip(CENTER,node)]
        integral += prod(y)*(1-sum(y)) * 4**node.count(0)
    m = integral*(a/3)**5*Q(factorial(11))/index
    assert m > 0
    print(f'odd mass {m} = {float(m):.13g}')
    d = (-96,-192,-192,-192,-96)
    spot = tuple(x-30*v for x,v in zip((31,50,75,82,47),d))
    assert trace(spot,271141) == (619,24)


def check_dump(path):
    lattices = defaultdict(set)
    counts, mins = Counter(), {}
    tails, endpoints, histories = set(), set(), set()
    with gzip.open(path, 'rt') as source:
        for line in source:
            entry = json.loads(line)
            point = divided(entry['point'])
            h = leaf_basis(tuple(map(tuple, entry['basis'])))
            rep = representative(point, h)
            assert rep not in lattices[h], 'overlap for equal lattices'
            lattices[h].add(rep)

            row = tuple(Q(v) for v in entry['cut'][:-1])
            offset = Q(entry['cut'][-1])
            scale, coeff, c = integer_clock(row, offset)
            assert (6*dot(coeff, point) + c) % scale == 0
            assert integral_increments(row, h)
            lo = corner_values(row)[0]
            counts[lo] += 1
            mins[lo] = min(mins.get(lo, offset), offset)
            history = tuple((n,tuple(d)) for n,d in entry['hist'])
            histories.add(history)
            length, delta = history[-1]
            assert (length,delta) == (entry['word_len'],tuple(entry['delta']))
            endpoint = tuple(entry['endpoint'])
            endpoints.add((length,delta,endpoint))
            spot = tuple(x-30*d for x,d in zip(endpoint, delta))
            assert spot == tuple(entry['spot']) and residual(spot)>504
            tails.add(spot)

    assert set(counts) == set(EXPECTED)
    assert sum(counts.values()) == 361560 and len(histories) == 325
    for lo, (count, offset) in EXPECTED.items():
        assert (counts[lo],mins[lo]) == (count,offset)
        assert lo + offset/CUTOFF > RAMP
    first = min(EXPECTED)
    for lo, (_,offset) in EXPECTED.items():
        assert lo*CUTOFF + offset >= first*CUTOFF + mins[first]

    density = sum(Q(len(reps), prod(h[i][i] for i in range(5)))
                  for h, reps in lattices.items())
    assert density == 1
    pairs = 0
    for left, right in combinations(lattices,2):
        # The two affine cosets intersect iff their offsets agree modulo L+M.
        common = basis_form((*left,*right))
        l, r = sorted((lattices[left],lattices[right]),key=len)
        images = {representative(p,common) for p in l}
        assert all(representative(p,common) not in images for p in r), 'cross-lattice overlap'
        pairs += 1
    # Every missing class in the finite common quotient has positive density.
    print(f'partition {sum(counts.values())} leaves, {len(lattices)} lattices, {pairs} pairs, density {density}')
    print(f'form counts {sorted(counts.items())}')
    indices = [prod(h[i][i] for i in range(5)) for h in lattices]
    print(f'indices relative to physical lattice {min(indices)}..{max(indices)}')
    clocks = [trace(spot,271140) for spot in tails]
    assert min(mu for mu,_ in clocks) > 0
    print(f'{len(histories)} routes, {len(endpoints)} terminal family-endpoint pairs, {len(tails)} preterminal states')
    print(f'preterminal residuals {min(map(residual,tails))}..{max(map(residual,tails))}')
    print(f'preterminal difference clocks {min(n for n,_ in clocks)}..{max(n for n,_ in clocks)}, periods {sorted({p for _,p in clocks})}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dump')
    parser.add_argument('lock')
    args = parser.parse_args()
    with open(args.lock, 'r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        check_odd()
        check_dump(args.dump)
