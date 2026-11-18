"""Small independent checks of hand-proved steps.
(1) Theorem 4 identity F_{B+jm}(q+jv)=F_B(q)+jv on the paired region, k=1..6, and
    that the translated source is paired; also how often it fails off the region.
(2) number map q -> K-output is injective (all ranks), n=5,7,9, small B.
(3) T_B(N) = 1 + h(q(N)) for transient q(N), and h<=T<=h+1 always: direct
    iteration of K_B on all 9-digit strings at B=4,5.
(4) multiplicity: #{N: q(N)=q} vs leading term B^{k+1} rho_n(q/B), exact counts n=5,B=14 and n=7,B=8.
"""
import random, itertools, math
from collections import Counter
from km import *

random.seed(4)


def paired(q, B):
    k = len(q)
    if q[-1] <= 0:
        return False
    a = [q[i] - (1 if i == k - 1 else 0) for i in range(k)]
    b = [B - q[k - 1 - i] - (1 if i < k - 1 else 0) for i in range(k)]
    return all(min(a[i], b[i]) >= max(a[i + 1], b[i + 1]) for i in range(k - 1))


ok = True
for k in range(1, 7):
    n = 2 * k + 1
    m = k + 1
    v = tuple(range(k, 0, -1))
    good = bad = tried = offbad = offtried = 0
    while good < 3000:
        B = random.choice([random.randint(5, 60), random.randint(100, 10 ** 6)])
        q = tuple(sorted(random.sample(range(1, B), k) if B > k else [1] * k, reverse=True))
        j = random.randint(0, 50)
        lhs = F(tuple(t + j * s for t, s in zip(q, v)), B + j * m, n)
        rhs = tuple(t + j * s for t, s in zip(F(q, B, n), v))
        if paired(q, B):
            good += 1
            if lhs != rhs or not paired(tuple(t + j * s for t, s in zip(q, v)), B + j * m):
                bad += 1
        elif j > 0:
            offtried += 1
            offbad += lhs != rhs
    print(f"(1) k={k}: paired sources {good}, covariance/pairing failures {bad};  unpaired sources {offtried}, identity fails on {offbad}")
    ok &= bad == 0


def numdigits(q, B, n):
    return tuple(kaprekar_digits(realize(q, n), B))


for n, Bs in ((5, range(2, 16)), (7, range(2, 10)), (9, range(2, 8))):
    k = n // 2
    coll = 0
    for B in Bs:
        seen = {}
        for q in itertools.combinations_with_replacement(range(B - 1, -1, -1), k):
            d = numdigits(q, B, n)
            if d in seen:
                coll += 1
            seen[d] = q
    print(f"(2) n={n}: number map injective on all q for B in {list(Bs)[0]}..{list(Bs)[-1]}: collisions={coll}")
    ok &= coll == 0

for B in (4, 5):
    n = 9
    Kmemo = {}

    def K(N):
        if N not in Kmemo:
            Kmemo[N] = tuple(kaprekar_digits(list(N), B))
        return Kmemo[N]

    # cyclic numbers: iterate from every string; collect cycles
    cyc = set()
    allN = list(itertools.product(range(B), repeat=n))
    status = {}
    for N in allN:
        path = []
        x = N
        while x not in status and x not in path:
            path.append(x)
            x = K(x)
        if x in path:
            i = path.index(x)
            for y in path[i:]:
                cyc.add(y)
                status[y] = 0
            path = path[:i]
        for y in reversed(path):
            status[y] = status[K(y)] + 1 if y not in cyc else 0
    # difference-level h
    hq = {}

    def h(q):
        if q in hq:
            return hq[q]
        seen = []
        x = q
        while x not in seen:
            seen.append(x)
            x = F(x, B, n)
        i = seen.index(x)
        cycq = set(seen[i:])
        x = q
        t = 0
        while x not in cycq:
            x = F(x, B, n)
            t += 1
        hq[q] = t
        return t

    viol = viol2 = 0
    for N in allN:
        if len(set(N)) == 1:
            continue
        q = qvec(N)
        T = status[N]
        hh = h(q)
        if not (hh <= T <= hh + 1):
            viol += 1
        if hh > 0 and T != hh + 1:
            viol2 += 1
    print(f"(3) B={B}: all {B**n} strings: h<=T<=h+1 violations={viol}; T=h+1 for transient q violations={viol2}")
    ok &= viol == 0 and viol2 == 0

for n, B in ((5, 14), (7, 8)):
    k = n // 2
    cnt = Counter(qvec(N) for N in itertools.product(range(B), repeat=n))
    maxrel = 0
    rows = []
    for q, c in cnt.items():
        g = gaps(q, B)
        if min(g) >= 2:
            lead = math.factorial(n) * math.prod(g)
            # exact count for strict q with all-distinct digit strings: (B-q1)*prod(g_i-1)*n!
            strict = math.factorial(n) * g[0] * math.prod(t - 1 for t in g[1:])
            rows.append((c, lead, strict))
    # count >= strict part, and |count-lead| = O(B^k)
    worst = max(abs(c - l) / B ** k / math.factorial(n) for c, l, s in rows)
    below = sum(1 for c, l, s in rows if c < s)
    print(f"(4) n={n},B={B}: {len(rows)} strict q; count>=distinct-digit term always: {below==0}; "
          f"max |count - n! prod g| / (n! B^k) = {worst:.3f} (bounded => O(B^k) remainder)")
print("MISC OK" if ok else "MISC FAILED")
