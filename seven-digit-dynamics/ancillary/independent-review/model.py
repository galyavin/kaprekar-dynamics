"""Reference model written from the paper's definitions only (Section 1-2).

Two independent routes:
  * word level: sort seven base-B digits, subtract integers, re-expand (kappa_B);
  * difference level: q=(u,v,w) -> digits of L_B(q) by integer arithmetic -> sort -> differences.
The borrow formula (eq. borrow) is re-derived by hand and checked against both.
"""
from itertools import product
import random


def to_digits(n, B, k=7):
    out = []
    for _ in range(k):
        n, r = divmod(n, B)
        out.append(r)
    assert n == 0
    return out[::-1]  # high to low


def from_digits(ds, B):
    n = 0
    for d in ds:
        n = n * B + d
    return n


def kappa(word, B):
    """One Kaprekar step on a 7-digit word (tuple of digits high->low)."""
    hi = sorted(word, reverse=True)
    lo = sorted(word)
    return tuple(to_digits(from_digits(hi, B) - from_digits(lo, B), B))


def diffs(word):
    d = sorted(word, reverse=True)
    return (d[0] - d[6], d[1] - d[5], d[2] - d[4])


def L(q, B):
    u, v, w = q
    return (B**6 - 1) * u + (B**5 - B) * v + (B**4 - B**2) * w


def F_direct(q, B):
    return diffs(to_digits(L(q, B), B))


def borrow_digits(q, B):
    """Hand-derived column subtraction (high -> low)."""
    u, v, w = q
    if w > 0:
        return (u, v, w - 1, B - 1, B - w - 1, B - v - 1, B - u)
    if v > 0:
        return (u, v - 1, B - 1, B - 1, B - 1, B - v - 1, B - u)
    return (u - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - u)


def F(q, B):
    return diffs(borrow_digits(q, B))


def states(B):
    for u in range(1, B):
        for v in range(0, u + 1):
            for w in range(0, v + 1):
                yield (u, v, w)


def V(q, B):
    u, v, w = q
    g = (B - u, u - v, v - w, w)
    return sum(abs(4 * x - B) for x in g)


def orbit_h(q, B):
    """First cyclic index of q under F_B (difference level)."""
    seen = {}
    x = q
    i = 0
    while x not in seen:
        seen[x] = i
        x = F(x, B)
        i += 1
    return seen[x]  # index of first cyclic element


if __name__ == "__main__":
    # 1. borrow formula vs integer subtraction, exhaustive small bases, random large bases
    for B in range(2, 41):
        for q in states(B):
            assert tuple(to_digits(L(q, B), B)) == borrow_digits(q, B), (B, q)
    rnd = random.Random(1)
    for _ in range(20000):
        B = rnd.randrange(2, 10**7)
        u = rnd.randrange(1, B)
        v = rnd.randrange(0, u + 1)
        w = rnd.randrange(0, v + 1)
        q = (u, v, w)
        assert tuple(to_digits(L(q, B), B)) == borrow_digits(q, B), (B, q)
    # 2. word-level step agrees with L_B(diffs(word))
    for _ in range(20000):
        B = rnd.randrange(2, 10**6)
        word = tuple(rnd.randrange(B) for _ in range(7))
        if len(set(word)) == 1:
            continue
        assert from_digits(kappa(word, B), B) == L(diffs(word), B)
        assert diffs(kappa(word, B)) == F(diffs(word), B)
    # 3. h <= T <= h+1 relation on all words of small bases (word-level brute force)
    for B in (3, 4, 5):
        words = [w for w in product(range(B), repeat=7) if len(set(w)) > 1]
        succ = {w: kappa(w, B) for w in words}
        # cyclic words: iterate succ |words| times from each? use standard: w cyclic iff returns
        cyc = set()
        for w in words:
            x = w
            path = []
            seen = set()
            while x not in seen:
                seen.add(x)
                path.append(x)
                x = succ[x]
            # x is first repeated -> cycle from x
            y = x
            while True:
                cyc.add(y)
                y = succ[y]
                if y == x:
                    break
        for w in words:
            T = 0
            x = w
            while x not in cyc:
                x = succ[x]
                T += 1
            h = orbit_h(diffs(w), B)
            assert h <= T <= h + 1
            if w not in cyc:
                assert T == h + 1, (B, w, T, h)
    print("model self-checks passed")
