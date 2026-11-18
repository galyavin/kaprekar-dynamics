"""Independent transition model for odd-digit Kaprekar dynamics.

Written from the paper's definitions only (Section 1 and eq. (borrow) is NOT
assumed: the base model is column subtraction of the ascending arrangement
from the descending one).  All helpers here are independent of the authors'
ancillary code.
"""
from fractions import Fraction as Fr
import itertools, fcntl, os, contextlib

LOCK = os.environ.get("KAPREKAR_COMPUTE_LOCK", os.path.join(os.path.dirname(__file__), "compute.lock"))


@contextlib.contextmanager
def compute_lock():
    f = open(LOCK, "a")
    fcntl.flock(f, fcntl.LOCK_EX)
    try:
        yield
    finally:
        fcntl.flock(f, fcntl.LOCK_UN)
        f.close()


# ---------------------------------------------------------------- integer model
def kaprekar_digits(digits, B):
    """K_B on a digit string (most significant first): desc - asc by
    column subtraction with borrows.  Returns the digit list."""
    n = len(digits)
    desc = sorted(digits, reverse=True)
    asc = desc[::-1]
    out = [0] * n
    borrow = 0
    for j in range(n - 1, -1, -1):
        t = desc[j] - asc[j] - borrow
        if t < 0:
            t += B
            borrow = 1
        else:
            borrow = 0
        out[j] = t
    assert borrow == 0
    return out


def qvec(digits):
    d = sorted(digits, reverse=True)
    n = len(d)
    k = n // 2
    return tuple(d[i] - d[n - 1 - i] for i in range(k))


def realize(q, n):
    """A digit string with difference vector q (sorted descending)."""
    k = n // 2
    return list(q) + [q[-1] if k else 0] + [0] * k


def F(q, B, n=9):
    return qvec(kaprekar_digits(realize(q, n), B))


def gaps(q, B):
    k = len(q)
    return [B - q[0]] + [q[i] - q[i + 1] for i in range(k - 1)] + [q[-1]]


def S(g, j):
    return sum(sorted(g, reverse=True)[:j])


def G(q, B):
    g = gaps(q, B)
    return sum(abs(g[i] - g[j]) for i in range(len(g)) for j in range(i + 1, len(g)))


# ---------------------------------------------------------------- real map
def f_real(x):
    """Continuous leading map: sort 1, x_i, 1-x_i and take symmetric
    differences."""
    vals = [Fr(1) if isinstance(x[0], Fr) else 1.0] + list(x) + [1 - t for t in x]
    d = sorted(vals, reverse=True)
    n = len(d)
    return tuple(d[i] - d[n - 1 - i] for i in range(len(x)))


def g1(x):
    one = Fr(1) if isinstance(x[0], Fr) else 1.0
    return [one - x[0]] + [x[i] - x[i + 1] for i in range(len(x) - 1)] + [x[-1]]


def majorized(h, g, tol=0):
    """h  <  g  (equal sums assumed): S_j(h) <= S_j(g) for all j."""
    hs = sorted(h, reverse=True)
    gs = sorted(g, reverse=True)
    a = b = 0
    for j in range(len(hs)):
        a += hs[j]
        b += gs[j]
        if a > b + tol:
            return False
    return True


# ---------------------------------------------------------------- nine-digit real objects
CORE = [
    (-2, 0, 2, 1, 1), (-1, -1, -1, -1, 2), (-1, -1, 1, 1, 1), (0, -1, -1, 0, 1),
    (0, 0, 0, 1, 0), (0, 0, 1, -1, 0), (0, 1, -1, 0, 0), (1, -1, 0, 0, 0),
    (1, 0, 1, 1, -1), (1, 2, 1, 0, -2), (2, 0, -1, 0, -1), (2, 1, -2, -2, -1),
]


def Cmap(r):
    a, b, c, d = r
    h = Fr(1, 2)
    return (h * (1 + a), h * (a + b), h * (b + c), h * (c + d))


def Cinv(y):
    u, v, w, z = y
    a = 2 * u - 1
    b = 2 * v - a
    c = 2 * w - b
    d = 2 * z - c
    return (a, b, c, d)


def Hmap(r):
    a, b, c, d = r
    return tuple(sorted([a, abs(a + b - 1), abs(b + c - 1), abs(c + d - 1)], reverse=True))


def Dvec(r):
    a, b, c, d = r
    return (1 - a, 1 - b, a - c, b - d, c + d)


def core_slacks(r):
    return [row[0] * r[0] + row[1] * r[1] + row[2] * r[2] + row[3] * r[3] + row[4] for row in CORE]


def in_core(r, strict=False):
    s = core_slacks(r)
    return all(t > 0 for t in s) if strict else all(t >= 0 for t in s)


def Rmap(y):
    u, v, w, z = y
    return (u, u - z, v - z, v - w)


def Q(s):
    return (s, 2 * s - 1, 3 * s - 2, 1 - s)


# ---------------------------------------------------------------- sorting networks
def bubble_network(n):
    net = []
    for i in range(n):
        for j in range(n - 1 - i):
            net.append((j, j + 1))
    return net


def check_network(net, n):
    """0-1 principle: network sorts descending all 2^n binary inputs."""
    for bits in itertools.product((0, 1), repeat=n):
        v = list(bits)
        for i, j in net:
            if v[i] < v[j]:
                v[i], v[j] = v[j], v[i]
        if v != sorted(bits, reverse=True):
            return False
    return True


def batcher_network(n):
    """Batcher odd-even mergesort comparators for n (any n, by padding)."""
    p = 1
    while p < n:
        p *= 2
    net = []

    def merge(lo, m, r):
        step = r * 2
        if step < m:
            merge(lo, m, step)
            merge(lo + r, m, step)
            for i in range(lo + r, lo + m - r, step):
                net.append((i, i + r))
        else:
            net.append((lo, lo + r))

    def sort(lo, m):
        if m > 1:
            h = m // 2
            sort(lo, h)
            sort(lo + h, h)
            merge(lo, m, 1)

    sort(0, p)
    return [(i, j) for (i, j) in net if i < n and j < n]
