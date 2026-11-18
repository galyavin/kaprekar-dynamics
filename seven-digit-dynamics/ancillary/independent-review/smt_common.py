"""Own SMT encoding of the integer difference map F_B (written from the paper's definitions).

Sorting networks differ from the authors': a 7-input odd-even transposition network
(21 comparators) for the general step, and a 12-comparator 6-input network
(Dobbelaere's list) for the positive-state step with the maximal digit B-1 removed.
Both are validated by the 0-1 principle and by all permutations below.
"""
from itertools import permutations, product
import z3

# 7-input odd-even transposition sort: 7 rounds alternating (0,1),(2,3),(4,5) / (1,2),(3,4),(5,6)
OET7 = []
for r in range(7):
    start = r % 2
    for i in range(start, 6, 2):
        OET7.append((i, i + 1))

# 6-input optimal network (12 comparators), different from the authors' listing.
NET6 = [(0, 5), (1, 3), (2, 4), (1, 2), (3, 4), (0, 3), (2, 5), (0, 1), (2, 3), (4, 5), (1, 2), (3, 4)]


def apply_net_py(net, xs):
    xs = list(xs)
    for i, j in net:
        if xs[i] < xs[j]:
            xs[i], xs[j] = xs[j], xs[i]
    return xs


def validate_networks():
    for n, net in ((7, OET7), (6, NET6)):
        for bits in product((0, 1), repeat=n):
            assert apply_net_py(net, bits) == sorted(bits, reverse=True)
        for p in permutations(range(n)):
            assert apply_net_py(net, p) == sorted(p, reverse=True)
    return True


def zmax(a, b):
    return z3.If(a >= b, a, b)


def zmin(a, b):
    return z3.If(a >= b, b, a)


def apply_net(net, xs):
    xs = list(xs)
    for i, j in net:
        hi, lo = zmax(xs[i], xs[j]), zmin(xs[i], xs[j])
        xs[i], xs[j] = hi, lo
    return xs


def diffs7(s):
    return (s[0] - s[6], s[1] - s[5], s[2] - s[4])


def step_pos(q, B, net6=True):
    """Image of a positive state (w>=1) via the hand-derived borrow digits."""
    u, v, w = q
    if net6:
        # B-1 is the maximum digit: drop it, sort the six others
        s = apply_net(NET6, [B - u, w - 1, B - v - 1, u, B - w - 1, v])
        return (B - 1 - s[5], s[0] - s[4], s[1] - s[3])
    digs = [B - 1, B - u, w - 1, v, B - w - 1, u, B - v - 1]
    return diffs7(apply_net(OET7, digs))


def step_general(q, B):
    """Image of any valid state, case split on w>0, v>0 (7-input network for each case)."""
    u, v, w = q
    pos = diffs7(apply_net(OET7, [u, v, w - 1, B - 1, B - w - 1, B - v - 1, B - u]))
    zer = diffs7(apply_net(OET7, [u, v - 1, B - 1, B - 1, B - 1, B - v - 1, B - u]))
    zz = diffs7(apply_net(OET7, [u - 1, B - 1, B - 1, B - 1, B - 1, B - 1, B - u]))
    return tuple(z3.If(w > 0, pos[i], z3.If(v > 0, zer[i], zz[i])) for i in range(3))


def valid(q, B):
    u, v, w = q
    return z3.And(B > u, u >= v, v >= w, w >= 0, u > 0)


def zabs(x):
    return z3.If(x >= 0, x, -x)


def Vz(q, B):
    u, v, w = q
    return zabs(4 * (B - u) - B) + zabs(4 * (u - v) - B) + zabs(4 * (v - w) - B) + zabs(4 * w - B)


if __name__ == "__main__":
    print("networks valid:", validate_networks(), len(OET7), len(NET6))
