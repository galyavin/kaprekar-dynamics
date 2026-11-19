"""Minimize the ten-step ray return on each feasible branch and unit-box face."""

from fractions import Fraction as F
from time import monotonic

from ray_lp import lp_min


Z = F(0)


def unit(n, i):
    return [F(1) if j == i else Z for j in range(n)]


def lin(*coef):
    return tuple(F(c) for c in coef)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def face_lp(cons, outs, i, s, nv=4):
    """Minimize t subject to the branch guards, x_i=s, |x_j|<=1 and |out_k|<=t."""
    A, b = [], []
    for f in cons:
        A.append(list(f) + [Z])
        b.append(Z)
    for j in range(nv):
        row = [Z] * (nv + 1)
        row[j] = F(1)
        A.extend((row, neg(row)))
        b.extend((F(1), F(1)))
    for out in outs:
        A.extend((list(out) + [F(-1)], list(neg(out)) + [F(-1)]))
        b.extend((Z, Z))
    return lp_min([Z] * nv + [F(1)], A, b, [unit(nv + 1, i)], [F(s)])


def feasible(cons, nv=4):
    A = list(cons)
    b = [Z] * len(cons)
    for j in range(nv):
        row = unit(nv, j)
        A.extend((row, neg(row)))
        b.extend((F(1), F(1)))
    for i in range(nv):
        for s in (1, -1):
            status, *_ = lp_min([Z] * nv, A, b, [unit(nv, i)], [F(s)])
            if status != "infeasible":
                return True
    return False


def step_branches(state, cons):
    v, c, d, z = state
    a, b = v, neg(add(c, d))
    first, last = neg(add(d, z)), neg(add(v, c))
    for high, low in ((a, b), (b, a)):
        yield [first, high, low, last], cons + [sub(low, high)]


def explore(n):
    start = [lin(1, 0, 0, 0), lin(0, 1, 0, 0), lin(0, 0, 1, 0), lin(0, 0, 0, 1)]
    frontier = [(start, [lin(0, -1, 1, 0)])]  # d - c <= 0
    for _ in range(n):
        next_frontier = []
        for state, cons in frontier:
            for next_state, next_cons in step_branches(state, cons):
                if feasible(next_cons):
                    next_frontier.append((next_state, next_cons))
        frontier = next_frontier
    return frontier


def main():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    start_time = monotonic()
    frontier = explore(10)
    assert len(frontier) == 79, len(frontier)
    print(f"Ten-step feasible leaves: {len(frontier)}.", flush=True)
    best, point = None, None
    for state, cons in frontier:
        for i in range(4):
            for s in (1, -1):
                status, *solution = face_lp(cons, state, i, s)
                if status == "infeasible":
                    continue
                assert status == "opt", status
                value, coordinates = solution
                if best is None or value < best:
                    best, point = value, coordinates[:4]
    assert best == F(295, 77), best
    assert point[1] >= point[2] and max(map(abs, point)) == 1, point
    image = point
    for _ in range(10):
        v, c, d, z = image
        image = (-d - z, max(v, -c - d), min(v, -c - d), -v - c)
    assert max(map(abs, image)) == best, image
    print("Exact minimum:", best)
    print("Attaining input:", *point)
    print("Ten-step image:", *image)
    print(f"Elapsed: {monotonic() - start_time:.1f} seconds.", flush=True)


if __name__ == "__main__":
    main()
