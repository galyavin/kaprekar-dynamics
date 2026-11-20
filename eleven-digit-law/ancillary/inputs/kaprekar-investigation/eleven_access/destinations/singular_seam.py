"""Eliminate every two-return cylinder of the one-dimensional +0- seam."""
from itertools import product
from math import gcd


X = (1, 0, 0, 0)
Z = (0, 1, 0, 0)
KAPPA = (0, 0, 1, 0)
ONE = (0, 0, 0, 1)
ZERO = (0, 0, 0, 0)


def add(*rows):
    return tuple(map(sum, zip(*rows)))


def mul(n, row):
    return tuple(n*x for x in row)


def sub(a, b):
    return add(a, mul(-1, b))


def weak_positive(row, branch):
    return (row, row) if branch else (ZERO, mul(-1, row))


def step(state, branch):
    x, y, z, kappa = state
    r, guard_r = weak_positive(z, branch[0])
    t, guard_t = weak_positive(sub(r, x), branch[1])
    u, guard_u = weak_positive(sub(sub(x, r), kappa), branch[2])
    k = add(x, mul(-1, r), t)
    out = (sub(add(k, t), u),
           add(mul(-1, y), z, k),
           add(mul(2, r), mul(-1, y), t),
           add(kappa, mul(-1, y), z, k, u, mul(2, t)))
    return out, (guard_r, guard_t, guard_u)


def region(state):
    x, y, z, kappa = state
    scale = add(y, mul(-1, x), mul(-1, z), mul(-1, kappa))
    return (scale, add(scale, x), add(mul(2, scale), z, mul(-1, kappa)))


def canonical(row):
    divisor = gcd(*row)
    return tuple(x//divisor for x in row) if divisor else row


def feasible(guards):
    """Fourier–Motzkin on rows a*x+b*z+c*kappa+d >= 0."""
    rows = set(map(canonical, guards))
    for j in range(3):
        positive = [row for row in rows if row[j] > 0]
        negative = [row for row in rows if row[j] < 0]
        next_rows = {row for row in rows if row[j] == 0}
        for upper in negative:
            for lower in positive:
                combined = add(mul(-upper[j], lower), mul(lower[j], upper))
                next_rows.add(canonical(combined))
        rows = next_rows
        if any(row[3] < 0 and not any(row[:3]) for row in rows):
            return False
    return True


def main():
    start = (X, add(X, Z, KAPPA, ONE), Z, KAPPA)
    choices = tuple(product((False, True), repeat=3))
    for first, second in product(choices, repeat=2):
        mid, first_guards = step(start, first)
        end, second_guards = step(mid, second)
        guards = (*region(start), *first_guards, *region(mid),
                  *second_guards, *region(end))
        assert not feasible(guards), (first, second)
    print("One-dimensional +0- seam: all 64 two-return cylinders infeasible at L=1.")


if __name__ == "__main__":
    main()
