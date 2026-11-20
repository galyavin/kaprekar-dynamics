"""Exact cylinder bounds for the positive-anchor edge of Two."""
from fractions import Fraction
from itertools import combinations, product
from math import gcd


FRONT_DEPTH = 3
LOOKAHEAD = 4
IDENTITY = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
BOUNDS = tuple(map(Fraction, ("3/5", "2/3", "3/8", "1/2",
                              "5/13", "2/3", "3/8", "1/2")))
BRANCH_MATRICES = (
    ((0, 0, 1, 0), (1, 0, -1, 0), (0, -1, -1, 1), (0, 1, 1, 0)),
    ((1, 0, 0, 0), (-1, 0, 1, 0), (0, -1, -1, 1), (0, 1, 1, 0)),
)


def step(state, branch):
    return tuple(tuple(sum(weight * row[j] for weight, row in zip(weights, state))
                       for j in range(4))
                 for weights in BRANCH_MATRICES[branch])


def triple(u, v, w):
    return dot(u, cross(v, w))


def dot(u, v):
    return sum(x*y for x, y in zip(u, v))


def cross(u, v):
    a, b, c = u
    d, e, f = v
    return (b*f-c*e, c*d-a*f, a*e-b*d)


def canonical(row):
    divisor = gcd(*row)
    return tuple(x//divisor for x in row) if divisor else row


def vertices(guards):
    """Enumerate vertices of the bounded polytope on the section W=1."""
    rows = sorted({canonical(guard) for guard in guards})
    if any(not any(row[:3]) and row[3] < 0 for row in rows):
        return
    rows = tuple((row[:3], row[3]) for row in rows if any(row[:3]))
    for (u, i), (v, j), (w, k) in combinations(rows, 3):
        denominator = triple(u, v, w)
        if denominator == 0:
            continue
        numerator = tuple(-i*a-j*b-k*c
                          for a, b, c in zip(cross(v, w), cross(w, u), cross(u, v)))
        if denominator < 0:
            denominator = -denominator
            numerator = tuple(-value for value in numerator)
        if all(dot(row, numerator)+offset*denominator >= 0 for row, offset in rows):
            yield (*numerator, denominator)


def cylinder_maximum(word):
    state = IDENTITY
    guards = list(state) + [(-1, 0, 0, 2)]  # a <= 2 W for an infinite survivor
    for t, branch in enumerate(word, start=1):
        state = step(state, branch)
        # The nonnegativity of the new |a-v| row encodes the branch guard.
        guards.extend(state)
        if t == FRONT_DEPTH:
            output_scale = state[-1]
    maximum = max((Fraction(dot(output_scale, vertex), vertex[-1])
                   for vertex in vertices(guards)), default=Fraction())
    return maximum


def main():
    from verify import determinant

    assert all(abs(determinant(m)) == 1 for m in BRANCH_MATRICES)
    for prefix, expected in zip(product((0, 1), repeat=FRONT_DEPTH), BOUNDS):
        maximum = max(cylinder_maximum(prefix + suffix)
                      for suffix in product((0, 1), repeat=LOOKAHEAD))
        assert maximum == expected, (prefix, maximum, expected)
        print("Edge cylinder", "".join(map(str, prefix)), "scale <=", maximum)
    rate = sum(bound**4 for bound in BOUNDS)
    assert rate == Fraction(2105699171233, 2961204480000) < 1
    print(f"Positive-anchor edge: verified all 128 exact cylinders, volume factor {rate}.")


if __name__ == "__main__":
    main()
