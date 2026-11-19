"""A guarded clock near the negative eleven-slot leading two-cycle cone.

The phase recurrence is defined on all integer triples. Its interpretation
as actual subtraction requires the two-step regime guards on every return.
No global entrance, phase absorption, or first-capture theorem is assumed.
"""
from odd_base_covariance import centered, centered_step, in_region

V = (5, 4, 3, 2, 1)
M = ((-1, 0, -2), (-2, 0, -3), (2, -1, 1))
C = (3, -1, -2)


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def phase_step(u):
    """Return new phase and (a,d) displacement divided by six."""
    p, q, r = u
    A = max(-1, -r)
    J = max(0, -r-1)-max(r-1, 0)
    K = max(0, r+1)
    D = max(A-1, q-p)
    E = max(A, q-p-1)-max(p-q-1, -A)
    F = min(A, q-p-1)+J+1
    return (E-2*D, F-3*D, -p-K+D), (-D, J-p)


def phase_prefix(u, h):
    a = d = 0
    for _ in range(h):
        u, delta = phase_step(u)
        a += delta[0]
        d += delta[1]
    return u, (a, d)


def second_shift(u, h):
    """Change of actual q_2, when all h two-step returns are guarded."""
    v, (a, _) = phase_prefix(u, h)
    return -2*a+v[0]-u[0]


def offsets(a, d, u):
    p, q, r = u
    return (-a, -2*a+6*p, -3*a+6*q, d, a+6*r)


def coordinates(B, q):
    if len(q) != 5:
        raise ValueError('Eleven slots require five differences')
    x = centered(B, q)
    nums = (x[1]-2*x[0], x[2]-3*x[0], x[4]+x[0])
    if any(t % 6 for t in nums):
        raise ValueError('Nonphysical centered phase')
    return -x[0], x[3], tuple(t//6 for t in nums)


def pair(B, a, d, u):
    nums = tuple(B*v+x for v, x in zip(V, offsets(a, d, u)))
    if any(t % 6 for t in nums):
        raise ValueError('Base and centered coordinates have incompatible cosets')
    return tuple(t//6 for t in nums)


def in_regime(B, q):
    """Both actual paired-level guards, negative then positive inner signs."""
    if len(q) != 5 or not in_region(B, q):
        return False
    x = centered(B, q)
    if x[1]+x[3] > -6 or x[2] > -3:
        return False
    y = centered_step(x)
    nums = tuple(B*v+t for v, t in zip(V, y))
    if any(t % 6 for t in nums):
        return False
    q1 = tuple(t//6 for t in nums)
    return (in_region(B, q1) and y[1]+y[3] >= -6 and y[2] >= -3)


def two_step(B, q):
    if not in_regime(B, q):
        raise ValueError('Two-step actual sorting regime fails')
    a, d, u = coordinates(B, q)
    v, (da, dd) = phase_step(u)
    return pair(B, a+6*da, d+6*dd, v)


def affine_phase(u):
    return tuple(sum(a*x for a, x in zip(row, u))+c for row, c in zip(M, C))


def positive_cell(u):
    p, q, r = u
    return min(r-1, p-q-2, p+r+2, p+2*r-1) >= 0


def in_phase_polytope(u):
    for _ in range(6):
        if not positive_cell(u):
            return False
        u = affine_phase(u)
    return True


def long_family(t, j=0):
    if not isinstance(t, int) or not isinstance(j, int) or t < 2 or not 0 <= j <= t:
        raise ValueError('Require integer t>=2 and 0<=j<=t')
    return 600*t, (497*t+3*j, 394*t-4+6*j, 291*t-10+9*j,
                   182*t+9*j, 103*t+5-3*j)


def long_family_base(B, j=0):
    """A guarded prefix of length 8 floor(B/600), in every B>=1200."""
    if not isinstance(B, int) or B < 1200 or not isinstance(j, int) or not 0 <= j <= B//600:
        raise ValueError('Require B>=1200 and integer 0<=j<=floor(B/600)')
    t = B//600
    a = 18*t+(5*B)%6-18*j
    d = -108*t+(-2*B)%6+54*j
    return pair(B, a, d, (-4, -10, 5))
