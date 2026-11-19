"""Project the seed's guarded remainders to the gap and base lattices."""
from math import gcd, lcm
from affine_routes import witness


def extended_gcd(a, b):
    x, y, u, v = 1, 0, 0, 1
    while b:
        q, r = divmod(a, b)
        a, b, x, y, u, v = b, r, u, v, x - q * u, y - q * v
    if a < 0:
        a, x, y = -a, -x, -y
    return a, x, y


def impose(basis, row, modulus):
    """Restrict column generators to row times the increment = 0 mod modulus."""
    values = [sum(x * y for x, y in zip(row, col)) for col in basis]
    basis = list(basis)
    for j in range(1, len(basis)):
        a, b = values[0], values[j]
        g, s, t = extended_gcd(a, b)
        if not g:
            continue
        first, other = basis[0], basis[j]
        basis[0] = tuple(s * x + t * y for x, y in zip(first, other))
        basis[j] = tuple((-b // g) * x + (a // g) * y for x, y in zip(first, other))
        values[0] = g
    multiplier = modulus // gcd(values[0], modulus)
    basis[0] = tuple(multiplier * c for c in basis[0])
    return multiplier, basis


def identity(size):
    return [tuple(int(i == j) for i in range(size)) for j in range(size)]


def gaps_form(form):
    return tuple(-6 * sum(form.row[j:]) for j in range(5)) + (
        sum((i + 1) * value for i, value in enumerate(form.row)),)


if __name__ == '__main__':
    congruences = [witness(odd)[0] for odd in range(2)]
    gap_index = 1
    base_period = 6
    gap_basis = identity(5)
    for stage, pair in enumerate(zip(*congruences), start=1):
        for numerator, divisor, remainder in pair:
            row = gaps_form(numerator)
            scale = lcm(*(c.denominator for c in (*row, numerator.constant)))
            row = tuple(int(scale * c) for c in row)
            modulus = divisor * scale
            factor, gap_basis = impose(gap_basis, row[:5], modulus)
            gap_index *= factor
            base_period = lcm(base_period, modulus // gcd(row[-1], modulus))
            print('gap increment congruence', row, 'mod', modulus,
                  'new index factor', factor)
        print('through stage', stage, 'gap index', gap_index,
              'period at fixed initial five gaps', base_period)
