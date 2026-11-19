"""Independent exact rational checks for the small gap rectangle."""

from fractions import Fraction as Q
from itertools import product
from math import factorial, prod
import subprocess
import sys

center = tuple(Q(i, 180000) for i in (26293, 27557, 28481, 40447, 29865, 27357))
radius = Q(1, 10**6)


def centered(y):
    return tuple(i - 6 * sum(y[:i]) for i in range(1, 6))


def mass_coefficients():
    # Expand the product of the six gaps in the five offsets. Integrating
    # kills every odd monomial; the product has only degree 0 or 2 survivors.
    coefficients = {}
    for select in product((0, 1), repeat=5):
        degree = sum(select)
        prefactor = prod(center[i] for i in range(5) if not select[i])
        if degree == 0:
            coefficients[5] = 2**5 * factorial(11) * prefactor * center[5]
        if degree == 1:
            coefficients[7] = coefficients.get(7, Q(0)) - 2**5 * factorial(11) * prefactor / 3
    return coefficients


def main():
    mass = mass_coefficients()
    index = prod((5, 6, 12, 12, 768))
    omega1 = sum(c * (radius / 2) ** k for k, c in mass.items()) / index
    phi_mass = sum(
        c * (radius / 4) ** k * (3 ** (k + 1) - 2 ** (k + 1)) / (k + 1)
        for k, c in mass.items()
    )
    odd_mass = Q(1571, 5184) + Q(1219, 5184)
    maximum = max(
        sum(a * b for a, b in zip(
            (Q(205, 1152), Q(715, 1152), -Q(65, 192), -Q(715, 1152), Q(1235, 1152)),
            centered(tuple(a + radius * b for a, b in zip(center, signs))),
        ))
        for signs in product((-1, 1), repeat=5)
    )
    print(f"T1 index {index}; lower bound {float(omega1):.14g}; odd leading maximum {maximum}")
    print(f"T2 odd mass {odd_mass}; phi integral {float(phi_mass):.14g}; contrast {float(odd_mass * phi_mass):.14g}")


def parity_witness(binary):
    rows = (
        (5, 5, 5, 11, 473),
        (0, 6, 6, 0, 642),
        (0, 0, 12, 0, 120),
        (0, 0, 0, 12, 480),
        (0, 0, 0, 0, 768),
    )

    def representatives(base, point, i=0):
        if i == 5:
            yield point
            return
        aim = base * center[i]
        low = (aim - radius * base / 2 - point[i]) / rows[i][i]
        high = (aim + radius * base / 2 - point[i]) / rows[i][i]
        for multiplier in range(-(-low.numerator // low.denominator), high.numerator // high.denominator + 1):
            yield from representatives(base, tuple(x + multiplier * v for x, v in zip(point, rows[i])), i + 1)

    chosen = 12_000_000
    first_base = chosen + (3360 - chosen) % 3840
    origin = tuple(int(a * 180000) for a in center[:5])
    for base in range(first_base, first_base + 3840 * 100, 3840):
        g = next(representatives(base, origin), None)
        if g is not None:
            break
    else:
        raise ValueError("no interior lattice point in the selected bases")
    x = tuple(i * base - 6 * sum(g[:i]) for i in range(1, 6))
    even = sum(a * b for a, b in zip((Q(3, 16), Q(55, 48), -Q(5, 8), -Q(55, 48), Q(95, 48)), x)) + 40
    odd = sum(a * b for a, b in zip(
        (Q(205, 1152), Q(715, 1152), -Q(65, 192), -Q(715, 1152), Q(1235, 1152)), x
    )) + Q(233, 24)
    for b, gaps, expected in ((base, g, even), (base + 1, g[:4] + (g[4] + 1,), odd)):
        q = tuple(b - sum(gaps[:i]) for i in range(1, 6))
        clock = int(subprocess.check_output([binary, "clock", str(b), *map(str, q)]))
        assert clock == expected, (b, q, clock, expected)
        print(f"T1 base {b}; q {q}; direct word clock {clock}; formula {expected}")


if __name__ == "__main__":
    main()
    if len(sys.argv) == 2:
        parity_witness(sys.argv[1])
