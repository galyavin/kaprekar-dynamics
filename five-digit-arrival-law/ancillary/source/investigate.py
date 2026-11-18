"""Reproduce the finite computations in report.md; Python 3, no dependencies."""
from collections import Counter
import argparse


def output_digits(base, pair):
    u, v = pair
    if v:
        return (u, v - 1, base - 1 - v, base - u)
    return (u - 1, base - 1, base - 1, base - u)


def transition(base, pair):
    a, b, c, d = sorted(output_digits(base, pair), reverse=True)
    return a - d, b - c


def states(base):
    return [(u, v) for u in range(1, base) for v in range(u + 1)]


def fixed_pairs(base):
    return [p for p in states(base) if transition(base, p) == p]


def integer_cube_root(n):
    # Integer Newton iteration from an upper bound; no floating point.
    x = 1 << ((n.bit_length() + 2) // 3)
    while True:
        y = (2 * x + n // (x * x)) // 3
        if y >= x:
            return x
        x = y


def basin_stats(base):
    target = (3 * base // 5, base // 5)
    depths = []
    cycles = set()
    for start in states(base):
        p, path, visited = start, [], {}
        while p != target and p not in visited:
            visited[p] = len(path)
            path.append(p)
            p = transition(base, p)
        if p == target:
            depths.append(len(path))
        else:
            cycles.add(tuple(sorted(path[visited[p]:])))
    return len(depths), len(states(base)) - len(depths), max(depths) + 1, sorted(map(len, cycles))


def published_maximum_distance(base):
    """Devlin--Zeng (2021), Theorem 1, for bases divisible by five.

    https://math.colgate.edu/~integers/v97/v97.pdf
    Maximum finite hitting time; excludes trajectories entering other cycles.
    """
    if base < 5 or base % 5:
        raise ValueError('This function requires a positive base divisible by five.')
    odd_part, exponent = base // 5, 0
    while odd_part % 2 == 0:
        odd_part //= 2
        exponent += 1
    if odd_part > 1:
        return exponent + 2
    if exponent < 3:
        return (4, 7, 10)[exponent]
    if exponent % 4 == 0:
        return 4 * exponent + 6
    if exponent % 4 == 3:
        return 5 * exponent + 6
    return 3 * exponent + 5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cube-limit', type=int, default=1_000_000)
    args = parser.parse_args()

    # Independent finite check of the fixed-point classification through base 200.
    for base in range(2, 201):
        expected = ([(1, 0), (1, 1)] if base == 2 else
                    [(3, 1)] if base == 4 else
                    [(3 * base // 5, base // 5)] if base % 5 == 0 else [])
        assert fixed_pairs(base) == expected, base
    print('Fixed-point classification checked through base 200.')

    depths = Counter()
    farthest = []
    for start in states(10):
        p, steps, seen = start, 0, set()
        while p != (6, 2):
            assert p not in seen
            seen.add(p)
            p = transition(10, p)
            steps += 1
        depths[steps] += 1
        if steps == 6:
            farthest.append(start)
    print('Decimal pair distances:', dict(sorted(depths.items())))
    print('Farthest decimal pairs:', farthest)

    # Direct check on all 9,990 nonconstant decimal strings.
    def kaprekar(n):
        s = sorted(f'{n:04d}')
        return int(''.join(reversed(s))) - int(''.join(s))

    counts = Counter()
    for n in range(10_000):
        if len(set(f'{n:04d}')) == 1:
            continue
        steps = 0
        while n != 6174:
            n = kaprekar(n)
            steps += 1
            assert steps <= 7
        counts[steps] += 1
    print('Decimal number distances:', dict(sorted(counts.items())))

    print('Base, pairs reaching family fixed point, other pairs, max steps in basin, other cycle lengths:')
    for base in (5, 10, 15, 20, 40, 80, 160):
        stats = basin_stats(base)
        assert stats[2] == published_maximum_distance(base)
        print(base, stats)
    print('All displayed maximum distances agree with Devlin--Zeng Theorem 1.')

    # The proof covers all larger bases by inequalities; these are its finite cases.
    equality_bases = []
    for base in (5, 10, 20, 40, 80):
        numerator = base * (3 * base * base + 4 * base + 3)
        if numerator == 10 * published_maximum_distance(base) ** 3:
            equality_bases.append(base)
    assert equality_bases == [10]
    print('Cube/hitting-time corollary: finite exceptional family cases leave only base 10.')

    # Refute the universal cube/depth identity in globally convergent bases.
    for base, expected_quotient, expected_steps in ((5, 49, 4), (40, 19852, 21)):
        reached, other, steps, cycles = basin_stats(base)
        assert other == 0 and not cycles
        assert reached == len(states(base)) and steps == expected_steps
        numerator = base * (3 * base * base + 4 * base + 3)
        assert numerator == 10 * expected_quotient
        assert expected_quotient != steps ** 3
    print('Universal cube/depth equality refuted in globally convergent bases 5 and 40.')

    cubes = []
    for k in range(1, args.cube_limit + 1):
        numerator = k * (75 * k * k + 20 * k + 3)
        assert numerator % 2 == 0
        quotient = numerator // 2
        y = integer_cube_root(quotient)
        if y ** 3 == quotient:
            cubes.append((k, 5 * k, y))
    print(f'Cube quotients through k={args.cube_limit}: (k, base, cube root) =', cubes)


if __name__ == '__main__':
    main()
