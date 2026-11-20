#!/usr/bin/env python3
"""Independent first-exit and marked-key checks with whole-word subtraction."""

import csv
import gzip
from pathlib import Path
import sys


PERMUTATION = (0, 6, 1, 5, 2, 4, 3)
KICK = (1, -2, 1, 0, 0, 0, 0)
DRIFT = (6, -1, -1, -1, -1, -1, -1)


def add(a, b, scale=1):
    return tuple(x + scale * y for x, y in zip(a, b))


def left(y):
    return tuple(y[j] + kick for j, kick in zip(PERMUTATION, KICK))


def sums(y):
    return tuple(sum(y[i:]) + sum(y[7 - i:]) for i in range(1, 7))


def differences(base, centered):
    assert all((x + base) % 7 == 0 for x in centered)
    gaps = tuple((x + base) // 7 for x in centered)
    assert sum(gaps) == base and min(gaps) >= 0 and gaps[0] > 0
    return tuple(sum(gaps[i:]) for i in range(1, 7))


def centered(base, q):
    return tuple(7 * (a - b) - base for a, b in zip((base, *q), (*q, 0)))


def advance(base, q):
    value = sum(x * (base ** (12 - i) - base ** i) for i, x in enumerate(q))
    digits = []
    for _ in range(13):
        value, digit = divmod(value, base)
        digits.append(digit)
    assert value == 0
    digits.sort()
    return tuple(digits[-i - 1] - digits[i] for i in range(6))


def check_exit():
    # Strict phase-zero centre-pair chart from the archived S1c nomination.
    direction = (-440, 183, 155, 78, -62, 29, 57)
    phase_directions = [direction]
    for _ in range(5):
        phase_directions.append(tuple(phase_directions[-1][j] for j in PERMUTATION))
    leading = tuple(sums(y) for y in phase_directions)
    assert all(value > 0 for row in leading for value in row)
    assert leading[0][2:4] == (126, 126)
    assert min(value for e, row in enumerate(leading) for i, value in enumerate(row)
               if (e, i) not in {(0, 2), (0, 3)}) == 343
    assert sums(DRIFT) == (-7,) * 6

    # In the six free coordinates, z_3 has coefficient one. This gives an
    # integral basis for ker(L_3), and translation along DRIFT subtracts 7.
    coeff = (0, 0, 1, 2, 2, 2)
    basis = [tuple((i == j) - coeff[j] * (i == 2) for i in range(6))
             for j in range(6) if j != 2]
    assert all(sum(a * b for a, b in zip(coeff, col)) == 0 for col in basis)
    assert tuple(tuple(col[j] for j in (0, 1, 3, 4, 5)) for col in basis) == tuple(
        tuple(int(i == j) for i in range(5)) for j in range(5))

    base = 8 * 14342
    beta = base % 7
    g_start = tuple(8 * x for x in direction)
    n = tuple((x + beta) // 7 for x in g_start)
    assert tuple(7 * x - beta for x in n) == g_start
    candidates = []
    y = n
    for phase in range(6):
        for pair, val in enumerate(sums(y)):
            h = val - beta + (1 if pair < 5 else -1)
            assert h >= 0
            candidates.append((6 * (1 + h // 7) + phase, pair, h % 7))
        y = left(y)
    assert y == add(n, DRIFT)
    clock, pair, eta = min(candidates)
    J, phase = divmod(clock, 6)
    assert (pair, phase) == (2, 0)
    exit_y = add(n, DRIFT, J)
    assert sum(exit_y) == beta and sums(exit_y)[pair] == beta + eta - 8
    assert add(exit_y, DRIFT, -J) == n
    q = differences(base, g_start)
    proposed = n
    for t in range(clock + 1):
        assert centered(base, q) == tuple(7 * x - beta for x in proposed)
        assert base >= 8 * -min(centered(base, q)) + 7
        if t < clock:
            q, proposed = advance(base, q), left(proposed)
    assert advance(base, q) != differences(base, tuple(7 * x - beta for x in left(proposed)))
    print(f'first-exit cell: B={base}, phase={phase}, pair={pair + 1}, J={J}, '
          f'eta={eta}, elapsed={clock}, L3={sums(exit_y)[pair]}, rank=5, '
          f'exit y={exit_y}; literal slide and first changed update agree')


def check_keys(archive):
    wanted = {(140000, 0), (140000, 81), (140006, 0)}
    with gzip.open(archive / 'encounters.tsv.gz', 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            key = int(row['base']), int(row['draw_id'])
            if key not in wanted or row['section'] != 'R2000' or row['event'] != 'entrance':
                continue
            wanted.remove(key)
            gap = tuple(map(int, row['tangent'].split(',')))
            for base in (key[0], 32775 + (key[0] - 32775) % 7):
                q = differences(base, gap)
                seen = {}
                inner = None
                peak = -min(gap)
                while q not in seen:
                    seen[q] = len(seen)
                    r = -min(centered(base, q))
                    peak = max(r, peak)
                    assert base >= 8 * r + 7 and len(seen) <= 32768 and peak <= 4096
                    if r <= 34 and inner is None:
                        inner = seen[q]
                    q = advance(base, q)
                period = len(seen) - seen[q]
                assert period == int(row['mark']) and inner is not None
                print(f'outer key {key} at B={base}: period={period}, first R<=34 offset={inner}, '
                      f'repeat={len(seen)}, peak={peak}, key={gap}')
    assert not wanted


if __name__ == '__main__':
    check_exit()
    check_keys(Path(sys.argv[1]))
