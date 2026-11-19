"""Exhaust even routes on the closed gap rectangle."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil, gcd, lcm

from affine_routes import (INITIAL_X, VARIABLES, ZERO, compare, guards,
                           normalized_interval, repeat, translate)
from cohort import extended_gcd
from primitive import CENTER, step

RAMP_MAX = Q(131, 5000)
CYCLE_ENVELOPE = 504
PRETERMINAL_RETURNS = 30
MAX_SHORT = 1000
MAX_BLOCKS = 12

base_threshold = ceil(66 / (min(Q(g, 180000) for g in CENTER)
                           - Q(5, 10**6) - Q(1, 7)))
metrics = Counter()
bounds = {}
families = defaultdict(set)
sheet_fibers = defaultdict(set)
terminal_fibers = defaultdict(set)
routes = set()


def positive(form):
    global base_threshold
    if not any(form.row):
        return form.constant >= 0
    lo, hi = normalized_interval(form.row)
    if lo > 0:
        margin = lo
        result = True
    elif hi < 0:
        margin = -hi
        result = False
    else:
        raise ArithmeticError(f'undecided guard {form} in [{lo}, {hi}]')
    base_threshold = max(base_threshold, ceil((abs(form.constant) + 1) / margin))
    return result


def choose(state):
    return sum(positive(guard) << i for i, guard in enumerate(compare(state)))


def probe_step(src, unit):
    w, y, z, t, u = src
    pairs = ((w, -u - unit), (y, -t - unit), (z, -z - unit),
             (t, -y - unit), (u - unit, -w))
    rank = [-unit]
    code = 0
    for i, (left, right) in enumerate(pairs):
        if left >= right:
            rank.extend((left, right))
            code |= 1 << i
        else:
            rank.extend((right, left))
    return tuple(rank[i] - rank[10 - i] for i in range(5)), code


def integer_probe(state):
    """Propose long returns using the center with a very large base."""
    x = [f.evaluate(INITIAL_X) * 10**12 - f.constant * (10**12 - 1) for f in state]
    m = lcm(*(value.denominator for value in x))
    start = tuple(int(v * m) for v in x)
    src, choices = start, []
    for _ in range(65):
        src, branch = probe_step(src, 6 * m)
        choices.append(branch)
        shift = tuple(v - w for v, w in zip(src, start))
        if max(abs(v) for v in shift) >= 10**5 * m:
            continue
        target = tuple(v + 128 * d for v, d in zip(start, shift))
        for code in choices:
            target, got = probe_step(target, 6 * m)
            if code != got:
                break
        else:
            yield tuple(choices)


@lru_cache(maxsize=100000)
def macro(state):
    cur, pulled, word = state, [], []
    for suggestion in integer_probe(state):
        matches = True
        for choice in suggestion[len(word):]:
            exact = choose(cur)
            if exact != choice:
                matches = False
                break
            pulled.extend(guards(cur, choice))
            word.append(choice)
            cur = step(cur, choice)[0]
        if not matches:
            break
        drift = tuple(f - g for f, g in zip(cur, state))
        if any(any(d.row) for d in drift):
            continue
        delta = tuple(int(d.constant) for d in drift)
        if not any(delta):
            raise ArithmeticError(f'unexpected nonzero-leading cycle {state}')
        slope, returned = translate(suggestion, delta)
        if returned != delta:
            continue
        declining = [(-d, g) for g, d in zip(pulled, slope) if d < 0]
        if not declining:
            raise ArithmeticError('unexpected unbounded translation')
        div, num = min(declining, key=lambda pair: (
            pair[1].leading.evaluate(INITIAL_X) / pair[0],
            pair[1].constant / pair[0]))
        for divisor, guard in declining:
            difference = guard / divisor - num / div
            if not positive(difference):
                raise ArithmeticError(f'earliest exit changes on rectangle {state}')
        if any(num.row):
            return suggestion, delta, div, num
    return (choose(state),), (), 0, ZERO


def dot(row, values):
    return sum(a * b for a, b in zip(row, values))


def remainders(num, modulus, x, basis):
    """Enumerate exact projected cosets, including each new parametrization."""
    value = num.evaluate(x)
    assert value.denominator == 1
    residue = int(value) % modulus
    v = [dot(num.row, c) for c in basis]
    assert all(t.denominator == 1 for t in v)
    v = [int(t) % modulus for t in v]
    basis = list(basis)
    for j in range(1, len(basis)):
        if not v[j]:
            continue
        a, b = v[0], v[j]
        g, s, t = extended_gcd(a, b)
        first, other = basis[0], basis[j]
        basis[0] = tuple(s * p + t * q for p, q in zip(first, other))
        basis[j] = tuple((-b // g) * p + (a // g) * q for p, q in zip(first, other))
        v[0] = g

    g = gcd(v[0], modulus)
    period = modulus // g
    shift = basis[0]
    basis[0] = tuple(k * period for k in shift)
    basis = tuple(basis)
    inv = pow(v[0] // g, -1, period) if period > 1 else 0
    for rem in range(residue % g, modulus, g):
        n = ((rem - residue) // g * inv) % period
        point = tuple(p + n * k for p, k in zip(x, shift))
        assert num.evaluate(point) % modulus == rem
        yield rem, point, basis


def plain(state):
    assert all(not any(s.row) and s.constant.denominator == 1 for s in state)
    return tuple(int(s.constant) for s in state)


def residual(state):
    gaps = (-state[0], *(a - b for a, b in zip(state, state[1:])), state[-1])
    return -min(gaps)


def record_terminal(next_state, elapsed, word, delta):
    global base_threshold
    endpoint = plain(next_state)
    spot = tuple(x - PRETERMINAL_RETURNS * d for x, d in zip(endpoint, delta))
    assert residual(spot) > CYCLE_ENVELOPE
    cut = elapsed - PRETERMINAL_RETURNS * len(word)
    lower, upper = normalized_interval(cut.row)
    assert lower > RAMP_MAX, (lower, upper, cut, spot)
    key = tuple(cut.row)
    previous = bounds.get(key, (cut.constant, cut.constant))
    bounds[key] = min(cut.constant, previous[0]), max(cut.constant, previous[1])
    base_threshold = max(base_threshold, ceil((abs(cut.constant) + 1) / (lower - RAMP_MAX)))
    terminal_fibers[(len(word), delta)].add(endpoint)
    metrics['leaf'] += 1


def explore(state, elapsed, x, basis, depth=0, history=()):
    metrics['entry', depth] += 1
    for gap in range(MAX_SHORT):
        word, delta, modulus, num = macro(state)
        if modulus:
            break
        state, elapsed = step(state, word[0])[0], elapsed + 1
    else:
        raise ArithmeticError(f'no substantial return at depth {depth}')
    assert depth < MAX_BLOCKS, history
    speed = (len(word), delta)
    metrics['blocks', depth, speed] += 1
    families[(depth, speed)].add((num.row, modulus))
    anchor, radial = state[:2]
    if (not any(state[2].row) and anchor.row == state[3].row
            and radial.row == state[4].row):
        fiber = plain((state[2], state[3] - anchor, state[4] - radial))
        sheet_fibers[speed].add(fiber)
    for remainder, point, next_basis in remainders(num, modulus, x, basis):
        assert positive((num - remainder) / modulus)
        next_state, time = repeat(state, elapsed, (word, delta, modulus, num), remainder)
        if all(not any(s.row) for s in next_state):
            assert positive((num - remainder) / modulus + 1 - PRETERMINAL_RETURNS)
            record_terminal(next_state, time, word, delta)
            routes.add(history + (speed,))
        else:
            explore(next_state, time, point, next_basis, depth + 1, history + (speed,))


if __name__ == '__main__':
    basis = tuple(tuple(6 * int(i == j) for j in range(5)) for i in range(5))
    explore(VARIABLES, ZERO, (0,) * 5, basis)
    for metric, value in sorted(metrics.items(), key=lambda p: repr(p[0])):
        print(metric, value)
    print('uniform base threshold', base_threshold)
    for row, offsets in sorted(bounds.items()):
        lo, hi = normalized_interval(row)
        print('leading time', row, 'range', (str(lo), str(hi)), 'constants', offsets)
    for key, values in sorted(sheet_fibers.items()):
        print('sheet', key, 'fibers', sorted(values))
    for key, values in sorted(terminal_fibers.items()):
        print('terminal', key, 'endpoints', sorted(values))
    for route in sorted(routes):
        print('route', route)
    print(macro.cache_info())
