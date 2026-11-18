"""Five-slot Kaprekar dynamics, exact input weights, and finite reference data.

The number map is defined by actual sorted subtraction.  Pair formulas follow
Prichett (1978); finite graph enumeration is a reference, not a limit theorem.
"""
from collections import Counter


def valid(base, pair):
    u, v = pair
    return isinstance(base, int) and base >= 2 and 1 <= u < base and 0 <= v <= u


def states(base):
    if not isinstance(base, int) or base < 2:
        raise ValueError('Base must be an integer at least two.')
    return ((u, v) for u in range(1, base) for v in range(u + 1))


def digits(base, number):
    if not 0 <= number < base**5:
        raise ValueError('Number must fit in five slots.')
    out = []
    for _ in range(5):
        number, digit = divmod(number, base)
        out.append(digit)
    return tuple(reversed(out))


def difference_pair(ds):
    a, b, _, d, e = sorted(ds, reverse=True)
    return a - e, b - d


def subtract(base, number):
    ds = sorted(digits(base, number))
    up = down = 0
    for lo, hi in zip(ds, reversed(ds)):
        up = base * up + lo
        down = base * down + hi
    return down - up


def number_from_pair(base, pair):
    u, v = pair
    return (base**4 - 1) * u + (base**3 - base) * v


def transition(base, pair):
    """Exact borrowed digits, retaining all zero and equality cases."""
    u, v = pair
    if not valid(base, pair):
        raise ValueError('Invalid base or pair.')
    if v:
        ds = (u, v - 1, base - 1, base - v - 1, base - u)
    else:
        ds = (u - 1, base - 1, base - 1, base - 1, base - u)
    return difference_pair(ds)


def multiplicity(base, pair):
    u, v = pair
    if not valid(base, pair):
        raise ValueError('Invalid base or pair.')
    if v == 0:
        factor = 20 * u - 10
    elif u == v:
        factor = 30 * u - 10
    else:
        factor = 120 * v * (u - v) - 20
    return (base - u) * factor


def canonical_cycle(cycle):
    cycle = tuple(cycle)
    i = cycle.index(min(cycle))
    return cycle[i:] + cycle[:i]


def classical_cycles(base):
    """Prichett Tables 4.1--4.3; T1--T5 labels, small irregular exceptions.

    The regular non-T2 cycles shift by (2,1) when the base increases by 3.
    Their labels are taken from the original table, not inferred from period.
    """
    seeds = {
        2: {'I1': [(1, 0)], 'I2': [(1, 1)]},
        3: {'T1': [(2, 1)]},
        4: {'T2': [(3, 0), (2, 1)]},
        5: {'T5': [(4, 1), (4, 3), (3, 3), (3, 1)]},
        6: {'T1': [(4, 2)], 'T2': [(3, 2), (4, 0)]},
        7: {'I': [(4, 3), (4, 1), (6, 2), (5, 5), (5, 3)]},
        8: {'T5': [(6, 4), (5, 3), (5, 2), (6, 2)], 'T2': [(4, 3), (5, 0)]},
        9: {'T1': [(6, 3)], 'I': [(7, 5), (6, 4), (5, 3), (6, 1), (8, 4)]},
        10: {'T4': [(7, 5), (6, 3), (7, 2), (8, 4)],
             'I': [(6, 4), (6, 2), (8, 3), (7, 6)], 'T2': [(5, 4), (6, 0)]},
        11: {'T5': [(7, 4), (7, 3), (8, 3), (8, 5)]},
        12: {'T1': [(8, 4)], 'T2': [(6, 5), (7, 0)]},
        13: {'T4': [(9, 6), (8, 4), (9, 3), (10, 5)]},
    }
    if base < 2:
        raise ValueError('Base must be at least two.')
    if base < 14:
        raw = seeds[base]
    else:
        raw = {}
        if base % 3 == 0:
            raw['T1'] = [(2 * base // 3, base // 3)]
        elif base % 3 == 1:
            shift = (base - 16) // 3
            raw['T3'] = [(u + 2 * shift, v + shift)
                         for u, v in [(10, 4), (12, 5), (11, 8)]]
            raw['T4'] = [(u + 2 * shift, v + shift)
                         for u, v in [(10, 5), (11, 4), (12, 6), (11, 7)]]
        else:
            shift = (base - 14) // 3
            raw['T5'] = [(u + 2 * shift, v + shift)
                         for u, v in [(10, 4), (10, 6), (9, 5), (9, 4)]]
        if base % 2 == 0:
            raw['T2'] = [(base // 2, base // 2 - 1), (base // 2 + 1, 0)]
    return {label: canonical_cycle(cycle) for label, cycle in raw.items()}


def finite_graph(base, step=transition):
    """Memoized complete functional graph: canonical cycles and pair distances."""
    fate, found = {}, set()
    for start in states(base):
        if start in fate:
            continue
        path, index, q = [], {}, start
        while q not in fate and q not in index:
            index[q] = len(path)
            path.append(q)
            q = step(base, q)
        if q in index:
            cut = index[q]
            cycle = canonical_cycle(path[cut:])
            found.add(cycle)
            for p in path[cut:]:
                fate[p] = cycle, 0
            path = path[:cut]
        for p in reversed(path):
            cycle, distance = fate[step(base, p)]
            fate[p] = cycle, distance + 1
    return tuple(sorted(found)), fate


def distribution(base):
    cycles, fate = finite_graph(base)
    labels = {c: label for label, c in classical_cycles(base).items()}
    if set(cycles) != set(labels):
        raise AssertionError('Direct graph disagrees with classical cycle table.')
    hist = {label: Counter() for label in labels.values()}
    for pair, (cycle, distance) in fate.items():
        hist[labels[cycle]][distance + 1] += multiplicity(base, pair)
    for cycle, label in labels.items():
        hist[label][1] -= len(cycle)
        hist[label][0] += len(cycle)
    counts = {label: dict(sorted((t, w) for t, w in h.items() if w))
              for label, h in sorted(hist.items())}
    total = sum(sum(h.values()) for h in counts.values())
    assert total == base**5 - base
    return {'base': base, 'residue_mod_6': base % 6, 'total_input_weight': total,
            'cycles': {k: list(v) for k, v in sorted(classical_cycles(base).items())},
            'input_arrival_counts': counts,
            'first_moment_sum': sum(t * w for h in counts.values() for t, w in h.items()),
            'second_moment_sum': sum(t*t * w for h in counts.values() for t, w in h.items()),
            'maximum_input_time': max(t for h in counts.values() for t in h)}
