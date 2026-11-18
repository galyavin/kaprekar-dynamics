"""Exact first-return model for doubling a four-digit Kaprekar base.

The general proof is in doubling_theorem.md. Run with Python 3's standard
library. Direct verification uses sorted output digits from investigate.py;
the compressed model uses separately derived formulas.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json

from investigate import states, transition, published_maximum_distance

ROOT = Path(__file__).resolve().parent


def ordered(a, b):
    return (max(a, b), min(a, b))


def boundary(base, pair):
    u, v = pair
    return v == 0 or u == v or u + v == base


def boundary_states(base):
    return sorted({(u, 0) for u in range(1, base)} |
                  {(u, u) for u in range(1, base)} |
                  {ordered(u, base-u) for u in range(1, base)})


def folded_step(base, pair):
    u, v = pair
    if v == 0:
        return ordered(u-1, base-u)
    if u == v or u+v == base:
        return ordered(abs(2*u-base+1), abs(2*u-base-1))
    return ordered(abs(2*u-base), abs(2*v-base))


def return_rule(base, pair):
    """First return of 2*pair in base 2*base: (smaller pair, elapsed steps)."""
    u, v = pair
    if not boundary(base, pair):
        return ordered(abs(2*u-base), abs(2*v-base)), 1
    if v == 0:
        return ordered(abs(2*u-base-1), abs(2*u-base)), 2
    h = 2*u-base
    if h == 0 or 2*abs(h) == base:
        return (base-1, base-3), 3
    return ordered(abs(abs(2*h+1)-base), abs(abs(2*h-1)-base)), 2


def enter_even(base, pair, direct=False):
    """First entry in the even-coordinate set of base 2*base; max 3 steps."""
    step = transition if direct else folded_step
    q, elapsed = pair, 0
    while q[0] % 2 or q[1] % 2:
        q = step(2*base, q)
        elapsed += 1
        assert elapsed <= 3, (base, pair, q)
    return (q[0]//2, q[1]//2), elapsed


def multiplicity(base, pair):
    """Number of nonconstant four-slot strings with this difference pair."""
    u, v = pair
    if v == 0:
        return (base-u)*(12*u-4)
    if u == v:
        return 6*(base-u)
    return 24*(base-u)*(u-v)


def analyse(nodes, edge, target):
    """Weighted functional graph: finite target distance or None, plus cycles."""
    distance = {target: 0}
    cycle_records = []
    for start in nodes:
        if start in distance:
            continue
        path, index, q = [], {}, start
        while q not in distance and q not in index:
            index[q] = len(path)
            path.append(q)
            q = edge(q)[0]
        if q in index:
            loop = path[index[q]:]
            cycle_records.append({'states': loop,
                                  'period': sum(edge(r)[1] for r in loop)})
            for r in loop:
                distance[r] = None
            path = path[:index[q]]
        for r in reversed(path):
            nxt, cost = edge(r)
            distance[r] = None if distance[nxt] is None else cost+distance[nxt]
    return distance, cycle_records


def entry_polynomials(base):
    """Return sparse pair-count and input-count injection polynomials.

    The four generic preimages are explicit. Only O(base) boundary states
    require separate handling. No full graph in base 2*base is constructed.
    """
    counts, weights = {}, {}
    for r in states(base):
        u, v = r
        counts[r] = Counter({0: 1})
        weights[r] = Counter({0: multiplicity(2*base, (2*u, 2*v))})
        if u == v:
            continue
        preimages = {ordered(base+sgn*u, base+sgn2*v)
                     for sgn, sgn2 in product((-1, 1), repeat=2)}
        expected = 0 if (u % 2 == v % 2 == base % 2) else (4 if v else 2)
        off_even = [q for q in preimages if q[0] % 2 or q[1] % 2]
        assert len(off_even) == expected
        for q in off_even:
            assert not boundary(2*base, q)
            counts[r][1] += 1
            weights[r][1] += multiplicity(2*base, q)
    for q in boundary_states(2*base):
        if not (q[0] % 2 or q[1] % 2):
            continue
        r, elapsed = enter_even(base, q)
        assert elapsed in (2, 3)
        counts[r][elapsed] += 1
        weights[r][elapsed] += multiplicity(2*base, q)
    assert sum(sum(a.values()) for a in counts.values()) == (2*base)*(2*base+1)//2-1
    assert sum(sum(a.values()) for a in weights.values()) == (2*base)**4-2*base
    return counts, weights


def shifted_distribution(distances, injections):
    total = Counter()
    never = 0
    for r, polynomial in injections.items():
        if distances[r] is None:
            never += sum(polynomial.values())
        else:
            for elapsed, count in polynomial.items():
                total[elapsed+distances[r]] += count
    return total, never


def compressed_distributions(base):
    target = (3*base//5, base//5)
    distances, cycles = analyse(states(base), lambda q: return_rule(base, q), target)
    entries, weights = entry_polynomials(base)
    pairs, nonconverging_pairs = shifted_distribution(distances, entries)
    weighted_pairs, nonconverging_inputs = shifted_distribution(distances, weights)
    inputs = Counter({t+1: count for t, count in weighted_pairs.items()})
    inputs[1] -= 1  # the fixed number itself is already there at time zero
    inputs[0] = 1
    return pairs, inputs, nonconverging_pairs, nonconverging_inputs, cycles


def main():
    # Formula-level verification across both parities and all boundary cases.
    for base in range(5, 81):
        for r in states(base):
            q = (2*r[0], 2*r[1])
            q = transition(2*base, q)
            elapsed = 1
            while q[0] % 2 or q[1] % 2:
                q = transition(2*base, q)
                elapsed += 1
                assert elapsed <= 3
            assert return_rule(base, r) == ((q[0]//2, q[1]//2), elapsed)
        for q in states(2*base):
            assert enter_even(base, q) == enter_even(base, q, direct=True)
        entries, weights = entry_polynomials(base)
        direct_counts = {q: Counter() for q in states(base)}
        direct_weights = {q: Counter() for q in states(base)}
        for q in states(2*base):
            r, elapsed = enter_even(base, q, direct=True)
            direct_counts[r][elapsed] += 1
            direct_weights[r][elapsed] += multiplicity(2*base, q)
        assert entries == direct_counts
        assert weights == direct_weights
    print('Return map, entry times, and injection polynomials verified for every base B=5..80.')

    rows = []
    for base in [5, 10, 20, 40, 80, 160, 320, 640]:
        pairs, inputs, other_pairs, other_inputs, loops = compressed_distributions(base)
        actual = 2*base
        target = (3*actual//5, actual//5)
        distances, direct_loops = analyse(states(actual), lambda q: (transition(actual, q), 1), target)
        observed = Counter(d for d in distances.values() if d is not None)
        assert pairs == observed
        assert other_pairs == sum(d is None for d in distances.values())
        observed_inputs = Counter()
        for q, d in distances.items():
            if d is not None:
                observed_inputs[d+1] += multiplicity(actual, q)
        observed_inputs[1] -= 1
        observed_inputs[0] = 1
        assert inputs == observed_inputs
        assert other_inputs == sum(multiplicity(actual, q) for q, d in distances.items() if d is None)
        assert sorted(c['period'] for c in loops) == sorted(c['period'] for c in direct_loops)
        assert max(inputs) == published_maximum_distance(actual)
        row = {'base': actual,
               'pair_hitting_counts': dict(sorted(pairs.items())),
               'input_hitting_counts': dict(sorted(inputs.items())),
               'nonconverging_pairs': other_pairs,
               'nonconverging_inputs': other_inputs,
               'other_cycle_periods': sorted(c['period'] for c in loops),
               'maximum_finite_input_time': max(inputs),
               'converging_input_fraction': [actual**4-actual-other_inputs, actual**4]}
        rows.append(row)
        print(actual, 'max time', max(inputs), 'other pairs', other_pairs,
              'other cycles', row['other_cycle_periods'])

    # Independent weighting check from all ordered strings in small bases.
    for base in [5, 10, 20]:
        observed = Counter()
        for ds in product(range(base), repeat=4):
            a, b, c, d = sorted(ds, reverse=True)
            if a != d:
                observed[a-d, b-c] += 1
        assert all(observed[q] == multiplicity(base, q) for q in states(base))
    print('Input multiplicities verified by direct four-slot enumeration at bases 5, 10, 20.')
    (ROOT/'doubling_results.json').write_text(json.dumps(rows, indent=2)+'\n')


if __name__ == '__main__':
    main()
