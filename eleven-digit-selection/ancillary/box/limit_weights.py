"""Pull guarded checkpoint leaves back to the canonical original-input lattice."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product
import json
from math import gcd, lcm
from pathlib import Path

import coverage
from sections import prefixes, q3
import structural
import verify as route


COHORT_ZERO_SHIFT = (-5, 0, 1, 2, -1)
LANES = {'S': 'S', '1': 'M1', '2': 'M2', 'L': 'L'}
RESIDUAL_SIGNED = (0, -21, -35, 0, -12, Q(191, 3), Q(13, 3))


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def coefficients(form):
    return tuple(Q(value, route.SCALE) for value in form.v)


def source_groups():
    prefix_rows = {label: (q2, clock, kappa, len(switch))
                   for label, _, switch, _, _, q2, clock, kappa in prefixes()}
    families = defaultdict(list)
    for label, guard, state, word, drift, delay in structural.sections():
        q2, clock, kappa, switch_length = prefix_rows[label[:3]]
        r2 = label[3]
        assert kappa % 2 == 1
        # The first coordinate is the exact elapsed clock through the checkpoint.
        assert q2[0] == route.box.rows_at_switch(*label[:3])[0] - switch_length
        assert state[0] == q3(q2, clock, r2)[0] - delay
        index, rows, source = coverage.intrinsic(state, guard, label)
        word = tuple(route.paired_order(phase) for phase in word)
        linear = tuple(row.linear for row in source.pivots)
        key = index, rows, word, tuple(6 * value for value in drift), linear
        families[key].append((source, kappa))
    return families


def count_original_classes(families):
    # The columns of the inverse prefix map in (p,r,N,n,m) order.
    p = (-1, 0, 1, -1, 0)
    r = (0, 1, 0, 0, 1)
    late = (20, 4, 2, -2, 4)
    initial = (5, 4, 3, 2, 1)
    prepare = (10, 2, 0, -2, 2)
    # The three transverse differences (q1+q3, q3+q4, q2-q5)
    # kill the p and r directions.
    # Together with q3 and q2 they are a unimodular change of coordinates.
    transverse = tuple((col[0] + col[2], col[2] + col[3], col[1] - col[4])
                       for col in (late, initial, prepare))
    assert (p[2], p[1], r[2], r[1]) == (1, 0, 0, 1)
    assert all(col[0] + col[2] == col[2] + col[3] == col[1] - col[4] == 0
               for col in (p, r))
    prefix_index = abs(route.det(*transverse))
    # The original r2 splits the prefix lattice by p-2r+B mod 36.
    # The added residue rho splits it by r mod 6, independently.
    label_index = prefix_index * 36
    primary_index = label_index * 6
    labels = Counter()
    for sources in families.values():
        for source, _ in sources:
            cohort, z, r1, low = source.label
            cohorts = (0, 1) if cohort == 1 else (2,)
            labels.update(product(cohorts, (z,), (r1,), (low, low + 36)))
    expected = set(product(range(3), range(1, 3), range(-4, 18), range(601, 672, 2)))
    assert set(labels) == expected and all(count == 1 for count in labels.values())
    assert len(labels) == label_index == 4752 and primary_index == 28512
    return {'prefix_index': prefix_index, 'label_index': label_index,
            'primary_index': primary_index}


def intersect_congruence(point, columns, linear, offset, modulus):
    """Intersect an integer coset with linear·point + offset = 0 mod modulus."""
    generators = list(columns)
    values = [dot(linear, column) for column in columns]
    for j in (1, 2):
        a, b = values[0], values[j]
        if b == 0:
            continue
        g, u, v = route.extended_gcd(a, b)
        first, other = generators[0], generators[j]
        generators[0] = tuple(u*x + v*y for x, y in zip(first, other))
        generators[j] = tuple((-b//g)*x + (a//g)*y for x, y in zip(first, other))
        values[0] = g
    rhs = -(dot(linear, point) + offset)
    divisor = gcd(values[0], modulus)
    if rhs % divisor:
        return None
    index = modulus // divisor
    if index == 1:
        return point, columns
    k = (rhs // divisor) * pow(values[0] // divisor, -1, index) % index
    first = generators[0]
    point = tuple(x + k*y for x, y in zip(point, first))
    generators[0] = tuple(index*x for x in first)
    return point, tuple(generators)


def leaf_congruences(leaf):
    """The three Cramer numerators must vanish modulo the lattice determinant."""
    determinant = abs(route.det(*leaf.columns))
    standard_basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for i in range(3):
        linear = tuple(route.det(*(leaf.columns[:i] + (e,) + leaf.columns[i+1:]))
                       for e in standard_basis)
        yield linear, -dot(linear, leaf.point), determinant


def pullback(congruences, pivots, r2, kappa):
    """Preimage in (B,p,r) restricted to the canonical parent lattice."""
    point = (5, (r2-kappa)//2 - 5, 0)
    columns = ((6, -6, 0), (0, 18, 0), (0, 2, 1))
    first, second = pivots
    for (b, a, d), offset, modulus in congruences:
        composed = tuple(a*x + d*y for x, y in zip(first, second))
        u, v, w, constant = composed
        row = (b+u, v, w, offset+constant)
        denominator = lcm(*(q.denominator for q in row))
        u, v, w, constant = (int(q * denominator) for q in row)
        result = intersect_congruence(point, columns, (u, v, w),
                                      constant, modulus * denominator)
        if result is None:
            return None
        point, columns = result
    g = gcd(*(column[0] for column in columns))
    determinant = abs(route.det(*columns))
    assert g > 0 and determinant % (18*g) == 0
    return point[0] % g, g, determinant


def check_cohort_zero_shift():
    a, b, c, d, e = COHORT_ZERO_SHIFT
    chart_shift = (a, b, c, 2*(b-e), -c-d)
    shifted = route.orbit.prefix_coordinates(Q(0), *map(Q, chart_shift), 1, 0)
    p, r = route.orbit.chart_constants(1, 0)
    assert a == -5 and shifted[:2] == (p, r)


def profile_and_error(end, source, guard, terminal):
    base, a, d, constant = coefficients(end)
    leading = base*route.BASE + a*guard.leading[0] + d*guard.leading[1]
    profile = tuple(coefficients(leading)[:3])
    a_bias, d_bias = (row.at((0, *source.source_intercept)) for row in source.pivots)
    intercept = constant + a*a_bias + d*d_bias
    errors = {source.label[0]: (intercept + terminal[0], intercept + terminal[1])}
    if source.label[0] == 1:
        # The shift fixes leading P and R. Its -5 change in a cancels
        # the five extra updates of the aligned cohort-0 trajectory.
        errors[0] = errors[1]
    return profile, errors


def add_range(ranges, key, lower, upper):
    old_lower, old_upper = ranges.get(key, (lower, upper))
    ranges[key] = min(old_lower, lower), max(old_upper, upper)


def collect_leaves(families):
    """Collect integral preimages, their base projections and clock bounds."""
    family_records = []
    graph_ranges = {}
    terminal_ranges = {}
    cutoff = max_depth = empty = 0
    for key, sources in families.items():
        index, entry, word, drift, _ = key
        guard = coverage.CheckpointGuards(key, [source for source, _ in sources])
        start = coverage.physical_lattice(index, entry)
        clock = -(entry[0] + 5*route.BASE)/6
        leaf_records = []
        forms = [tuple(coefficients(row) for row in source.pivots) for source, _ in sources]
        for end, jumps, tail, leaf in route.walk_macro(entry, clock, word, drift, start, guard, ()):
            terminal_word = jumps[-1][0]
            terminal = (-len(terminal_word)*5042, 1 + max(t[0] for t in tail))
            add_range(terminal_ranges, 'physical_leaves', *terminal)
            max_depth = max(max_depth, len(jumps))
            congruences = tuple(leaf_congruences(leaf))
            intersections = []
            profiles = set()
            for source_index, ((source, kappa), pivots) in enumerate(zip(sources, forms)):
                r2 = source.label[3]
                original = pullback(congruences, pivots, r2, kappa)
                if original is None:
                    empty += 1
                    continue
                profile, errors = profile_and_error(end, source, guard, terminal)
                profiles.add(profile)
                for cohort, bounds in errors.items():
                    add_range(graph_ranges, cohort, *bounds)
                intersections.append((source_index, *original))
            assert len(profiles) <= 1
            leaf_records.append((next(iter(profiles), None), intersections))
        cutoff = max(cutoff, guard.minimum_base)
        family_records.append(([source.label for source, _ in sources], leaf_records))
    return family_records, cutoff, max_depth, empty, terminal_ranges, graph_ranges


def aggregate(family_records, primary_index):
    profiles = {profile for _, leaves in family_records for profile, _ in leaves if profile}
    ordered = sorted(profiles, key=lambda profile: (
        route.LANES.index(route.profile(route.Form.create(profile))[0]), profile[0]
    ))
    # Keep the established order of the two Long graphs (decreasing intercept).
    ordered[-2:] = reversed(ordered[-2:])
    assert len(ordered) == 7
    modulus = lcm(12, *(g for _, leaves in family_records for _, intersections in leaves
                        for _, _, g, _ in intersections))
    # Each source at a fixed base comprises 12 primary classes (2 r2 x 6 rho).
    # A viable leaf contributes primary_index * g / (132 * determinant).
    counts = defaultdict(Counter)
    for labels, leaves in family_records:
        source_totals = [Counter() for _ in labels]
        for profile, intersections in leaves:
            for source_index, base_residue, step, determinant in intersections:
                copies = 2 if labels[source_index][0] == 1 else 1
                units = Q(primary_index * step, 132 * determinant)
                for base in range(base_residue, modulus, step):
                    assert base % 6 == 5
                    source_totals[source_index][base] += units
                    counts[base][profile] += copies*units
        for total in source_totals:
            assert total == {base: 12 for base in range(5, modulus, 6)}, total
    arms = {arm: tuple(counts[arm][profile] for profile in ordered) for arm in (5, 11)}
    assert len(counts) == modulus//6
    for base, row in counts.items():
        assert tuple(row[profile] for profile in ordered) == arms[base % 12], base
        assert sum(row.values()) == primary_index
    signed = tuple(left - right for left, right in zip(arms[5], arms[11]))
    assert signed == RESIDUAL_SIGNED
    pos = sum(weight for weight in signed if weight > 0)
    neg = sum(weight for weight in signed if weight < 0)
    assert (pos, neg) == (68, -68)
    table = []
    for profile, difference in zip(ordered, signed):
        lane, _, slack, _ = route.profile(route.Form.create(profile))
        assert slack > 0
        table.append({'lane': LANES[lane], 's': list(map(str, profile)),
                      'arm_5_units': str(counts[5][profile]), 'arm_11_units': str(counts[11][profile]),
                      'signed_units': str(difference), 'frozen_test': int(lane == 'L')})
    assert sum(Q(row['signed_units']) * row['frozen_test'] for row in table) == pos
    for arm, expected in ((5, (15042, 1328, 2688, 9454)),
                          (11, (15098, 1328, 2700, 9386))):
        lane_totals = tuple(sum(Q(row[f'arm_{arm}_units']) for row in table
                                if row['lane'] == lane) for lane in LANES.values())
        assert lane_totals == expected
    return table, modulus, {'positive_units': int(pos), 'negative_units': int(neg),
                            'frozen_test_contrast': str(Q(pos, primary_index))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage_record', type=Path)
    parser.add_argument('--leaves', action='store_true', help='include the per-source leaf ledger')
    args = parser.parse_args()
    archived = json.loads(args.coverage_record.read_text())
    check_cohort_zero_shift()
    families = source_groups()
    census = count_original_classes(families)
    records, cutoff, depth, empty, terminal, graph = collect_leaves(families)
    table, modulus, contrast = aggregate(records, census['primary_index'])
    assert len(families) == archived['closed_families'] == archived['families']
    assert sum(len(leaves) for _, leaves in records) == archived['leaves']
    assert (cutoff, depth) == (archived['cutoff'], archived['max_depth'])
    assert sorted(row['s'] for row in table) == archived['leading_clock_forms']
    pullback_types = Counter((step, determinant // step)
                             for _, leaves in records for _, intersections in leaves
                             for _, _, step, determinant in intersections)
    result = {
        'census': census,
        'intrinsic_families': len(families),
        'intrinsic_leaves': archived['leaves'],
        'nonempty_pullbacks': sum(len(xs) for _, leaves in records for _, xs in leaves),
        'empty_pullbacks': empty,
        'nonempty_pullback_types': [
            {'base_step': step, 'fixed_base_covolume': covolume, 'count': count}
            for (step, covolume), count in sorted(pullback_types.items())
        ],
        'cutoff_cohorts_1_and_2': cutoff,
        'max_depth': depth,
        'base_modulus': modulus,
        'tested_base_residues': list(range(5, modulus, 6)),
        'terminal_clock_error_bounds': list(terminal['physical_leaves']),
        'graph_error_bounds_by_cohort': {cohort: list(map(str, graph[cohort]))
                                         for cohort in sorted(graph)},
        'graphs': table,
        'difference': contrast,
    }
    if args.leaves:
        result['leaf_ledger_columns'] = ['source_index', 'base_residue', 'base_step',
                                         'pullback_determinant']
        result['leaf_ledger'] = [
            {'source_labels': labels,
             'leaves': [{'s': list(map(str, profile)) if profile else None, 'pullbacks': rows}
                        for profile, rows in leaves]}
            for labels, leaves in records
        ]
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
