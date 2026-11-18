from collections import defaultdict
from fractions import Fraction
from symbolic_arrivals import residue_data


def order_two(core):
    value, order = 2 % core, 1
    while value != 1:
        value = 2*value % core
        order += 1
    return order


def phase_coefficients(core, phase, data):
    s, half = pow(2, phase, core), (core+1)//2
    coefficients, edges = {}, {}
    for name, residue, kind in (('S', s, 'Z'), ('E', 2*s % core, 'Z'),
                                ('A', half*s % core, 'A')):
        seeds = data['zero'] if kind == 'Z' else data['diagonal']
        if residue in seeds:
            coefficients[name] = 1
        elif kind == 'Z' and residue == half:
            edges[name] = 'A'
        elif kind == 'Z' and residue == 1:
            edges[name] = 'S'
        elif kind == 'A' and residue in ((core-1)//2, half):
            edges[name] = 'E'
        else:
            coefficients[name] = None
    for start in ('S', 'E', 'A'):
        path, seen, q = [], set(), start
        while q not in coefficients and q not in seen:
            path.append(q)
            seen.add(q)
            q = edges[q]
        value = coefficients.get(q)
        for q in reversed(path):
            if value is not None:
                value += 1
            coefficients[q] = value
    assert all(k is None or 1 <= k <= 3 for k in coefficients.values())
    return coefficients


def power_coefficient(core, power, kind, data, coefficients):
    residue = pow(2, power, core)
    seeds = data['zero'] if kind == 'Z' else data['diagonal']
    if residue in seeds:
        return 1
    parent = None
    if kind == 'Z':
        if residue == (core+1)//2:
            parent = 'A'
        elif residue == 1:
            parent = 'S'
    elif residue in ((core-1)//2, (core+1)//2):
        parent = 'E'
    if parent is None or coefficients[parent] is None:
        return None
    return coefficients[parent]+1


def limiting_law(core, phase, target, data=None):
    if data is None:
        data = residue_data(core, target)
    period = order_two(core)
    phase %= period
    k = phase_coefficients(core, phase, data)
    masses = defaultdict(Fraction)
    masses[1] += Fraction(len(data['residues']), core*core)
    masses[2] += Fraction(2*(len(data['zero'])+len(data['diagonal'])), core*core)
    for name, numerator in (('A', 2), ('E', 4), ('S', 2)):
        if k[name] is not None:
            masses[k[name]+2] += Fraction(numerator, core*core)
    denominator = 2**period-1
    for power in range(1, period+1):
        coefficient = power_coefficient(core, power, 'Z', data, k)
        if coefficient is not None:
            masses[coefficient+1] += Fraction(2*2**(period-power),
                                             3*core*core*denominator)
    for power in range(period):
        coefficient = power_coefficient(core, power, 'A', data, k)
        if coefficient is not None:
            masses[coefficient+1] += Fraction(2**(period-power),
                                             6*core*core*denominator)
    masses = dict(sorted((t, mass) for t, mass in masses.items() if mass))
    total = sum(masses.values(), Fraction())
    assert 0 <= total <= 1
    assert all(1 <= t <= 5 for t in masses)
    first = sum(t*mass for t, mass in masses.items())
    second = sum(t*t*mass for t, mass in masses.items())
    mean, second_moment = (first/total, second/total) if total else (None, None)
    return {'core': core, 'phase': phase, 'period': period,
            'critical_coefficients': k, 'atoms': masses,
            'basin_limit': total, 'complement_limit': 1-total,
            'conditional_mean_coefficient': mean,
            'conditional_second_moment_coefficient': second_moment,
            'conditional_variance_coefficient':
                second_moment-mean*mean if total else None}
