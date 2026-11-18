from collections import Counter
from dataclasses import dataclass
from basin_formula import moments, v2
from cycle_model import summarize
from doubling_checks import multiplicity
from investigate import transition


def ordered(a, b):
    return max(a, b), min(a, b)


def odd_representative(residue, core):
    residue %= core
    return residue if residue % 2 else core-residue


def residue_data(core, target):
    """Signed unordered residue representatives of the stable core target."""
    if core < 3 or core % 2 != 1:
        raise ValueError('Core must be odd and at least 3.')
    target = frozenset(target)
    if any(not (0 < v < u < core and u % 2 and v % 2) for u, v in target):
        raise ValueError('Target must lie in the odd, distinct, positive stable region.')
    if {transition(core, pair) for pair in target} != target:
        raise ValueError('Core target must be a union of cycles.')
    residues = {(a, b) for a in range(1, core) for b in range(1, core)
                if a != b and (a+b) % core
                and ordered(odd_representative(a, core),
                            odd_representative(b, core)) in target}
    assert len(residues) == 8*len(target)
    assert {((2*a) % core, (2*b) % core) for a, b in residues} == residues
    diagonal = {u for u in range(core) if ((2*u+1) % core, (2*u-1) % core) in residues}
    zero = {u for u in range(core) if ((u-1) % core, (-u) % core) in residues}
    assert not (zero & {0, 1, (core+1)//2})
    assert not (diagonal & {0, (core-1)//2, (core+1)//2})
    coefficient4 = len(residues)*core*core
    coefficient2 = 6*sum(core*abs(a-b)-(a-b)**2 for a, b in residues)-coefficient4
    target_weight = sum(multiplicity(core, pair) for pair in target)
    return {'residues': residues, 'diagonal': diagonal, 'zero': zero,
            'coefficient4': coefficient4, 'coefficient2': coefficient2,
            'target_weight': target_weight, 'target_count': len(target)}


@dataclass(frozen=True)
class Rule:
    target: tuple
    family: str
    parameter: int
    elapsed: int


def boundary_rules(core, exponent):
    if exponent < 1:
        raise ValueError('Symbolic boundary rules require n>=1.')
    scale, base = 1 << exponent, core << exponent
    rules = [Rule(('D', 1), 'midpoint', 0, 1),
             Rule(('Z', base-1), 'one', 0, 1)]
    rules.extend(Rule(('Z', 1 << t), 'zero_block', t, t+1)
                 for t in range(1, exponent+1))
    rules.append(Rule(('Z', 2*scale), 'diagonal_residues', 0, exponent+1))
    rules.extend(Rule(('A', base//2+(1 << t)), 'diagonal_valuation', t, t+1)
                 for t in range(exponent-1))
    rules.append(Rule(('A', (core+1)*scale//2), 'zero_residue', 0, exponent+1))
    return rules


def in_rule(core, exponent, point, rule):
    kind, u = point
    base = core << exponent
    if rule.family == 'midpoint':
        return kind == 'D' and u == base//2
    if rule.family == 'one':
        return kind == 'Z' and u == 1
    if rule.family == 'zero_block':
        if kind != 'Z':
            return False
        step = base >> rule.parameter
        return any(x > 0 and x % step == 0 and (x//step) % 2
                   for x in (u, u-1))
    if rule.family == 'diagonal_residues':
        return kind in ('D', 'A') and u % core in ((core-1)//2, (core+1)//2)
    if rule.family == 'diagonal_valuation':
        return (kind in ('D', 'A') and u % core == 0
                and v2(u) == exponent-rule.parameter-2)
    if rule.family == 'zero_residue':
        return kind == 'Z' and u % core == (core+1)//2
    raise AssertionError(rule)


def critical_times(core, exponent, data):
    """Solve the finite deterministic graph of critical boundary targets."""
    rules = boundary_rules(core, exponent)
    times, edges = {}, {}
    for rule in rules:
        q = rule.target
        kind, u = q
        seeds = data['zero'] if kind == 'Z' else data['diagonal']
        parents = [r for r in rules if in_rule(core, exponent, q, r)]
        assert len(parents) <= 1
        assert not (u % core in seeds and parents)
        if u % core in seeds:
            times[q] = exponent+2
        elif parents:
            edges[q] = (parents[0].target, parents[0].elapsed)
        else:
            times[q] = None
    for start in (r.target for r in rules):
        path, seen, q = [], set(), start
        while q not in times and q not in seen:
            seen.add(q)
            path.append(q)
            q = edges[q][0]
        value = times.get(q)
        for q in reversed(path):
            if value is not None:
                value += edges[q][1]
            times[q] = value
    return rules, times


def residue_stats(base, core, residue, valuation, lower=1):
    power = 1 << valuation
    odd = (pow(power, -1, core)*residue) % core
    if not odd % 2:
        odd += core
    return moments(power*odd, 2*core*power, lower, base-1)


def add_tree(powers, counts, base, kind, depth, stats, arrival):
    """Sum a family of generic predecessor trees at a fixed boundary time."""
    c, u, u2 = stats
    if not c:
        return
    if kind == 'D':
        counts[arrival] += 6*(base*c-u)
        return
    if kind == 'Z':
        counts[arrival] += -4*base*c+(12*base+4)*u-12*u2
        for t in range(1, depth+1):
            counts[arrival+t] += base*base*(powers[t]-4)*c+12*base*u-6*u2
        return
    if kind == 'A':
        counts[arrival] += -24*base*base*c+72*base*u-48*u2
        for t in range(1, depth+1):
            counts[arrival+t] += 2*base*base*(powers[t]-10)*c+48*base*u-24*u2
        return
    raise ValueError(kind)


def symbolic_distribution(core, exponent, target, data=None):
    """Exact input distribution without enumerating pair states at base B."""
    if not isinstance(exponent, int) or exponent < 1:
        raise ValueError('Requires an integer exponent n>=1.')
    if data is None:
        data = residue_data(core, target)
    scale, base = 1 << exponent, core << exponent
    a, d, w = data['coefficient4'], data['coefficient2'], data['target_weight']
    powers = [1]
    for _ in range(2*exponent+2):
        powers.append(4*powers[-1])
    counts = Counter({0: scale*scale*w, 1: scale*scale*(a+d-w)})
    for t in range(2, exponent+2):
        counts[t] += 3*a*powers[exponent+t-2]
    rules, times = critical_times(core, exponent, data)
    # Every seeded boundary reaches the selected stable cycles in n+2 pair steps.
    for kind in ('Z', 'D', 'A'):
        residues = data['zero'] if kind == 'Z' else data['diagonal']
        for residue in residues:
            for valuation in range((base-1).bit_length()):
                stats = residue_stats(base, core, residue, valuation,
                                      base//2+1 if kind == 'A' else 1)
                add_tree(powers, counts, base, kind, min(valuation, exponent),
                         stats, exponent+2)
    for rule in rules:
        target_time = times[rule.target]
        if target_time is None:
            continue
        arrival = target_time+rule.elapsed
        if rule.family == 'midpoint':
            u = base//2
            add_tree(powers, counts, base, 'D', 0, (1, u, u*u), arrival)
        elif rule.family == 'one':
            add_tree(powers, counts, base, 'Z', 0, (1, 1, 1), arrival)
        elif rule.family == 'zero_block':
            depth = exponent-rule.parameter
            step = core << depth
            add_tree(powers, counts, base, 'Z', depth,
                     moments(step, 2*step, 1, base-1), arrival)
            if depth:
                add_tree(powers, counts, base, 'Z', 0,
                         moments(step+1, 2*step, 1, base-1), arrival)
            else:
                for valuation in range(1, (base-1).bit_length()):
                    add_tree(powers, counts, base, 'Z', min(valuation, exponent),
                             residue_stats(base, core, 1, valuation), arrival)
        elif rule.family == 'diagonal_residues':
            for residue in ((core-1)//2, (core+1)//2):
                for kind in ('D', 'A'):
                    for valuation in range((base-1).bit_length()):
                        add_tree(powers, counts, base, kind, min(valuation, exponent),
                                 residue_stats(base, core, residue, valuation,
                                               base//2+1 if kind == 'A' else 1), arrival)
        elif rule.family == 'diagonal_valuation':
            valuation = exponent-rule.parameter-2
            step = core << valuation
            for kind in ('D', 'A'):
                add_tree(powers, counts, base, kind, valuation,
                         moments(step, 2*step, base//2+1 if kind == 'A' else 1, base-1),
                         arrival)
        elif rule.family == 'zero_residue':
            for valuation in range((base-1).bit_length()):
                add_tree(powers, counts, base, 'Z', min(valuation, exponent),
                         residue_stats(base, core, (core+1)//2, valuation), arrival)
    counts = {t: value for t, value in sorted(counts.items()) if value}
    result = summarize(base, data['target_count'], counts)
    result.update(core=core, exponent=exponent,
                  critical_pair_times={f'{kind}({u})': value
                                       for (kind, u), value in times.items()})
    return result
