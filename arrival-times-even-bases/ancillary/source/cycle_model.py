from collections import Counter
from fractions import Fraction
from investigate import transition
from predecessor_model import predecessors, valid
from doubling_checks import multiplicity


def number_from_pair(base, pair):
    u, v = pair
    return (base-1)*((base*base+base+1)*u+base*v)


def summarize(base, target_count, pair_counts):
    """Convert weighted pair-arrival counts to actual-input counts."""
    counts = Counter({t+1: w for t, w in pair_counts.items() if w})
    counts[1] -= target_count
    counts[0] += target_count
    counts = dict(sorted((t, w) for t, w in counts.items() if w))
    assert all(w > 0 for w in counts.values())
    total = sum(counts.values())
    first = sum(t*w for t, w in counts.items())
    second = sum(t*t*w for t, w in counts.items())
    assert 0 <= total <= base**4-base
    return {
        'base': base, 'target_pair_count': target_count,
        'pair_arrival_counts': dict(sorted(pair_counts.items())),
        'input_arrival_counts': counts,
        'basin_input_count': total, 'complement_input_count': base**4-base-total,
        'first_moment_sum': first, 'second_moment_sum': second,
        'conditional_mean': str(Fraction(first, total)) if total else None,
        'conditional_second_moment': str(Fraction(second, total)) if total else None,
        'maximum_finite_input_time': max(counts) if counts else None,
    }


def weighted_cycle_layers(base, target, retain_states=False):
    if not isinstance(base, int) or base < 2:
        raise ValueError('Base must be an integer at least 2.')
    target = frozenset(target)
    if any(not valid(base, pair) for pair in target):
        raise ValueError('Invalid target pair.')
    if {transition(base, pair) for pair in target} != target:
        raise ValueError('Target must be an invariant union of cycles.')
    seen, frontier, layers = set(target), set(target), []
    pair_counts = {}
    while frontier:
        weight = sum(multiplicity(base, pair) for pair in frontier)
        item = {'pair_count': len(frontier), 'input_weight': weight}
        if retain_states:
            item['states'] = sorted(frontier)
        pair_counts[len(layers)] = weight
        layers.append(item)
        following = set()
        for pair in frontier:
            following.update(predecessors(base, pair))
        frontier = following - seen
        seen.update(frontier)
    result = summarize(base, len(target), pair_counts)
    result.update(layers=layers, basin_pair_count=len(seen))
    return result
