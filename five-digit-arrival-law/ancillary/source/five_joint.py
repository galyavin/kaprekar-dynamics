"""Five-digit joint limiting laws and exact cycle labels after core entry.

See five_digit_limit.md. Cycle labels are Prichett's T1--T5 labels; T4 and
T5 are different four-cycles in different base congruence classes.
"""
from fractions import Fraction

from five_clock import core_region
from five_digit import valid, transition
from five_limit_law import cdf, moment_bounds


def cycle_weights(base_residue):
    """Limiting cycle probabilities along B == base_residue (mod 3)."""
    if base_residue not in (0, 1, 2):
        raise ValueError('Require residue 0, 1, or 2.')
    return ({'T1': Fraction(1)},
            {'T3': Fraction(2, 3), 'T4': Fraction(1, 3)},
            {'T5': Fraction(1)})[base_residue]


def joint_cdf(t, label, base_residue):
    """Limit of P(T_B/B <= t, terminal cycle == label)."""
    if label not in ('T1', 'T2', 'T3', 'T4', 'T5'):
        raise ValueError('Unknown classical cycle label.')
    weight = cycle_weights(base_residue).get(label, Fraction(0))
    return weight * cdf(t) if weight else Fraction(0)


def joint_moment_bounds(order, label, base_residue, cutoff=300):
    """Certified limiting mixed moment E[L^order 1{C=label}]."""
    if label not in ('T1', 'T2', 'T3', 'T4', 'T5'):
        raise ValueError('Unknown classical cycle label.')
    weight = cycle_weights(base_residue).get(label, Fraction(0))
    if not weight:
        return Fraction(0), Fraction(0)
    return tuple(weight * x for x in moment_bounds(order, cutoff))


def rare_pairs(base):
    if not isinstance(base, int) or base < 30:
        raise ValueError('Exact rare-basin theorem requires integer B >= 30.')
    if base % 2:
        return frozenset()
    m = base // 2
    return frozenset(((m, m - 1), (m + 1, 0), (m, 0)))


def rare_arrival_counts(base):
    """Exact ordered-input histogram for the entire T2 basin, including t=0."""
    if not rare_pairs(base):
        return {}
    return {0: 2, 1: 35 * base**2 - 75 * base - 12,
            2: 5 * base**2 - 5 * base}


def core_cycle_label(base, pair):
    if not core_region(base, pair):
        raise ValueError('Require a pair in the certified five-digit core.')
    if base % 3 == 0:
        return 'T1'
    if base % 3 == 2:
        return 'T5'
    return 'T4' if sum(pair) % 3 == 0 else 'T3'


def terminal_label(base, pair):
    """Exact label via core entry, without following the O(B) core drift.

    The existing global entry theorem bounds this loop by O(sqrt(B)).
    This function does not compute an exact arrival time.
    """
    if base < 30 or not valid(base, pair):
        raise ValueError('Require integer B >= 30 and a valid pair.')
    rare = rare_pairs(base)
    if pair in rare:
        return 'T2'
    while not core_region(base, pair):
        pair = transition(base, pair)
    return core_cycle_label(base, pair)
