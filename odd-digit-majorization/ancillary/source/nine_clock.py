"""Numerical quantities for the nine-digit clock certificate and capture bound."""
from math import comb
from nine_clock_reduction import dispersion, terminal_residual, zero_preparation_bound

HORIZON = 9
THRESHOLD = 100
BASE_MINIMUM = 40


def exceptional_count_bound(base):
    if not isinstance(base,int) or base < BASE_MINIMUM:
        raise ValueError('Require integer B >= 40.')
    return 5*(base+1)*comb(THRESHOLD+4,4)


def pair_clock_bound(base):
    return (zero_preparation_bound(base)+4*HORIZON*base
            +(8*HORIZON+1)*exceptional_count_bound(base)+HORIZON)


def number_clock_bound(base):
    return pair_clock_bound(base)+1
