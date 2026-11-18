"""Certified seven-digit edge clocks, conditional on an explicit entry strip.

This is not a global seven-digit arrival law. The small offset graph retains
two different average speeds that a continuous state alone cannot determine.
"""
from fractions import Fraction as Q

from odd_gap import subtract_pair, valid


def offset_return(parity, state):
    """Base-free three-step edge return, certified by seven_boundary_audit.py."""
    m, s = state
    if parity not in (0, 1) or m % 2 != parity:
        raise ValueError('Offset parity must match the base.')
    base, u = 100000+parity, 80000
    q = (u, (base+m)//2, base-u+s)
    for _ in range(3):
        q = subtract_pair(base, q)
    return (2*q[1]-base, q[0]+q[2]-base), q[0]-u


def phase_data(parity, state):
    seen, path, shifts = {}, [], []
    while state not in seen:
        seen[state] = len(path); path.append(state)
        state, shift = offset_return(parity, state)
        shifts.append(shift)
        if len(path) > 100:
            raise ValueError('State is outside the certified offset graph.')
    cut = seen[state]
    delta = sum(shifts[cut:])
    return {'transient_blocks': cut, 'transient_u_shift': sum(shifts[:cut]),
            'cycle_state': path[cut], 'cycle_blocks': len(path)-cut,
            'cycle_u_shift': delta, 'time_per_unit_u': Q(3*(len(path)-cut), -delta)}


def edge_domain(base, q):
    if not valid(base, q) or len(q) != 3 or base < 100000:
        return False
    u, v, w = q
    return (abs(2*v-base) <= 20 and abs(u+w-base) <= 20
            and 4*u >= 3*base+4000 and 6*u <= 5*base-6000)


def edge_clock(base, q):
    """Pair arrival estimate with audited absolute error at most 200.

    The rate is either 3/2 or 2, selected by the integer offset graph.
    """
    if not edge_domain(base, q):
        raise ValueError('Outside the certified edge-entry domain.')
    u, v, w = q
    data = phase_data(base % 2, (2*v-base, u+w-base))
    return data['time_per_unit_u']*(Q(u)-Q(3*base, 4))


def edge_arrival(base, q):
    """Exact pair-cycle arrival in the certified edge domain.

    Skip complete offset cycles, then resolve first capture by direct digits.
    The skipped region is disjoint from the audited terminal cycles.
    """
    if not edge_domain(base, q):
        raise ValueError('Outside the certified edge-entry domain.')
    u, v, w = q
    data = phase_data(base % 2, (2*v-base, u+w-base))
    elapsed = 3*data['transient_blocks']
    for _ in range(elapsed):
        q = subtract_pair(base, q)
    decrement = -data['cycle_u_shift']
    threshold = (3*base+3)//4+256
    repeats = max(0, (q[0]-threshold)//decrement+1)
    q = (q[0]-repeats*decrement, q[1], q[2]+repeats*decrement)
    elapsed += repeats*3*data['cycle_blocks']
    seen = {}
    while q not in seen:
        seen[q] = len(seen)
        q = subtract_pair(base, q)
        assert len(seen) < 5000
    return elapsed+seen[q], len(seen)-seen[q]
