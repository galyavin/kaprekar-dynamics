"""Finite leading access through a proved high-region neighborhood of K.

Provides the transferred five-class joint kernel and guarded actual capture.
This extends the solved source, not the global entry theorem.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
from pathlib import Path
import json
from nine_entry import CORE_FACETS, H, C, inverse_C, evaluate, in_core_radius
from odd_gap import subtract_pair
from nine_core import exact_arrival as core_arrival

ROOT = Path(__file__).resolve().parent
ORDER = tuple(CORE_FACETS[i] for i in (4, 5, 6, 7))
OTHER = tuple(g for i, g in enumerate(CORE_FACETS) if i not in (4, 5, 6, 7))
UPPER, LOWER = (-1, 0, 0, 0, 1), (4, 0, 0, 0, -3)


def primitive(row):
    divisor = gcd(*row)
    return tuple(v//divisor for v in row) if divisor else row


COLLAR = tuple(sorted(set((*ORDER, UPPER, LOWER,
    *(primitive(tuple(2*a+b for a, b in zip(g, h))) for g in OTHER for h in (UPPER, LOWER))))))
OVERLAP = tuple(sorted(set((*CORE_FACETS, LOWER))))
CORE_MASS = Q(1973075435497331, 24981193359360000)


def in_collar(r, strict=True):
    values = [evaluate(g, (*r, 1)) for g in COLLAR]
    return min(values) > 0 if strict else min(values) >= 0


def leading_target(x):
    x = tuple(map(Q, x))
    if len(x) != 4 or not 1 > x[0] > x[1] > x[2] > x[3] > 0:
        raise ValueError('Require strict normalized original differences.')
    r = tuple(sorted((abs(2*v-1) for v in x), reverse=True))
    if in_core_radius(r, strict=True): return C(r), 1
    if not in_collar(r): raise ValueError('Outside the solved core/collar source.')
    for _ in range(6):
        raw = (r[0], *(abs(a+b-1) for a, b in zip(r, r[1:])))
        if min(raw) <= 0 or len(set(raw)) != 4:
            raise ValueError('Leading prefix lies on a sorting boundary.')
        r = H(r)
    if not in_core_radius(r, strict=True):
        raise ValueError('Leading prefix reaches a core boundary.')
    return C(r), 7


@lru_cache(None)
def route_table():
    from nine_core_checks import route_domains
    records = json.loads((ROOT/'nine_kernel_results.json').read_text())['kernels']
    return route_domains(), {(tuple(r['route']), r['base_residue']): r for r in records}


def limit_kernel(x, residue):
    if residue not in range(5): raise ValueError('Base residue must be 0,...,4.')
    y, prefix = leading_target(x); r = inverse_C(y)
    domains, table = route_table()
    routes = [tuple(s['route']) for s in domains if min(evaluate(g, (*r, 1)) for g in s['domain']) > 0]
    if len(routes) != 1: raise ValueError('Leading core route lies on an exit tie.')
    record = table[routes[0], residue]
    atoms = [{'weight': Q(a['weight']), 'time': evaluate(tuple(map(Q, a['clock'])), (*y, 1)),
              'terminal': tuple(evaluate(tuple(map(Q, g)), (*y, 1)) for g in a['terminal'])}
             for a in record['atoms']]
    return {'source': tuple(map(Q, x)), 'core_source': y, 'prefix_steps': prefix,
            'route': routes[0], 'atoms': atoms}


def exact_arrival(B, q):
    q = tuple(q)
    if not isinstance(B, int) or B < 1000000 or any(not isinstance(v, int) for v in q):
        raise ValueError('Require integer differences and B>=1000000.')
    _, prefix = leading_target(tuple(Q(v, B) for v in q))
    target = q
    for _ in range(prefix): target = subtract_pair(B, target)
    tail = core_arrival(B, target, first_image=False)
    if tail['pair_time'] <= 0: raise ArithmeticError('Core target must be transient.')
    # A predecessor of a transient state cannot be cyclic; no endpoint-only
    # assumption is needed to exclude the entire finite incoming prefix.
    return {'source': q, 'prefix_steps': prefix, 'core_source': target, 'core_result': tail,
            'pair_time': prefix+tail['pair_time'], 'actual_time': prefix+tail['actual_time'],
            'period': tail['period'], 'first_cycle_state': tail['first_cycle_state'],
            'route': tail['route']}
