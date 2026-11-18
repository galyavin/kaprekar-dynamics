"""Exact recurrent words and drift on the B0 facet."""
from math import gcd
from fractions import Fraction as Q

from nine_entry import CORE_FACETS,ORDER,compose,evaluate
from odd_affine import word,progression

INVERSE_C=((2,0,0,0,-1),(-2,2,0,0,1),(2,-2,2,0,-1),
           (-2,2,-2,2,1),(0,0,0,0,1))
B0=(0,1,1,0,-1)
FACETS=tuple(tuple(x//gcd(*r) for x in r)
             for r in (compose(row,INVERSE_C) for row in CORE_FACETS))
NON_B=tuple(r for r in FACETS if r!=B0)
PHASES=(-4,-3,-2,-1)
DELTA_B=(-6,-1,1,-2)
LOW_ORDER=(4,0,5,6,1,7,2,3,8)
THREE_ORDERS={'high':(ORDER,ORDER,ORDER),'low':(LOW_ORDER,ORDER,ORDER)}


def phase(B,q):
    return q[1]+q[2]-B


def phase_next(e):
    return e-3 if e>=-1 else -e-5


def three_formula(B,q):
    u,v,w,z=q
    if phase(B,q)>=-1:return u-3,v-2,w-1,z-1
    return u-3,B-w-3,B-v-2,z-1


def six_orders(e):
    if e not in PHASES:raise ValueError('Require one of the four recurrent phases.')
    return THREE_ORDERS['high' if e>=-1 else 'low']+THREE_ORDERS[
        'high' if phase_next(e)>=-1 else 'low']


def guard_domain(key):
    """Nonnegative affine rows in (u,v,w,z,B,1)."""
    if key not in THREE_ORDERS:raise ValueError('Unknown branch.')
    return (*[(*r,-10000) for r in NON_B],(0,0,0,0,1,-100000),
            (*B0,512),(*[-x for x in B0],512),
            (*B0,1) if key=='high' else (*[-x for x in B0],-2))


def absorb(B,q):
    """At most 513 actual steps to a recurrent integer phase.

    The 14000 margin guarantees every three-step guard throughout absorption.
    It is a convenient sufficient domain, not a sharp boundary threshold.
    """
    q=tuple(q)
    if B<100000 or abs(phase(B,q))>512 or min(
            evaluate(r,(*q,B)) for r in NON_B)<14000:
        raise ValueError('Outside the certified B0 absorption domain.')
    elapsed=0
    while phase(B,q) not in PHASES:
        q=three_formula(B,q);elapsed+=3
        if elapsed>513:raise ArithmeticError('Phase absorption bound failed.')
    return {'elapsed':elapsed,'target':q,'phase':phase(B,q)}


def facet_exit(B,q):
    """Maximal exact repetitions of the recurrent six-step word."""
    q=tuple(q);e=phase(B,q)
    expected=six_orders(e)
    record=word(B,q,6)
    if record['orders']!=expected:
        raise ValueError('Outside the specified B0 six-step word.')
    count,delta=progression(B,q,record)
    if delta!=DELTA_B:raise ArithmeticError('B0 phase displacement failed.')
    return {'elapsed':6*count,'repetitions':count,
            'target':tuple(a+count*d for a,d in zip(q,delta)),'phase':e}


def leading_exit(x):
    """Leading facet time and location; x must lie strictly on B0."""
    x=tuple(map(Q,x))
    if evaluate(B0,(*x,1))!=0 or min(evaluate(r,(*x,1)) for r in NON_B)<=0:
        raise ValueError('Require the relative interior of B0 in C(K).')
    choices=[(evaluate(r,(*x,1))/(-evaluate(r[:4],DELTA_B)),r)
             for r in NON_B if evaluate(r[:4],DELTA_B)<0]
    t=min(t for t,r in choices)
    return {'elapsed_over_B':6*t,'target':tuple(a+t*d for a,d in zip(x,DELTA_B)),
            'faces':tuple(r for s,r in choices if s==t)}
