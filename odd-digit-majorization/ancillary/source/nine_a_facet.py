"""Exact nine-digit A2-facet segments; no terminal-arrival claim.

Together with nine_b_facet, this covers the first facet drift after every
strict core bulk exit. All elapsed values are timed paths, not first capture.
"""
from fractions import Fraction as Q
from nine_entry import ORDER,evaluate
from nine_b_facet import FACETS
from odd_affine import word,progression

A2=(2,-1,0,0,-1)
NON_A=tuple(r for r in FACETS if r!=A2)
PHASES=(-2,-1,0,2,3,4)
DELTA_A=(-6,-12,-6,-4)
THREE_ORDERS={
    'low':(ORDER,ORDER,(4,5,0,1,6,2,7,8,3)),
    'middle':(ORDER,ORDER,(4,0,5,1,6,2,7,8,3)),
    'high':(ORDER,ORDER,ORDER),
}


def phase(B,q):return 2*q[0]-q[1]-B


def branch(e):return 'low' if e<0 else 'middle' if e<2 else 'high'


def phase_next(e):return 2-e if e<0 else 2-2*e if e<2 else e-4


def three_formula(B,q):
    u,v,w,z=q;e=phase(B,q)
    if e<0:return -u+v+B-1,v-4,w+e-1,z-1
    if e<2:return -u+v+B-1,2*u-B-4,w-1,z-1
    return u-3,v-2,w-1,z-1


def twelve_orders(e):
    if e not in PHASES:raise ValueError('Require a recurrent A2 phase.')
    orders=()
    for _ in range(4):orders+=THREE_ORDERS[branch(e)];e=phase_next(e)
    return orders


def guard_domain(key):
    common=(*[(*r,-10000) for r in NON_A],(0,0,0,0,1,-100000),
            (*A2,512),(*[-x for x in A2],514))
    if key=='low':return (*common,(*[-x for x in A2],-1))
    if key=='middle':return (*common,(*A2,0),(*[-x for x in A2],1))
    if key=='high':return (*common,(*A2,-2))
    raise ValueError('Unknown A2 branch.')


def absorb(B,q):
    q=tuple(q)
    if B<100000 or abs(phase(B,q))>512 or min(evaluate(r,(*q,B)) for r in NON_A)<14000:
        raise ValueError('Outside the certified A2 absorption domain.')
    elapsed=0
    while phase(B,q) not in PHASES:
        q=three_formula(B,q);elapsed+=3
        if elapsed>390:raise ArithmeticError('A2 phase absorption bound failed.')
    return {'elapsed':elapsed,'target':q,'phase':phase(B,q)}


def facet_exit(B,q):
    q=tuple(q);e=phase(B,q);expected=twelve_orders(e);record=word(B,q,12)
    if record['orders']!=expected:raise ValueError('Outside the A2 twelve-step word.')
    count,delta=progression(B,q,record)
    if delta!=DELTA_A:raise ArithmeticError('A2 phase displacement failed.')
    return {'elapsed':12*count,'repetitions':count,
            'target':tuple(a+count*d for a,d in zip(q,delta)),'phase':e}


def leading_exit(x):
    x=tuple(map(Q,x))
    if evaluate(A2,(*x,1))!=0 or min(evaluate(r,(*x,1)) for r in NON_A)<=0:
        raise ValueError('Require the relative interior of A2 in C(K).')
    choices=[(evaluate(r,(*x,1))/(-evaluate(r[:4],DELTA_A)),r)
             for r in NON_A if evaluate(r[:4],DELTA_A)<0]
    t=min(t for t,r in choices)
    return {'elapsed_over_B':12*t,'target':tuple(a+t*d for a,d in zip(x,DELTA_A)),
            'faces':tuple(r for s,r in choices if s==t)}
