"""Exact A2/B1 corner transport; subsequent capture remains separate."""
from fractions import Fraction as Q
from nine_entry import ORDER,evaluate
from nine_b_facet import FACETS,LOW_ORDER
from nine_a_facet import A2,branch as a_branch,phase_next,THREE_ORDERS as A_ORDERS
from nine_bb import B1
from odd_affine import word,progression

A1=(1,1,-1,0,-1)
NON_AB=tuple(row for row in FACETS if row not in (A2,B1))
DELTA_AB=(-12,-24,-12,-18)
PHASE_CYCLES=(((-2,-5),(4,0),(0,-3),(2,-4),(-2,-3),(4,-2),(0,-5),(2,-2)),
              ((-2,-4),(4,-1),(0,-4),(2,-3)),
              ((-1,-5),(3,-1),(-1,-4),(3,-2)))
PHASES=tuple(sorted(p for cycle in PHASE_CYCLES for p in cycle))
KEYS=tuple((a,b) for a in ('low','middle','high') for b in ('high','low'))


def decode(B,q):
    u,v,w,z=q
    return u,w,(2*u-v-B,u+v-2*z-B)


def encode(B,U,W,phase):
    p,r=phase
    if (3*U-p-r)%2:raise ValueError('Integer AB chart requires U=p+r modulo two.')
    return U,2*U-B-p,W,(3*U-2*B-p-r)//2


def phase_step(phase):
    p,r=phase
    return (phase_next(p),abs(r+2)-5+max(-p,0)),(-p-1 if p<2 else -3),min(p,0)-1


def phase_entry_bound(phase):
    p,r=phase
    return 9*abs(p)+abs(r)+70


def branch(phase):
    p,r=phase
    return a_branch(p),'high' if r>=-2 else 'low'


def phase_domain(key):
    a,b=key
    pg={'low':((-1,0,-1),),'middle':((1,0,0),(-1,0,1)),
        'high':((1,0,-2),)}[a]
    return (*pg,(0,1,2) if b=='high' else (0,-1,-3))


def three_orders(key):
    a,b=key
    return ORDER,ORDER if b=='high' else LOW_ORDER,A_ORDERS[a][-1]


def return_orders(p):
    if p not in PHASES:raise ValueError('Require a recurrent A-B phase.')
    orders=()
    for _ in range(8):orders+=three_orders(branch(p));p,_,_=phase_step(p)
    return orders


def absorb(B,q,max_returns=100000):
    q=tuple(q);U,W,p=decode(B,q);count=0;bound=phase_entry_bound(p)
    while p not in PHASES:
        record=word(B,q,3)
        if record['orders']!=three_orders(branch(p)):raise ValueError('Outside the AB corner word.')
        p,du,dw=phase_step(p);U+=du;W+=dw;q=encode(B,U,W,p);count+=1
        if count>bound:raise ArithmeticError('AB phase absorption bound failed.')
        if count>=max_returns:raise ValueError('Phase evaluation cap, not a nonarrival conclusion.')
    return {'elapsed':3*count,'target':q,'phase':p}


def corner_exit(B,q):
    q=tuple(q);_,_,p=decode(B,q);orders=return_orders(p);record=word(B,q,24)
    if record['orders']!=orders:raise ValueError('Outside the recurrent AB corner word.')
    count,delta=progression(B,q,record)
    if delta!=DELTA_AB:raise ArithmeticError('AB corner displacement failed.')
    return {'elapsed':24*count,'repetitions':count,
            'target':tuple(a+count*d for a,d in zip(q,delta)),'phase':p}


def leading_exit(x):
    x=tuple(map(Q,x))
    if any(evaluate(row,(*x,1)) for row in (A2,B1)) or min(
            evaluate(row,(*x,1)) for row in NON_AB)<=0:
        raise ValueError('Require the relative interior of A2=B1=0.')
    t=evaluate(A1,(*x,1))/24
    return {'elapsed_over_B':24*t,'target':tuple(a+t*d for a,d in zip(x,DELTA_AB)),
            'faces':(A1,)}
