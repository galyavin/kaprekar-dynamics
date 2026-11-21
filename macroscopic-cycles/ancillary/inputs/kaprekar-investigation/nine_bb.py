"""Exact nine-digit B0/B1 corner transport, retaining integer phases.

All initial integer phases absorb; actual words require their checked guards.
The twelve-step drift is shared by three phase cycles. Terminal capture is a
separate obligation unless this exit is composed with nine_bbb.exact_arrival.
"""
from fractions import Fraction as Q
from nine_entry import ORDER,evaluate
from nine_b_facet import FACETS,B0,LOW_ORDER
from odd_affine import word,progression

B1=(1,1,0,-2,-1)
B2=(2,-2,2,-1,-1)
NON_BB=tuple(r for r in FACETS if r not in (B0,B1))
DELTA_BB=(-12,-2,2,-7)
THREE_ORDERS=((ORDER,ORDER,ORDER),(LOW_ORDER,ORDER,ORDER),
              (LOW_ORDER,LOW_ORDER,ORDER),(ORDER,LOW_ORDER,ORDER))
PHASE_CYCLES=(((-3,-4),(-2,-5),(-3,-3),(-2,-4)),
              ((-1,-2),(-4,-5),(-1,-5),(-4,-2)),
              ((-1,-3),(-4,-4),(-1,-4),(-4,-3)))
PHASES=tuple(sorted(p for cycle in PHASE_CYCLES for p in cycle))


def decode(B,q):
    u,v,w,z=q
    return u,z,(v+w-B,u+v-2*z-B)


def encode(B,U,Z,phase):
    p,r=phase
    return U,2*Z+B-U+r,U-2*Z-r+p,Z


def phase_step(phase):
    p,r=phase;h=r+max(-p-1,0)
    return (abs(p+1)-4,abs(h+2)-5),min(h+1,-1)


def phase_entry_bound(phase):
    """A coarse bound in three-step returns from the two section folds."""
    p,r=phase
    return 17*abs(p)+2*abs(r)+73


def branch(phase):
    p,r=phase
    if p>=-1:return 0 if r>=-2 else 3
    return 1 if r>=p-1 else 2


def phase_domain(key):
    return (((0,1,2),(1,0,1)),((-1,0,-2),(-1,1,1)),
            ((-1,0,-2),(1,-1,-2)),((0,-1,-3),(1,0,1)))[key]


def twelve_orders(p):
    if p not in PHASES:raise ValueError('Require a recurrent BB phase.')
    orders=()
    for _ in range(4):orders+=THREE_ORDERS[branch(p)];p,_=phase_step(p)
    return orders


def absorb(B,q,max_returns=100000):
    q=tuple(q);U,Z,p=decode(B,q);count=0;bound=phase_entry_bound(p)
    while p not in PHASES:
        rec=word(B,q,3)
        if rec['orders']!=THREE_ORDERS[branch(p)]:raise ValueError('Outside the BB corner word.')
        p,dz=phase_step(p);U-=3;Z+=dz;q=encode(B,U,Z,p);count+=1
        if count>bound:raise ArithmeticError('BB phase absorption bound failed.')
        if count>=max_returns:raise ValueError('Phase evaluation cap, not a nonarrival conclusion.')
    return {'elapsed':3*count,'target':q,'phase':p}


def corner_exit(B,q):
    q=tuple(q);_,_,p=decode(B,q);orders=twelve_orders(p);record=word(B,q,12)
    if record['orders']!=orders:raise ValueError('Outside the recurrent BB corner word.')
    count,delta=progression(B,q,record)
    if delta!=DELTA_BB:raise ArithmeticError('BB corner displacement failed.')
    return {'elapsed':12*count,'repetitions':count,
            'target':tuple(a+count*d for a,d in zip(q,delta)),'phase':p}


def leading_exit(x):
    x=tuple(map(Q,x))
    if any(evaluate(row,(*x,1)) for row in (B0,B1)) or min(
            evaluate(row,(*x,1)) for row in NON_BB)<=0:
        raise ValueError('Require the relative interior of B0=B1=0.')
    choices=[(evaluate(row,(*x,1))/(-evaluate(row[:4],DELTA_BB)),row)
             for row in NON_BB if evaluate(row[:4],DELTA_BB)<0]
    t=min(t for t,row in choices)
    return {'elapsed_over_B':12*t,'target':tuple(a+t*d for a,d in zip(x,DELTA_BB)),
            'faces':tuple(row for s,row in choices if s==t)}
