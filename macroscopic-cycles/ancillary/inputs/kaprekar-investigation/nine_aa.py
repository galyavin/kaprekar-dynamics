"""Phase-dependent A1/A2 corner transport; no terminal-capture claim."""
from fractions import Fraction as Q
from nine_entry import ORDER,evaluate
from nine_b_facet import FACETS
from nine_a_facet import A2,branch as a_branch,THREE_ORDERS as A_ORDERS
from nine_ab import A1
from odd_affine import word,progression

A0=(1,0,0,1,-1)
NON_AA=tuple(row for row in FACETS if row not in (A1,A2))
PHASE_CYCLES=(((-2,-2),(-2,4),(4,4),(0,0),(-2,2),(4,2),(0,-2),
               (0,4),(2,0),(0,2),(2,-2),(2,4),(-2,0),(2,2)),
              ((-1,-2),(-1,4),(3,2)))
PHASES=tuple(sorted(p for cycle in PHASE_CYCLES for p in cycle))
DELTAS=((-10,-20,-30,-20),(-2,-4,-6,-3))
KEYS=tuple((a,b) for a in ('low','middle','high') for b in ('low','middle','high'))


def decode(B,q):
    u,v,w,z=q
    return u,z,(2*u-v-B,u+v-w-B)


def encode(B,U,Z,phase):
    p,r=phase
    return U,2*U-B-p,3*U-2*B-p-r,Z


def branch(phase):
    p,r=phase;a=a_branch(r)
    b=a_branch(p) if r>=2 else 'low' if p-r<=-3 else 'middle' if p-r<0 else 'high'
    return a,b


def phase_domain(key):
    a,b=key
    rd={'low':((0,-1,-1),),'middle':((0,1,0),(0,-1,1)),
        'high':((0,1,-2),)}[a]
    if a=='high':pd={'low':((-1,0,-1),),'middle':((1,0,0),(-1,0,1)),
                     'high':((1,0,-2),)}[b]
    else:pd={'low':((-1,1,-3),),'middle':((1,-1,2),(-1,1,-1)),
             'high':((1,-1,0),)}[b]
    return (*rd,*pd)


def phase_matrix(key):
    a,b=key
    if a=='high':
        pn={'low':(-1,0,2),'middle':(-2,0,2),'high':(1,0,-4)}[b]
        rn=(-2 if b=='low' else 0,1,-4)
        du=(0,0,-3) if b=='high' else (-1,0,-1);dz=(0,0,-1)
    else:
        pn={'low':(-1,1,0),'middle':(-2,2,-2),'high':(1,-1,-2)}[b]
        rn=(-2,1 if a=='low' else 0,-2) if b=='low' else (0,-1 if a=='low' else -2,2)
        du=(0,-1,-1) if b=='high' else (-1,0,-1)
        dz=(0,-1 if a=='low' else 1,-3)
    return pn,rn,du,dz


def phase_step(phase):
    return tuple(evaluate(row,(*phase,1)) for row in phase_matrix(branch(phase)))


def phase_entry_bound(phase):
    p,r,_,_=phase_step(phase)
    return 1+2*max(0,max(p,r)-4)+9


def three_orders(key):
    a,b=key
    return ORDER,A_ORDERS[a][-1],A_ORDERS[b][-1]


def return_spec(p):
    if p not in PHASES:raise ValueError('Require a recurrent AA phase.')
    kind=next(i for i,c in enumerate(PHASE_CYCLES) if p in c);orders=()
    for _ in PHASE_CYCLES[kind]:
        orders+=three_orders(branch(p));pn,rn,_,_=phase_step(p);p=(pn,rn)
    return orders,DELTAS[kind]


def absorb(B,q,max_returns=100000):
    q=tuple(q);U,Z,p=decode(B,q);count=0;bound=phase_entry_bound(p)
    while p not in PHASES:
        record=word(B,q,3)
        if record['orders']!=three_orders(branch(p)):raise ValueError('Outside the AA corner word.')
        pn,rn,du,dz=phase_step(p);p=(pn,rn);U+=du;Z+=dz;q=encode(B,U,Z,p);count+=1
        if count>bound:raise ArithmeticError('AA phase absorption bound failed.')
        if count>=max_returns:raise ValueError('Phase evaluation cap, not a nonarrival conclusion.')
    return {'elapsed':3*count,'target':q,'phase':p}


def corner_exit(B,q):
    q=tuple(q);_,_,p=decode(B,q);orders,expected=return_spec(p);record=word(B,q,len(orders))
    if record['orders']!=orders:raise ValueError('Outside the recurrent AA corner word.')
    count,delta=progression(B,q,record)
    if delta!=expected:raise ArithmeticError('AA corner displacement failed.')
    return {'elapsed':len(orders)*count,'repetitions':count,'length':len(orders),
            'target':tuple(a+count*d for a,d in zip(q,delta)),'phase':p}


def leading_exit(x,kind):
    x=tuple(map(Q,x))
    if kind not in (0,1):raise ValueError('Require phase kind zero or one.')
    if any(evaluate(row,(*x,1)) for row in (A1,A2)) or min(
            evaluate(row,(*x,1)) for row in NON_AA)<=0:
        raise ValueError('Require the relative interior of A1=A2=0.')
    delta=DELTAS[kind];length=3*len(PHASE_CYCLES[kind])
    t=evaluate(A0,(*x,1))/(-evaluate(A0[:4],delta))
    return {'elapsed_over_B':length*t,'target':tuple(a+t*d for a,d in zip(x,delta)),
            'faces':(A0,)}
