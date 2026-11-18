"""Exact phase contraction and first capture on the triple-B terminal strip."""
from nine_entry import ORDER,evaluate
from nine_b_facet import LOW_ORDER
from odd_affine import word
from odd_gap import valid

CANONICAL=(-3,2,0)
# B=5b+r; entry V=3b+256+s. Each pair is (further pair time, period).
TERMINALS=(
    ((780,21),(793,21),(778,21)),
    ((777,9),(786,1),(784,9)),
    ((770,14),(796,2),(777,14)),
    ((784,10),(776,2),(769,10)),
    ((767,23),(770,23),(765,23)),
)
# First cyclic pair minus (4b,3b,2b,b), in the same base/entry order.
CAPTURE_OFFSETS=(
    ((2,1,3,2),(0,0,0,-1),(2,2,3,2)),
    ((1,1,1,1),(1,1,1,0),(1,1,1,1)),
    ((2,2,2,1),(2,2,0,-1),(2,2,2,1)),
    ((3,3,1,0),(4,3,4,2),(3,1,2,1)),
    ((4,2,3,1),(5,3,4,3),(3,2,3,2)),
)


def phase_step(p):
    a,b,c=p;l=a-c+1
    return (-2*a+3*c-4,a-2*c-min(b,-1),-b-2*a+2*c-2),l


def coordinates(p):
    a,b,c=p
    return b+1,a-2*c+3,2*a+2*b-c+4


def inverse_coordinates(q):
    x,y,z=q
    if sum(q)%3!=2:raise ValueError('Unrealizable integer strip phase.')
    return (2*z-y-4*x-1)//3,x-1,(z-2*y-2*x+4)//3


def folded_step(q):
    x,y,z=q
    return y-min(x,0)-1,z-1,abs(x)-1


def decode(B,q):
    u,V,w,z=q
    return V,(u-3*V+B,w+V-B,z-2*V+B)


def encode(B,V,p):
    a,b,c=p
    return 3*V-B+a,V,B-V+b,2*V-B+c


def edge_state(B,V):return encode(B,V,CANONICAL)


def strip_step(B,q):
    """One exact step when the two possible terminal-strip orders apply."""
    V,p=decode(B,q);record=word(B,q,1)
    expected=ORDER if p[1]>=-1 else LOW_ORDER
    if record['orders']!=(expected,):raise ValueError('Outside the triple-B strip word.')
    new,delta_V=phase_step(p)
    return encode(B,V+delta_V,new)


def phase_entry_bound(p):
    q=coordinates(p)
    for _ in range(3):q=folded_step(q)
    return 3+6*max(0,max(q)-8)+30


def exact_arrival(B,q,max_phase_steps=100000):
    """First pair-cycle arrival when the checked strip/edge contract applies.

    Arbitrary fixed bounded offsets satisfy this contract eventually on every
    compact subinterval of 3/5<V/B<2/3. Rejection at finite B is not nonarrival.
    """
    q=tuple(q)
    if (not isinstance(B,int) or B<1000000 or len(q)!=4
            or any(not isinstance(v,int) for v in q) or not valid(B,q)):
        raise ValueError('Require a valid nine-slot state and B>=1000000.')
    initial=q;V,p=decode(B,q);bound=phase_entry_bound(p)
    if bound>max_phase_steps:raise ValueError('Phase evaluation cap, not a dynamical obstruction.')
    prefix=0;base_block,base_residue=divmod(B,5)
    while p!=CANONICAL:
        if q[0]<=4*base_block+6:raise ValueError('Missing safe-capture buffer.')
        q=strip_step(B,q);V,p=decode(B,q);prefix+=1
        if prefix>bound:raise ArithmeticError('Certified strip phase entry bound failed.')
    if not (5*V>=3*B+1295 and 3*V<=2*B-768):
        raise ValueError('Outside the buffered nine-step edge interval.')
    residue=(V-3*base_block-256)%3
    count=(V-3*base_block-256-residue)//3
    h,period=TERMINALS[base_residue][residue]
    pair=prefix+9*count+h
    target=edge_state(B,3*base_block+256+residue)
    if pair!=initial[0]-target[0]+h:raise ArithmeticError('Longitudinal clock identity failed.')
    return {'source':initial,'pair_time':pair,'actual_time':pair+1,'period':period,
            'phase_prefix':prefix,'edge_repetitions':count,'terminal_seed':target,
            'first_cycle_state':tuple(k*base_block+c for k,c in zip(
                (4,3,2,1),CAPTURE_OFFSETS[base_residue][residue])),
            'base_residue':base_residue,'entry_residue':residue}
