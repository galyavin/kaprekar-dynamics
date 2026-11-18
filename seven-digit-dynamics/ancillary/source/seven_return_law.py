"""Certified return construction on the seven-digit leading period-three region.

Integer facet transit and rational certificate replay for both source signs.
"""
from fractions import Fraction as Q
from itertools import product

from odd_gap import valid, subtract_pair
from odd_affine import word, progression
from seven_facet import phase_iterate, tetra_weight
from seven_boundary import edge_domain, edge_arrival, phase_data

# Homogeneous affine rows in (u,v,w,B,1), all interpreted as >=0.
RETURN_DOMAIN = (
    (0,0,0,0,1),(0,0,0,1,-100000),(-1,0,0,1,-1000),
    (0,0,1,0,-1000),(0,2,0,-1,6),
    (1,0,1,-1,-16),(1,-2,1,0,-16),
    (2,0,-2,-1,-16),(2,-2,-2,0,-16),
    (-1,-1,2,1,-8),(-2,2,-1,1,-8),(0,-1,-1,1,-8),
)
SOURCE_COMMON = (
    (0,0,0,1,-100000),(-1,0,0,1,-1000),(0,0,1,0,-1000),
    (0,2,0,-1,0),(1,-2,1,0,-1000),(2,-2,-2,0,-1000),
    (-1,-1,2,1,-1000),(-2,2,-1,1,-1000),(0,-1,-1,1,-1000),
)
SOURCE_POSITIVE = (*SOURCE_COMMON,(1,0,-3,0,-8000),(-3,0,9,1,-24000))
SOURCE_NEGATIVE = (*SOURCE_COMMON,(-1,0,3,0,-8000),(6,0,-18,1,-24000))
NEGATIVE_EXIT_DOMAIN = (
    (0,0,0,0,1),(0,0,0,1,-100000),
    (-1,0,3,0,-8000),(6,0,-18,1,-24000),
)

A=(Q(3,4),Q(1,2),Q(1,4))
B0=(Q(5,6),Q(1,2),Q(1,6))
G=(Q(6,7),Q(1,2),Q(2,7))
D=(Q(9,10),Q(1,2),Q(1,5))
E0=(Q(1),Q(2,3),Q(1,3))
F=(Q(5,6),Q(1,2),Q(1,3))
POSITIVE_TETRAHEDRA=((E0,A,B0,D),(E0,A,D,G))
NEGATIVE_TETRAHEDRA=((E0,A,G,F),)


def evaluate(row,base,q):
    return sum(a*b for a,b in zip(row,(*q,base,1)))


def source_side(base,q):
    if not valid(base,q) or len(q)!=3:
        return None
    for side,domain in (('positive',SOURCE_POSITIVE),('negative',SOURCE_NEGATIVE)):
        if all(evaluate(row,base,q)>=0 for row in domain):
            return side
    return None


def bulk_exit(base,q):
    """Exact maximal repetitions of the strict interleaving three-step word.

    The returned endpoint is the first anchor at which that word's complete
    tie-broken guard conjunction fails. It is not first arrival to a cycle.
    """
    certificate=word(base,q,3)
    expected=((1,0,0,0,-3),(0,1,0,0,-2),(0,0,1,0,-1))
    order=(3,0,4,1,5,2,6)
    if certificate['matrix']!=expected or certificate['orders']!=(order,)*3:
        raise ValueError('Outside the declared interleaving bulk word.')
    repetitions,delta=progression(base,q,certificate)
    target=tuple(x+repetitions*d for x,d in zip(q,delta))
    return {'repetitions':repetitions,'elapsed':3*repetitions,'target':target,
            'guards':certificate['guards'],'offsets':certificate['offsets']}


def negative_offsets(e,m):
    if e not in range(12,16) or m not in (-5,-4,-3,-2,-1) or (e-m)%2:
        raise ValueError('Not a certified negative exit phase.')
    return e+2, (e+m)//2-1 if m>=-1 else (e-m)//2-2


def source_entry(base,q):
    side=source_side(base,q)
    if side is None:
        raise ValueError('Outside the buffered canonical source.')
    u,v,w=q
    S,m,d=u+w-base,2*v-base,u-3*w
    active=S if side=='positive' else 2*u-2*w-base
    threshold=16 if side=='positive' else 12
    blocks=(active-threshold)//4
    final_m=phase_iterate(m,blocks)
    target=(u-3*blocks,(base+final_m)//2,w-blocks)
    elapsed=3*blocks
    if side=='negative':
        e=active-4*blocks
        expected_offsets=negative_offsets(e,final_m)
        target=subtract_pair(base,target)
        elapsed+=1
        if (2*target[1]-base,target[0]+target[2]-base)!=expected_offsets:
            raise ArithmeticError('Negative exit phase mismatch.')
    if not edge_domain(base,target):
        raise ArithmeticError('Certified return did not reach the edge.')
    offsets=(2*target[1]-base,target[0]+target[2]-base)
    rate=phase_data(base%2,offsets)['time_per_unit_u']
    return {'side':side,'blocks':blocks,'elapsed':elapsed,'target':target,
            'rate':rate,'offsets':offsets,'d':d}


def limiting_times(x,z,side):
    if side=='positive':
        return (9*x-3*z-6)/8,(5*x-3*z-3)/4
    if side=='negative':
        return 3*(x+z-1)/4,(2*x+6*z-3)/4
    raise ValueError('Unknown source side.')


def source_clock(base,q):
    entry=source_entry(base,q)
    times=limiting_times(Q(q[0],base),Q(q[2],base),entry['side'])
    return {**entry,'limiting_time':times[entry['rate']==2],
            'actual_error_bound':216}


def source_arrival(base,q):
    entry=source_entry(base,q)
    remaining,period=edge_arrival(base,entry['target'])
    return {**entry,'pair_time':entry['elapsed']+remaining,'period':period}


def period_region_arrival(base,q,first_image=False):
    """Try the <=2 canonical rotations, optionally after one initial step.

    A successful result is exact. Rejection means the buffered source is not
    reached by this prefix; it does not assert nonentry or a different scale.
    """
    p=subtract_pair(base,q) if first_image else q
    elapsed=int(first_image)
    for rotation in range(3):
        if source_side(base,p):
            answer=source_arrival(base,p)
            return {**answer,'pair_time':elapsed+answer['pair_time'],
                    'prefix':elapsed,'canonical_state':p}
        p=subtract_pair(base,p)
        elapsed+=1
    raise ValueError('No buffered canonical source in the declared prefix.')


def leading_bulk(q):
    u,v,w=q
    return u,u-w,v-w


def region_tetrahedra():
    pieces=[*POSITIVE_TETRAHEDRA,*NEGATIVE_TETRAHEDRA]
    result=[]
    for _ in range(3):
        result.extend(pieces)
        pieces=[tuple(leading_bulk(q) for q in t) for t in pieces]
    return tuple(result)


def inverse_leading(q,signs):
    u,v,w=q
    a=2*u-1;b=2*v-a;c=2*w-b
    return tuple(sorted(((1+s*r)/2 for s,r in zip(signs,(a,b,c))),reverse=True))


def population_masses():
    pieces=region_tetrahedra()
    direct=sum(tetra_weight(t) for t in pieces)
    predecessor=sum(tetra_weight(tuple(inverse_leading(q,signs) for q in t))
                    for signs in product((-1,1),repeat=3) for t in pieces)
    return {'positive':sum(tetra_weight(t) for t in POSITIVE_TETRAHEDRA),
            'negative':sum(tetra_weight(t) for t in NEGATIVE_TETRAHEDRA),
            'period_region':direct,'one_step_predecessor':predecessor}
