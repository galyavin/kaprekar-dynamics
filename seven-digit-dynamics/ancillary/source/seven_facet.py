"""Positive-source edge transport. The rational witnesses are in
seven_facet_certificates.json, with exact replay in seven_facet_checks.py.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial

from odd_gap import valid
from seven_boundary import edge_domain, edge_arrival, phase_data

# Rows are in (u,v,w,B,1), interpreted as >=0. ONE permits positive slack in
# Farkas certificates and is itself an unconditional nonnegative row.
RETURN_DOMAIN = (
    (0,0,0,0,1),(0,0,0,1,-100000),(-1,0,0,1,-1000),
    (0,0,1,0,-1000),(0,2,0,-1,6),(1,0,1,-1,-20),
    (1,-2,1,0,-20),(1,0,-3,0,-8000),(-3,0,9,1,-24000),
    (0,1,-2,0,-20),(1,-1,-1,0,-20),(-1,-1,2,1,-20),
    (-2,2,-1,1,-20),
)

SOURCE_DOMAIN = (
    (0,0,0,1,-100000),(-1,0,0,1,-1000),(0,0,1,0,-1000),
    (0,2,0,-1,0),(1,-2,1,0,-1000),(1,0,-3,0,-8000),
    (-3,0,9,1,-24000),(0,1,-2,0,-1000),(1,-1,-1,0,-1000),
    (-1,-1,2,1,-1000),(-2,2,-1,1,-1000),
)

VERTICES = (
    (Q(3,4),Q(1,2),Q(1,4)),
    (Q(5,6),Q(1,2),Q(1,6)),
    (Q(7,8),Q(1,2),Q(1,4)),
    (Q(9,10),Q(1,2),Q(1,5)),
    (Q(1),Q(2,3),Q(1,3)),
)
# The limiting source is a pyramid, split along the A--D base diagonal.
TETRAHEDRA = ((4,0,1,3),(4,0,3,2))


def ev(row,base,q):
    return sum(a*b for a,b in zip(row,(*q,base,1)))


def source_domain(base,q):
    return valid(base,q) and len(q)==3 and all(ev(row,base,q)>=0 for row in SOURCE_DOMAIN)


def phase_iterate(m,blocks):
    """Exact iterate of m -> |m+1|-5 in O(1) integer operations."""
    if not isinstance(blocks,int) or blocks<0 or not isinstance(m,int) or m < -6:
        raise ValueError('Require integer m>=-6 and blocks>=0.')
    if m>=0:
        count=min(blocks,m//4+1)
        m-=4*count
        blocks-=count
    seen,path={},[]
    while m not in seen:
        if blocks==0:
            return m
        seen[m]=len(path);path.append(m)
        m=abs(m+1)-5
        blocks-=1
    cycle=path[seen[m]:]
    return cycle[blocks%len(cycle)]


def edge_entry(base,q):
    if not source_domain(base,q):
        raise ValueError('Outside the positive-weight facet source domain.')
    u,v,w=q
    S,m=u+w-base,2*v-base
    blocks=(S-16)//4
    final_m=phase_iterate(m,blocks)
    target=(u-3*blocks,(base+final_m)//2,w-blocks)
    if not edge_domain(base,target):
        raise ArithmeticError('Proved facet transport did not reach the edge.')
    rate=phase_data(base%2,(final_m,S-4*blocks))['time_per_unit_u']
    return {'blocks':blocks,'elapsed':3*blocks,'target':target,
            'rate':rate,'offsets':(final_m,S-4*blocks)}


def limiting_times(x,z):
    """Fast and slow clocks on the limiting source; y controls only phase."""
    return (9*x-3*z-6)/8,(5*x-3*z-3)/4


def clock(base,q):
    entry=edge_entry(base,q)
    u,_,w=q
    fast,slow=limiting_times(Q(u,base),Q(w,base))
    tau=fast if entry['rate']==Q(3,2) else slow
    return {**entry,'limiting_time':tau,'actual_error_bound':216}


def exact_arrival(base,q):
    """Exact pair and actual times for any input in this source domain."""
    entry=edge_entry(base,q)
    h,period=edge_arrival(base,entry['target'])
    pair_time=entry['elapsed']+h
    return {'pair_time':pair_time,'actual_time':pair_time+1,'period':period,**entry}


@lru_cache(None)
def residue_rate(parity,S_mod8,m_mod8):
    if parity not in (0,1) or m_mod8%2!=parity or not 0<=S_mod8<8 or not 0<=m_mod8<8:
        raise ValueError('Invalid residue data.')
    S,m=8000+S_mod8,800+m_mod8
    blocks=(S-16)//4
    return phase_data(parity,(phase_iterate(m,blocks),S-4*blocks))['time_per_unit_u']


def conditional_rate_weights(parity):
    if parity not in (0,1):
        raise ValueError('Require base parity 0 or 1.')
    return {Q(3,2):Q(1,2) if parity==0 else Q(3,8),
            Q(2):Q(1,2) if parity==0 else Q(5,8)}


def determinant(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def tetra_weight(vertices):
    """Exact integral of 5040(1-x)(x-y)(y-z)z over a tetrahedron."""
    jac=abs(determinant(*(tuple(a-b for a,b in zip(v,vertices[0])) for v in vertices[1:])))
    forms=[(1-x,x-y,y-z,z) for x,y,z in vertices]
    total=Q(0)
    # Barycentric monomial integral: product(k_i!)/(3+4)!; 5040 cancels.
    for indices in product(range(4),repeat=4):
        term=Q(1)
        for j,i in enumerate(indices):
            term*=forms[i][j]
        for i in range(4):
            term*=factorial(indices.count(i))
        total+=term
    return jac*total


def source_mass():
    return sum(tetra_weight(tuple(VERTICES[i] for i in ids)) for ids in TETRAHEDRA)
