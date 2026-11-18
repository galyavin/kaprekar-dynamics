"""Nine-slot leading return region and exact timed bulk segments.

This module does not give terminal arrival times. In particular, a timed bulk
segment is not asserted transient until a terminal/first-capture proof supplies
that fact. Leading entry and its weight are distinct from a global entry law.
"""
from fractions import Fraction as Q

from odd_gap import subtract_pair, valid

# Homogeneous rows act on (a,b,c,d,B), with nonnegative slack.
CORE_FACETS = (
    (-2,0,2,1,1), (-1,-1,-1,-1,2), (-1,-1,1,1,1),
    (0,-1,-1,0,1), (0,0,0,1,0), (0,0,1,-1,0),
    (0,1,-1,0,0), (1,-1,0,0,0), (1,0,1,1,-1),
    (1,2,1,0,-2), (2,0,-1,0,-1), (2,1,-2,-2,-1),
)
RETURN_MATRIX = (
    (1,0,0,0,0), (0,0,-1,-1,1), (1,1,0,0,-1),
    (0,-1,-1,0,1), (0,0,0,0,1),
)
BRANCH_GUARDS = (
    (1,-1,0,0,0),(0,1,-1,0,0),(0,0,1,-1,0),(0,0,0,1,0),
    (-1,0,0,0,1),(1,1,0,0,-1),(0,-1,-1,0,1),
    (1,0,1,1,-1),(-1,-1,-1,-1,2),(1,2,1,0,-2),
)
DELTA = (-3,-2,-1,-1)
ORDER = (4,0,5,1,6,2,7,3,8)
# The six decreasing core facets, in difference coordinates (u,v,w,z,B).
EXIT_FACETS = (
    ('A0',(1,0,0,1,-1),4),
    ('A1',(1,1,-1,0,-1),4),
    ('A2',(2,-1,0,0,-1),4),
    ('B0',(0,1,1,0,-1),3),
    ('B1',(1,1,0,-2,-1),3),
    ('B2',(2,-2,2,-1,-1),3),
)


def evaluate(row,values):
    return sum(a*b for a,b in zip(row,values))


def compose(row,matrix):
    return tuple(sum(row[i]*matrix[i][j] for i in range(len(row)))
                 for j in range(len(matrix[0])))


def H(r,B=1):
    return tuple(sorted((r[0],*[abs(a+b-B) for a,b in zip(r,r[1:])]),reverse=True))


def C(r,B=1):
    return (Q(B+r[0],2),*[Q(a+b,2) for a,b in zip(r,r[1:])])


def inverse_C(q,B=1):
    r=[];previous=B
    for x in q:
        previous=2*x-previous;r.append(previous)
    return tuple(r)


def in_core_radius(r,B=1,strict=False):
    slacks=[evaluate(g,(*r,B)) for g in CORE_FACETS]
    return min(slacks)>0 if strict else min(slacks)>=0


def in_core(q,B=1,strict=False):
    return in_core_radius(inverse_C(q,B),B,strict)


def leading_bulk(q,B=1):
    u,v,w,z=q
    return u,u-z,v-z,v-w


def inverse_leading(r,signs,B=1):
    if len(signs)!=4 or any(s not in (-1,1) for s in signs):
        raise ValueError('Four signs are required.')
    return tuple(sorted((Q(B+s*x,2) for s,x in zip(signs,r)),reverse=True))


def bulk_word():
    """Reconstruct three strict interleaving steps and every affine guard."""
    matrix=tuple(tuple(int(i==j) for j in range(6)) for i in range(4))
    base=(0,0,0,0,1,0);one=(0,0,0,0,0,1)
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    guards=[];offsets=[]
    for _ in range(3):
        offsets.append(matrix)
        guards.extend((sub(matrix[-1],one),sub(sub(base,one),matrix[0]),
                       *[sub(a,b) for a,b in zip(matrix,matrix[1:])]))
        ds=(*matrix[:-1],sub(matrix[-1],one),sub(base,one),
            *[sub(sub(base,r),one) for r in reversed(matrix[1:])],sub(base,matrix[0]))
        guards.extend(sub(sub(ds[i],ds[j]),one) for i,j in zip(ORDER,ORDER[1:]))
        matrix=tuple(sub(ds[ORDER[i]],ds[ORDER[-1-i]]) for i in range(4))
    return matrix,tuple(guards),tuple(offsets)


BULK_MATRIX,BULK_GUARDS,BULK_OFFSETS=bulk_word()


def bulk_exit(B,q):
    """Maximal number of consecutive three-step words with strict guards.

    'elapsed' is time to exit this word contract, not terminal arrival time.
    A tie-resolved actual word can persist a bounded number of further steps.
    """
    if not valid(B,q) or len(q)!=4:
        raise ValueError('Require valid nine-slot difference coordinates.')
    limits=[]
    for guard in BULK_GUARDS:
        value=evaluate(guard,(*q,B,1))
        if value<0:raise ValueError('Initial strict three-step guards fail.')
        slope=evaluate(guard[:4],DELTA)
        if slope<0:limits.append(value//(-slope)+1)
    if not limits:raise ArithmeticError('Missing finite bulk exit guard.')
    count=min(limits)
    target=tuple(x+count*d for x,d in zip(q,DELTA))
    return {'repetitions':count,'elapsed':3*count,'target':target}


def limiting_exit(x):
    """Leading elapsed bulk time and boundary location, on C(K) only."""
    x=tuple(map(Q,x))
    if not in_core(x,strict=True):raise ValueError('Require the strict leading core.')
    candidates=[(Q(evaluate(row,(*x,1)),speed),key)
                for key,row,speed in EXIT_FACETS]
    t=min(v for v,key in candidates)
    return {'elapsed_over_B':3*t,'target':tuple(v+t*d for v,d in zip(x,DELTA)),
            'faces':tuple(key for v,key in candidates if v==t)}


def first_image_segment(B,q):
    """One actual subtraction followed by a certified bulk segment.

    Success supplies an exact timed path. It does not certify that no state
    along the path is already cyclic. Small bases or boundary inputs may fail
    the strict guards even when their limiting image lies in the core.
    """
    if not valid(B,q) or len(q)!=4:raise ValueError('Invalid nine-slot input.')
    r=tuple(sorted((abs(2*x-B) for x in q),reverse=True))
    if not in_core_radius(r,B,strict=True):
        raise ValueError('Leading first image is outside the strict core.')
    y=subtract_pair(B,q);segment=bulk_exit(B,y)
    return {**segment,'entry':y,'elapsed':1+segment['elapsed'],
            'leading_entry':C(tuple(Q(x,B) for x in r))}
