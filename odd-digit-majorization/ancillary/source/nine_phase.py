"""A certified positive-population nine-slot phase obstruction.

Scope: the parameter box below, bases divisible by five, B >= 1,000,000.
The certificate replay is in nine_phase_checks.py; this is not a global law.
"""
from fractions import Fraction as Q
from math import factorial

MIN_BASE = 1000000
BOX = ((Q(401, 500), Q(403, 500)), (Q(1, 10000), Q(4, 10000)),
       (Q(3, 10000), Q(10, 10000)), (Q(2, 10000), Q(20, 10000)))
# Rows act on (U,M,N,P,B,1).
TERMINAL = ((1,0,0,0,0,0), (2,0,0,0,-1,0),
            (3,0,0,0,-2,0), (-1,0,0,0,1,0))
INCREMENTS = ((0,10,6,3,0,0), (0,20,12,2,0,0),
              (0,30,6,1,0,0), (0,20,4,1,0,0))
BASE_ROW, ONE = (0,0,0,0,1,0), (0,0,0,0,0,1)
SEED = (803612,389,833,1000,1000000,1)
STAGES = ((3,3,(-3,-2,-1,-1)), (2,12,(-6,-12,-6,-4)),
          (1,42,(-10,-20,-30,-20)))


def evaluate(row, parameters):
    return sum(a*b for a,b in zip(row, parameters))


def stage(k, slow=False):
    return tuple(tuple(x + (y if 1 <= j <= k else 0)
                       - (4 if slow and i == 2 and j == 5 else 0)
                       for j,(x,y) in enumerate(zip(t,v)))
                 for i,(t,v) in enumerate(zip(TERMINAL,INCREMENTS)))


def source(base, parameters, slow=False):
    parameters = tuple(parameters)
    if (not isinstance(base,int) or base < MIN_BASE or base % 5
            or len(parameters) != 4 or any(not isinstance(p,int) for p in parameters)
            or any(not lo*base <= p <= hi*base for p,(lo,hi) in zip(parameters,BOX))):
        raise ValueError('Require integral parameters in the box and B >= 10^6 divisible by 5.')
    return tuple(evaluate(r,(*parameters,base,1)) for r in stage(3,slow))


def exact_arrival(base, parameters, slow=False):
    q = source(base,parameters,slow)
    U,M,N,P = parameters
    h = 3*P+12*N+42*M
    if slow:
        h += 16*(U-4*(base//5)) + 16*((U-1-4*(base//5)) % 2)
    return {'source':q, 'pair_time':h, 'actual_time':h+1,
            'period':21 if slow else 1}


def source_mass():
    """Exact original-input limiting mass of EACH of the two lattice families.

    Integration is in parameter coordinates, hence no determinant factor.
    Each family occupies one coset of an index-1440 sublattice in state space.
    """
    rows = [tuple(Q(r[i]) for i in range(5)) for r in stage(3)]
    one = (Q(0),)*4+(Q(1),)
    sub = lambda a,b: tuple(x-y for x,y in zip(a,b))
    gaps = [sub(one,rows[0]), *[sub(a,b) for a,b in zip(rows,rows[1:])],rows[-1]]
    poly = {(0,0,0,0): Q(factorial(9))}
    for gap in gaps:
        out = {}
        for powers, coefficient in poly.items():
            for i,g in enumerate(gap):
                if not g:
                    continue
                p = tuple(v+int(i == j) for j,v in enumerate(powers))
                out[p] = out.get(p,Q(0)) + coefficient*g
        poly = out
    answer = Q(0)
    for powers, coefficient in poly.items():
        term = coefficient
        for p,(lo,hi) in zip(powers,BOX):
            term *= (hi**(p+1)-lo**(p+1))/(p+1)
        answer += term
    return answer
