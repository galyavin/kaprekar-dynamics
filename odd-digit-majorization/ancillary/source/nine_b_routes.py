"""All-base weighted nine-digit routes through B, BB and BBB.

The six original-input cells retain every integer phase. All have a common
leading time formula and terminal position. This is a conditional population
law; it is not full weighted exhaustion or global leading entry.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd,lcm

from nine_entry import CORE_FACETS,EXIT_FACETS,DELTA,compose,evaluate,bulk_exit
import nine_b_facet as facet
import nine_bb as corner
import nine_bbb as terminal
from odd_gap import subtract_pair
from odd_affine import word

I=tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
CM=((Q(1,2),0,0,0,Q(1,2)),(Q(1,2),Q(1,2),0,0,0),
    (0,Q(1,2),Q(1,2),0,0),(0,0,Q(1,2),Q(1,2),0),(0,0,0,0,1))
L=((1,0,0,0,0),(1,0,0,-1,0),(0,1,0,-1,0),(0,1,-1,0,0),(0,0,0,0,1))
MASS=Q(234391778,38443359375)


def normalized(row):
    d=lcm(*(Q(v).denominator for v in row));values=[int(Q(v)*d) for v in row];g=gcd(*values)
    return tuple(v//g for v in values) if g else tuple(values)


def multiply(a,b):return tuple(compose(row,b) for row in a)


def flow(row,speed,delta):
    return tuple(tuple(I[i][j]+(Q(delta[i]*row[j],speed) if i<4 else 0)
                       for j in range(5)) for i in range(5))


def choose(domain,matrix,selected,speed,candidates):
    for row,s in candidates:
        if row==selected:continue
        comparison=tuple(speed*a-s*b for a,b in zip(row,selected))
        r=normalized(compose(comparison,matrix))
        if not any(r):raise ArithmeticError('Forced simultaneous exit needs a separate route.')
        domain.append(r)


@lru_cache(None)
def route_specs():
    records=[]
    for first,row,speed in EXIT_FACETS:
        if not first.startswith('B'):continue
        for second,srow,sspeed in (('B1',corner.B1,3),('B2',corner.B2,6)):
            domain=list(CORE_FACETS);matrix=CM
            choose(domain,matrix,row,speed,[(g,s) for name,g,s in EXIT_FACETS])
            matrix=multiply(flow(row,speed,DELTA),matrix)
            for _ in range(int(first[1])):matrix=multiply(L,matrix)
            choose(domain,matrix,srow,sspeed,[(g,-evaluate(g[:4],facet.DELTA_B))
                    for g in facet.NON_B if evaluate(g[:4],facet.DELTA_B)<0])
            matrix=multiply(flow(srow,sspeed,facet.DELTA_B),matrix)
            for _ in range(2 if second=='B2' else 0):matrix=multiply(L,matrix)
            choose(domain,matrix,corner.B2,9,[(g,-evaluate(g[:4],corner.DELTA_BB))
                    for g in corner.NON_BB if evaluate(g[:4],corner.DELTA_BB)<0])
            target=multiply(flow(corner.B2,9,corner.DELTA_BB),matrix)
            records.append({'first':first,'second':second,'domain':tuple(sorted(set(domain))),
                            'bb_entry_matrix':matrix,'bbb_entry_matrix':target})
    return tuple(records)


def limit_summary(x):
    x=tuple(map(Q,x))
    if len(x)!=4 or not 1>x[0]>x[1]>x[2]>x[3]>0:
        raise ValueError('Require strict interior normalized differences.')
    r=tuple(sorted((abs(2*v-1) for v in x),reverse=True))
    matches=[s for s in route_specs() if min(evaluate(g,(*r,1)) for g in s['domain'])>0]
    if len(matches)!=1:raise ValueError('Outside the six strict B-to-BBB route cells.')
    s=matches[0]
    return {'route':(s['first'],s['second']),'radius':r,
            'time':Q(1+r[0],2)-Q(4,5),'terminal':(Q(4,5),Q(3,5),Q(2,5),Q(1,5))}


def guarded_facet_absorb(B,q):
    """Replay the already proved B phase formulas using their actual guards.

    This avoids a conservative 14000 sufficient buffer at finite bases. No
    step is accepted solely from the limiting geometry or from a phase count.
    """
    q=tuple(q);elapsed=0
    if abs(facet.phase(B,q))>512:raise ValueError('Outside the bounded incoming B phase range.')
    while facet.phase(B,q) not in facet.PHASES:
        key='high' if facet.phase(B,q)>=-1 else 'low';record=word(B,q,3)
        if record['orders']!=facet.THREE_ORDERS[key]:raise ValueError('B absorption guards fail.')
        q=facet.three_formula(B,q);elapsed+=3
        if elapsed>513:raise ArithmeticError('B phase absorption bound failed.')
    return {'elapsed':elapsed,'target':q}


def exact_arrival(B,q):
    q=tuple(q)
    if not isinstance(B,int) or B<1000000 or any(not isinstance(v,int) for v in q):
        raise ValueError('Require integer differences and B>=1000000.')
    limiting=limit_summary(tuple(Q(v,B) for v in q));first,second=limiting['route']
    initial=q;q=subtract_pair(B,q);entry=q;elapsed=1;stages=[]
    if entry[0]!=max(B-initial[3],initial[0]-1):raise ArithmeticError('First-image exterior identity failed.')
    def append(label,start,segment):
        if segment['target'][0]!=start[0]-segment['elapsed']:
            raise ArithmeticError('Longitudinal clock segment failed.')
        stages.append({'stage':label,**segment})
        return segment['target'],segment['elapsed']
    segment=bulk_exit(B,q);q,n=append('bulk',q,segment);elapsed+=n
    for _ in range(int(first[1])):
        y=subtract_pair(B,q)
        if y[0]!=q[0]-1:raise ValueError('First rotation leaves the B clock contract.')
        q=y;elapsed+=1
    segment=guarded_facet_absorb(B,q);q,n=append('B absorption',q,segment);elapsed+=n
    segment=facet.facet_exit(B,q);q,n=append('B facet',q,segment);elapsed+=n
    for _ in range(2 if second=='B2' else 0):
        y=subtract_pair(B,q)
        if y[0]!=q[0]-1:raise ValueError('Second rotation leaves the B clock contract.')
        q=y;elapsed+=1
    segment=corner.absorb(B,q);q,n=append('BB absorption',q,segment);elapsed+=n
    segment=corner.corner_exit(B,q);q,n=append('BB corner',q,segment);elapsed+=n
    tail=terminal.exact_arrival(B,q)
    if q[0]<=4*(B//5)+6:raise ValueError('Missing terminal exclusion for earlier stages.')
    pair=elapsed+tail['pair_time']
    if elapsed-1!=entry[0]-q[0]:raise ArithmeticError('Total prefix clock failed.')
    return {'source':initial,'entry':entry,'route':(first,second),'pair_time':pair,
            'actual_time':pair+1,'period':tail['period'],
            'first_cycle_state':tail['first_cycle_state'],'stages':stages,'tail':tail,
            'limiting_time':str(limiting['time'])}
