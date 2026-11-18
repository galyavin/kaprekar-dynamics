"""Leading-entry representation of the global seven-digit limit kernel.

Exact rational evaluation is supported away from the theorem's null exceptional
set. A finite iteration cap is an evaluation limit, not an entry disproof.
"""
from fractions import Fraction as Q

from odd_gap import leading_pair, subtract_pair
from seven_return_law import (leading_bulk,limiting_times,source_side,source_arrival)


def H(q):
    a,b,c=q
    return tuple(sorted((a,abs(a+b-1),abs(b+c-1)),reverse=True))


def defect(q):
    a,b,c=q
    return sum(max(0,t) for t in (1-a-b,b+c-a,2*a-1-b-2*c,2*b+c+a-2))


def C(q):
    a,b,c=q
    return (1+a)/2,(a+b)/2,(b+c)/2


def inverse_C(q):
    u,v,w=q
    a=2*u-1;b=2*v-a;c=2*w-b
    return a,b,c


def in_K(q,strict=False):
    a,b,c=q
    slacks=(1-a,a-b,b-c,c,2*a-1,a+b-1,a-b-c,
            b+2*c-2*a+1,2-a-2*b-c)
    return min(slacks)>0 if strict else min(slacks)>=0


def canonical(q):
    """Unique core rotation off the primary-distance and d=0 seams."""
    y=q
    for j in range(3):
        u,v,w=y
        if u+w>2*v and u>v+w:
            d=u-3*w
            if not d:raise ValueError('Excluded zero-separation seam.')
            return y,j,'positive' if d>0 else 'negative'
        y=leading_bulk(y)
    raise ValueError('Excluded core rotation seam.')


def limit_kernel(x,parity,max_steps=100000):
    """Return the two limiting clocks and weights at a normalized initial state."""
    x=tuple(Q(v) for v in x)
    if parity not in (0,1) or not (1>x[0]>=x[1]>=x[2]>0):
        raise ValueError('Require an interior ordered state and a base parity.')
    r=tuple(sorted((abs(2*v-1) for v in x),reverse=True))
    seen=set();steps=1
    while not in_K(r):
        if r in seen:raise ValueError('Exceptional periodic leading orbit.')
        if steps>=max_steps:raise RuntimeError('Leading iteration cap; no nonentry conclusion.')
        seen.add(r);r=H(r);steps+=1
    if not in_K(r,strict=True):
        raise ValueError('Excluded core boundary.')
    y,rotation,side=canonical(C(r))
    fast,slow=limiting_times(y[0],y[2],side)
    p=Q(3,8) if parity else Q(1,2) if side=='positive' else Q(3,4)
    return {'fast':fast,'slow':slow,'fast_weight':p,'side':side,
            'canonical_state':y,'leading_prefix':steps+rotation}


def capture_arrival(base,q,max_steps=100000):
    """An exact actual-map evaluator when the predicted prefix reaches a source.

    This is not a complete finite-base solver. Null/seam inputs and insufficient
    buffers are rejected; the limiting theorem handles their population by
    compact exhaustion, not by asserting this evaluator succeeds everywhere.
    """
    kernel=limit_kernel(tuple(Q(v,base) for v in q),base%2,max_steps)
    p=q
    for _ in range(kernel['leading_prefix']):p=subtract_pair(base,p)
    if source_side(base,p) is None:
        raise ValueError('Finite base does not reach the buffered predicted source.')
    answer=source_arrival(base,p)
    return {**kernel,'pair_time':kernel['leading_prefix']+answer['pair_time'],
            'period':answer['period'],'actual_canonical_state':p}
