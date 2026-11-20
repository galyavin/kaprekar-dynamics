"""Exact, finite expression for the independent phase work at positive anchor."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from faces import layout

def coeffs(g):
    return tuple(2*(b-a) for a,b in zip(g,g[1:]))

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

# For each digit group, record a vector of source q offsets and its ordered
# digit offsets. Shifts are those of the positive-innermost borrowed list.
@lru_cache(maxsize=None)
def options(k, group, center):
    ids=tuple(i for i,s in group)
    signs=tuple(s for i,s in group)
    n=len(ids)
    if center:
        grids=(tuple(2*x for x in range(-n,n+1)),
               tuple(2*x+1 for x in range(-n,n)))
    else:
        grids=(tuple(2*x for x in range(n)),)
    result=[]
    for grid in grids:
        for xs in product(grid,repeat=n):
            if not center and 0 not in xs:
                continue
            ds=xs if center else tuple(s*x for s,x in zip(signs,xs))
            if any(ds[j]<ds[j+1] for j in range(n-1) if signs[j]==signs[j+1]):
                continue
            if center:
                shifts=tuple(d-2*(i==k-1) for i,d in zip(ids,ds)) + tuple(-d-2*(i!=0) for i,d in zip(ids,ds))
                hi=tuple(sorted(shifts,reverse=True))
                lo=()
            else:
                hi=tuple(sorted((x-2*(i==k-1 if s==1 else i!=0) for i,s,x in zip(ids,signs,xs)),reverse=True))
                lo=tuple(sorted((-x-2*(i!=0 if s==1 else i==k-1) for i,s,x in zip(ids,signs,xs)),reverse=True))
            result.append((hi+lo,ds))
    return tuple(sorted(set(result)))

@lru_cache(maxsize=None)
def forms(spec,k):
    middle,groups,ep=spec
    assert not ep
    positions=layout(spec)
    high_from=1
    low_until=2*k+1
    result=[]
    for t,(h,l) in enumerate(groups):
        n=h+l
        indices=tuple((i,s) for i,(r,s) in enumerate(positions) if r==t)
        a=tuple(range(high_from,high_from+n))
        b=tuple(range(low_until-n,low_until))
        result.append((tuple(i for i,s in indices), a+b, options(k,indices,False)))
        high_from+=n
        low_until-=n
    if middle:
        indices=tuple((i,s) for i,(r,s) in enumerate(positions) if r is None)
        result.append((tuple(i for i,s in indices),tuple(range(high_from,low_until)),options(k,indices,True)))
    return tuple(result)

def w_phase(g,h,spec):
    k=len(g)-1
    source=coeffs(g)
    target=coeffs(h)
    digits=target+(0,)+tuple(-x for x in target[::-1])
    total=-digits[0]
    for ids, idx, opts in forms(spec,k):
        c=tuple(source[i] for i in ids)
        d=tuple(digits[j] for j in idx)
        total+=Q(1,2)*max(dot(d,u)-dot(c,v) for u,v in opts)
    return total
