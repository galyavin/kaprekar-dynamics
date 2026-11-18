"""Raw-sort reference and exceptional cycles."""
from fractions import Fraction as Q
from nine_entry import CORE_FACETS,evaluate


def reference(r):
    # Conjugate a leading nine-digit subtraction back to radius coordinates;
    # this does not use the four-value radius map or its sorting network.
    q=(Q(1+r[0],2),*(Q(a+b,2) for a,b in zip(r,r[1:])))
    digits=sorted((*q,Q(1),*(1-v for v in reversed(q))),reverse=True)
    pair=tuple(digits[i]-digits[-1-i] for i in range(4))
    out=[];previous=Q(1)
    for value in pair:previous=2*value-previous;out.append(previous)
    return tuple(out)


def outside(q):return min(evaluate(g,(*q,1)) for g in CORE_FACETS)<0


def known(q,n):
    a,b,c,d=q
    fixed=all(v==Q(1,3) for v in q)
    if n==2:
        return fixed or (a==Q(3,5) and d==Q(1,5) and ((b,c)==(a,a) or (b,c)==(d,d)))
    return fixed or (Q(5,7)<=a<=Q(3,4) and (
        (b==c==2*a-1 and 0<d<=1-a) or
        (c==3*a-2 and d==3-4*a and 1-a<=b<2-2*a) or
        (b==a and c+d==2-2*a and 0<=d<1-a)))
