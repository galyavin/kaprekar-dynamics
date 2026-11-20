from collections import Counter
from fractions import Fraction as Q
import itertools


def templates(k):
    """Label a radius group (h,l) by its counts above/below 1/2."""
    def groups(n):
        if n == 0:
            yield ()
        for size in range(1, n + 1):
            for h in range(size + 1):
                for tail in groups(n - size):
                    yield ((h, size - h),) + tail
    for middle in range(k + 1):
        for layers in groups(k - middle):
            for endpoint in (False, True) if layers else (False,):
                yield middle, layers, endpoint


def form(k, c=0, j=None, a=0):
    return tuple([2*c] + [2*a if i == j else 0 for i in range(k)])


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def q_forms(k, spec):
    middle, groups, endpoint = spec
    upper, lower = [], []
    for i,(h,l) in enumerate(groups):
        hi = form(k, 1) if endpoint and i==0 else form(k,Q(1,2),i,1)
        lo = form(k, 0) if endpoint and i==0 else form(k,Q(1,2),i,-1)
        upper.extend([hi]*h)
        lower.extend([lo]*l)
    return upper + [form(k,Q(1,2))]*middle + lower[::-1]


def digits_forms(k,spec):
    middle, groups, endpoint=spec
    upper=[]; lower=[]
    for i,(h,l) in enumerate(groups):
        hi = form(k, 1) if endpoint and i==0 else form(k,Q(1,2),i,1)
        lo = form(k, 0) if endpoint and i==0 else form(k,Q(1,2),i,-1)
        upper.extend([hi]*(h+l))
        lower.extend([lo]*(h+l))
    return [form(k,1)] + upper + [form(k,Q(1,2))]*(2*middle) + lower[::-1]


def order_zero_loss(spec):
    middle, groups, endpoint = spec
    if endpoint:
        groups = groups[1:]
    circle = [(1,1)] + list(groups)
    if middle:
        circle.append((middle,middle))
    circle.extend((l,h) for h,l in groups[::-1])
    if any((a==0 and b>1) or (b==0 and a>1) for a,b in circle):
        return False
    return all(not ((a==c==0) or (b==d==0)) for (a,b),(c,d) in zip(circle,circle[1:]+circle[:1]))


def gaps(k,q):
    return [sub(a,b) for a,b in zip([form(k,1)]+q,q+[form(k)])]


def permutation(k,spec):
    q=q_forms(k,spec)
    lead=digits_forms(k,spec)
    nxt=[sub(lead[i],lead[-1-i]) for i in range(k)]
    a,b=gaps(k,q),gaps(k,nxt)
    if Counter(a)!=Counter(b):
        return None
    avail=list(enumerate(a))
    p=[]
    for x in b:
        n=next(i for i,(_,y) in enumerate(avail) if x==y)
        p.append(avail.pop(n)[0])
    return tuple(p)


if __name__=='__main__':
    for k in range(2,6):
        specs=list(templates(k))
        pcounts=Counter()
        for spec in specs:
            p=permutation(k,spec)
            assert (p is not None)==order_zero_loss(spec),spec
            if p is not None:
                pcounts[p]+=1
        print(k,len(specs),sum(pcounts.values()),len(pcounts))
        print(pcounts)


def layout(spec):
    middle, groups, _ = spec
    upper, lower = [], []
    for j, (h,l) in enumerate(groups):
        upper.extend([(j,1)]*h)
        lower.extend([(j,-1)]*l)
    return upper + [(None,0)]*middle + lower[::-1]


def spec_of(q):
    radii={}
    for x in q:
        radius=abs(2*x-1)
        counts=radii.setdefault(radius,[0,0])
        counts[0 if 2*x>1 else 1]+=1
    middle=sum(radii.pop(0,[0,0]))
    return middle,tuple(tuple(radii[r]) for r in sorted(radii,reverse=True)),1 in radii
