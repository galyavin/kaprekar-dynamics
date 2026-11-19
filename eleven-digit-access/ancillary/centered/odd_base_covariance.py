"""Exact base translation on the paired-level region for odd digit counts.

This region is not forward invariant. Iterated conclusions require the
whole orbit to remain in it, as supplied by a separate capture theorem.
"""
from itertools import product


def anchor(k):
    if k < 1:
        raise ValueError('Positive number of difference coordinates required')
    return tuple(range(k,0,-1))


def blocks(B,q):
    k=len(q)
    return tuple((q[i]-(i==k-1),B-q[k-1-i]-(i<k-1)) for i in range(k))


def in_region(B,q):
    if not q or not B>q[0] or q[-1]<1 or any(a<b for a,b in zip(q,q[1:])):
        return False
    pairs=blocks(B,q)
    return all(min(a)>=max(b) for a,b in zip(pairs,pairs[1:]))


def paired_step(B,q):
    if not in_region(B,q):
        raise ValueError('Outside the paired-level region')
    digits=[B-1]
    for pair in blocks(B,q):
        digits.extend(sorted(pair,reverse=True))
    return tuple(digits[i]-digits[-1-i] for i in range(len(q)))


def centered(B,q):
    k=len(q)
    return tuple((k+1)*v-a*B for v,a in zip(q,anchor(k)))


def centered_step(x):
    """Autonomous integer map; finite-base validity remains a guard question."""
    k=len(x);m=k+1
    offsets=[-m]
    for i in range(k):
        pair=(x[i]-m*(i==k-1),-x[k-1-i]-m*(i<k-1))
        offsets.extend(sorted(pair,reverse=True))
    return tuple(offsets[i]-offsets[-1-i] for i in range(k))


def paired_orders(k):
    for switches in product((0,1),repeat=k):
        order=[k]
        for i,swap in enumerate(switches):
            pair=(i,k+1+i)
            order.extend(reversed(pair) if swap else pair)
        yield tuple(order)


def order_is_paired(order):
    if len(order)%2!=1 or sorted(order)!=list(range(len(order))):
        return False
    k=len(order)//2
    return order[0]==k and all(set(order[2*i+1:2*i+3])=={i,k+1+i} for i in range(k))
