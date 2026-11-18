"""Exact odd-slot difference states, gap transport, and input multiplicities.

No global higher-digit arrival theorem is assumed. Direct subtraction below
uses digit columns; the chamber formula is a separate, guarded identity.
"""
from itertools import product
from math import factorial


def valid(base, q):
    return (isinstance(base, int) and base >= 2 and len(q) >= 1
            and all(isinstance(x, int) for x in q)
            and base > q[0] >= 1
            and all(a >= b for a, b in zip(q, q[1:])) and q[-1] >= 0)


def subtract_pair(base, q):
    """Independent digit-column subtraction of a sorted representative."""
    if not valid(base, q):
        raise ValueError('Require valid nonconstant odd-slot differences.')
    n = 2*len(q)+1
    digits = [*q, *([0]*(n-len(q)))]
    out, borrow = [0]*n, 0
    for j in range(n-1, -1, -1):
        z = digits[j]-digits[n-1-j]-borrow
        borrow = int(z < 0)
        out[j] = z+borrow*base
    assert borrow == 0
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(len(q)))


def leading_pair(q, denominator=1):
    digits = sorted([*q, *(denominator-x for x in q), denominator], reverse=True)
    return tuple(digits[i]-digits[-1-i] for i in range(len(q)))


def gaps(base, q):
    return (base-q[0], *(q[i]-q[i+1] for i in range(len(q)-1)), q[-1])


def permutation(n):
    if n < 3 or n % 2 == 0:
        raise ValueError('Require odd digit count at least three.')
    return tuple(j//2 if j % 2 == 0 else (n-j)//2 for j in range(1, (n+1)//2))


def permutation_cycles(n):
    p = permutation(n)
    todo, cycles = set(range(1, len(p)+1)), []
    while todo:
        start = min(todo); cycle = []; j = start
        while j not in cycle:
            cycle.append(j); todo.remove(j); j = p[j-1]
        cycles.append(tuple(cycle))
    return tuple(cycles)


def chamber_digits(base, q):
    """Claimed descending interleaving; the caller must test its order."""
    ds = [base-1]
    for i in range(len(q)):
        ds.append(q[i] if i < len(q)-1 else q[i]-1)
        complement = base-q[-1-i]
        ds.append(complement-1 if i < len(q)-1 else complement)
    return tuple(ds)


def in_chamber(base, q):
    if not valid(base, q) or q[-1] == 0:
        return False
    ds = chamber_digits(base, q)
    return all(a > b for a, b in zip(ds, ds[1:]))


def chamber_gap_step(base, q):
    if not in_chamber(base, q):
        raise ValueError('Strict chamber guards fail.')
    g = gaps(base, q)
    c = (-1,) if len(q) == 1 else (-2, 1, *([0]*(len(q)-2)))
    return (g[0]+1, *(g[p]+d for p, d in zip(permutation(2*len(q)+1), c)))


def ordered_multiplicity(base, q):
    """Exact ordered-input count in O(3^k), for 2k+1 slots.

    Split each inner gap into the two adjacent sorted-digit gaps it combines.
    Only zero/nonzero patterns affect the multinomial coefficient. There are
    g-1 interior splits and two endpoint splits when g>0.
    """
    if not valid(base, q):
        raise ValueError('Require valid differences.')
    gs = gaps(base, q)
    choices = []
    for g in gs[1:]:
        choices.append([(False, False, 1)] if not g else
                       [(False, True, 1), (True, False, 1)] +
                       ([(True, True, g-1)] if g > 1 else []))
    total, n = 0, 2*len(q)+1
    for choice in product(*choices):
        cuts = [c[0] for c in choice]+[c[1] for c in reversed(choice)]
        factor, run, denominator = 1, 1, 1
        for _, _, count in choice:
            factor *= count
        for cut in cuts:
            if cut:
                denominator *= factorial(run); run = 1
            else:
                run += 1
        denominator *= factorial(run)
        total += factor*(factorial(n)//denominator)
    return gs[0]*total
