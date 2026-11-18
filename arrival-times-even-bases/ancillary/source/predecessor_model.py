from itertools import product


def ordered(x, y):
    return (max(x, y), min(x, y))


def valid(base, p):
    return 1 <= p[0] < base and 0 <= p[1] <= p[0]


def predecessors(base, q):
    """Complete inverse image, including a fixed point's self predecessor.

    Separate input-branch guards avoid relying on a forward transition test.
    Sets remove sign symmetries and overlap between diagonal/anti-diagonal
    descriptions of the same input.
    """
    if not valid(base, q):
        raise ValueError('Invalid difference pair')
    a, b = q
    result = set()
    for s, t in product((-1, 1), repeat=2):
        x, y = base+s*a, base+t*b
        if x % 2 or y % 2:
            continue
        p = ordered(x//2, y//2)
        u, v = p
        if valid(base, p) and u > v > 0 and u+v != base:
            result.add(p)
    if a+b == base-1:
        for u in {a+1, b+1}:
            if 1 <= u < base:
                result.add((u, 0))
    hs = {0} if q == (1, 1) else ({a-1, 1-a} if a-b == 2 else set())
    for h in hs:
        if (base+h) % 2:
            continue
        u = (base+h)//2
        if 1 <= u < base:
            result.add((u, u))
            result.add(ordered(u, base-u))
    return result
