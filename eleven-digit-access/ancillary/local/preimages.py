"""Complete difference preimages via the borrowed digits' mirror sums."""

from itertools import permutations


MACROS = ((4,3,3,2,1), (4,3,2,1,1), (4,3,3,2,1), (4,3,2,1,1))
RESIDUALS = ((1,0,-1,0,2), (0,0,-2,0,-1), (1,0,0,1,2), (0,0,-2,1,-1))


def cycle(m):
    return {tuple(m * a + b for a, b in zip(p, e))
            for p, e in zip(MACROS, RESIDUALS)}


def from_digits(digits):
    return tuple(digits[j] - digits[10-j] for j in range(5))


def borrow(q, base):
    h = sum(x > 0 for x in q)
    if not h:
        return (0,)*11
    return (*[base-1]*(11-2*h), *q[:h-1], q[h-1]-1,
            base-q[0], *[base-x-1 for x in q[1:h]])


def step(q, base):
    return from_digits(sorted(borrow(q, base), reverse=True))


def whole_step(q, base):
    word = (*q, 0, 0, 0, 0, 0, 0)
    hi = lo = 0
    for left, right in zip(word, reversed(word)):
        hi = hi*base + left
        lo = lo*base + right
    n = hi-lo
    digits = []
    for _ in word:
        n, d = divmod(n, base)
        digits.append(d)
    assert n == 0
    return from_digits(sorted(digits, reverse=True))


def mirror_sums(h):
    zero = (0,)*h
    yield zero
    if h == 1:
        # The sole pair is (q1-1, B-q1).
        return
    for high, low in permutations(range(h), 2):
        row = list(zero)
        row[high], row[low] = 1, -1
        yield tuple(row)


def infer_digits(q, base, h, deviations):
    """Recover sorted digits from their differences and core mirror sums."""
    copies = 11-2*h
    span = base-1
    digits = [span]*copies + [None]*(2*h)

    def set_digit(i, x):
        if digits[i] is not None:
            return digits[i] == x
        digits[i] = x
        return True

    for i in range(5):
        bottom = 10-i
        assert digits[i] is not None
        if not set_digit(bottom, digits[i]-q[i]):
            return None
        if i < h and not set_digit(copies+i, span+deviations[i]-digits[bottom]):
            return None
    assert all(x is not None for x in digits)
    if digits[-1] < 0 or any(x < y for x, y in zip(digits, digits[1:])):
        return None
    if any(digits[copies+i]+digits[10-i] != span+d
           for i, d in enumerate(deviations)):
        return None
    return digits[copies:]


def remove_pair(digits, x, y):
    remaining = list(digits)
    try:
        remaining.remove(x)
        remaining.remove(y)
    except ValueError:
        return
    return remaining


def inner_choices(digits, base, count):
    """Use each remaining sum-(B-1) pair and choose its direct member."""
    if count == 0:
        yield ()
        return
    x = digits[0]
    y = base-1-x
    remaining = remove_pair(digits, x, y)
    if remaining is not None:
        for chosen in inner_choices(remaining, base, count-1):
            for p in {x, y}:
                yield (*chosen, p)


def decode(digits, base, h):
    if h == 1:
        for x in set(digits):
            yield (x+1, 0, 0, 0, 0)
        return
    for q1 in set(digits):
        fewer = remove_pair(digits, q1, base-q1)
        if fewer is None:
            continue
        for last_digit in set(fewer):
            rest = remove_pair(fewer, last_digit, base-2-last_digit)
            if rest is None:
                continue
            qh = last_digit+1
            for selected in inner_choices(rest, base, h-2):
                yield (q1, *sorted(selected, reverse=True), qh, *[0]*(5-h))


def preimages(q, base):
    found = set()
    for h in range(1,6):
        for deviations in mirror_sums(h):
            digits = infer_digits(q, base, h, deviations)
            if digits is None:
                continue
            for source in decode(digits, base, h):
                if (source[0] < base and source[h-1] > 0
                        and all(a >= b for a, b in zip(source, source[1:]))
                        and step(source, base) == q):
                    found.add(source)
    return found


def basin(m):
    base = 5*m
    seen = layer = cycle(m)
    levels = [len(layer)]
    while layer:
        layer = {p for q in layer for p in preimages(q, base)} - seen
        levels.append(len(layer))
        seen = seen | layer
    return seen, levels
