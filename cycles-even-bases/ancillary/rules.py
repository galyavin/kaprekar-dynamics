"""Compare the boundary tables and examples with positional four-digit subtraction."""
import json
import sys
from collections import Counter


def subtract(base, digits):
    descending = sorted(digits, reverse=True)
    # Work from the units place so leading zeroes and every borrow are retained.
    borrow = 0
    result = []
    for small, large in zip(reversed(descending), descending):
        digit = small - large - borrow
        borrow = int(digit < 0)
        result.append(digit + borrow * base)
    assert borrow == 0
    return tuple(reversed(result))


def pair(digits):
    a, b, c, d = sorted(digits, reverse=True)
    return a - d, b - c


def step(base, state):
    u, v = state
    return pair(subtract(base, (u, v, 0, 0)))


def boundary(base, state):
    u, v = state
    return v == 0 or u == v or u + v == base


def regular(m, state):
    u, v = state
    return all(x % m for x in (u, v, u - v, u + v))


def check_examples():
    examples = [
        (10, [(6, 1, 7, 4)]),
        (20, [(12, 3, 15, 8)]),
        (20, [(7, 19, 19, 12), (12, 6, 12, 8), (6, 3, 15, 14),
              (12, 7, 11, 8), (5, 2, 16, 15), (14, 9, 9, 6)]),
        (56, [(24, 7, 47, 32), (40, 7, 47, 16), (40, 23, 31, 16)]),
        (56, [(7, 55, 55, 48), (48, 6, 48, 8),
              (42, 39, 15, 14), (28, 23, 31, 28)])
    ]
    for base, cycle in examples:
        assert len(set(cycle)) == len(cycle)
        for i, digits in enumerate(cycle):
            assert subtract(base, digits) == cycle[(i + 1) % len(cycle)], (base, digits)
    assert pair(examples[0][1][0]) == (6, 2)
    return len(examples)


def boundaries(base):
    for u in range(1, base):
        yield (u, 0)
        yield (u, u)
        if u > base // 2:
            yield (u, base - u)


def table(m, n):
    s, h = 2**n, (m + 1) // 2
    base = m * s
    entries = {}

    def enter(source, target, elapsed):
        assert source not in entries, (m, n, source)
        entries[source] = (target, elapsed)

    enter((base // 2, base // 2), (1, 1), 1)
    enter((1, 0), (base - 1, 0), 1)
    for t in range(1, n + 1):
        length = base // 2**t
        for i in range(1, 2**t, 2):
            for u in (length * i, length * i + 1):
                enter((u, 0), (2**t, 0), t + 1)
    for u, v in boundaries(base):
        if v == 0:
            if u % m == h:
                enter((u, v), (h * s, base - h * s), n + 1)
        elif u % m in (h - 1, h):
            enter((u, v), (2 * s, 0), n + 1)
        elif u % m == 0:
            order = (u & -u).bit_length() - 1
            t = n - order - 2
            if t >= 0:
                a = base // 2 + 2**t
                enter((u, v), (a, base - a), t + 1)
    return entries


def first_return(base, start, cache):
    q = step(base, start)
    path, seen = [], set()
    while not boundary(base, q) and q not in cache and q not in seen:
        seen.add(q)
        path.append(q)
        q = step(base, q)
    # A repeated generic state proves there will be no boundary visit.
    answer = (q, 0) if boundary(base, q) else cache.get(q)
    for p in reversed(path):
        if answer is not None:
            answer = (answer[0], answer[1] + 1)
        cache[p] = answer
    return (answer[0], answer[1] + 1) if answer is not None else None


def check_base(m, n):
    s, h = 2**n, (m + 1) // 2
    base, sigma = m * s, s % m
    S, E, U = (s, 0), (2*s, 0), (h*s, base - h*s)
    named = {S: 'S', E: 'E', U: 'U'}
    critical = {(1, 1), (base - 1, 0), E, U}
    critical.update((2**t, 0) for t in range(1, n + 1))
    critical.update((base//2 + 2**t, base//2 - 2**t) for t in range(n - 1))
    declared = table(m, n)
    actual, cache = {}, {}
    for q in boundaries(base):
        ret = first_return(base, q, cache)
        assert ret == declared.get(q), (m, n, q, ret, declared.get(q))
        actual[q] = ret
    assert {r[0] for r in actual.values() if r} == critical, (m, n)
    arrows = {
        S: S if sigma == 1 else U if sigma == h else None,
        E: S if sigma == h else U if (4*sigma) % m == 1 else None,
        U: E if sigma in (1, m - 1) else None
    }
    critical_seeded = 0
    for q in sorted(critical):
        u, v = q
        residue = u % m
        assert residue != 0
        if v == 0:
            target = S if residue == 1 else U if residue == h else None
        else:
            target = E if residue in (h - 1, h) else None
        assert regular(m, step(base, q)) == (target is None), (m, n, q)
        critical_seeded += target is None
        assert actual[q] == ((target, n + 1) if target else None), (m, n, q)
        if q in named:
            assert arrows[q] == target, (m, n, q)
    return dict(m=m, n=n, base=base, boundaries=len(actual),
                returns=sum(ret is not None for ret in actual.values()),
                critical=len(critical), critical_seeded=critical_seeded,
                arrows={named[q]: named.get(target) for q, target in arrows.items()})


def main():
    if not __debug__:
        raise SystemExit("Run Python with assertions enabled.")
    examples = check_examples()
    cases = [check_base(m, n) for m in range(3, 42, 2) for n in range(1, 7)]
    counts = Counter()
    for row in cases:
        for field in ('boundaries', 'returns', 'critical', 'critical_seeded'):
            counts[field] += row[field]
    original = sum(r['boundaries'] for r in cases if r['m'] <= 31 and r['n'] <= 4)
    json.dump(dict(digit_examples=examples, paper_range_boundaries=original,
                   totals=dict(counts), cases=cases), sys.stdout, indent=2)
    print()


if __name__ == '__main__':
    main()
