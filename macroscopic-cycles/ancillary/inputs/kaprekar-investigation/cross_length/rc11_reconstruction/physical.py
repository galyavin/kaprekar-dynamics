"""Exact whole-word controls for the RC-11 reconstructions."""

from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from math import factorial, isqrt, prod


def sign(z):
    return (z > 0) - (z < 0)


def whole(word, base):
    """Subtract the complete sorted ascending word from the descending word."""
    digits = sorted(word)
    borrow = 0
    reversed_output = []
    for high, low in zip(digits, reversed(digits)):
        value = high - low - borrow
        borrow = int(value < 0)
        reversed_output.append(value + base * borrow)
    assert borrow == 0
    return tuple(reversed(reversed_output))


def differences(word):
    digits = sorted(word, reverse=True)
    return tuple(digits[i] - digits[-1 - i] for i in range(len(word) // 2))


def canonical(q):
    return q + (0,) * (len(q) + 1)


def physical(q, base):
    return differences(whole(canonical(q), base))


def columns(q, base):
    """The integer-linear trace used in the Presburger support description."""
    carry = 0
    out = []
    for a in tuple(-x for x in q) + (0,) + q[::-1]:
        z = a - carry
        carry = int(z < 0)
        out.append(z + carry * base)
    assert not carry
    return tuple(reversed(out))


def gaps(q, base):
    return (base - q[0],) + tuple(a - b for a, b in zip(q, q[1:] + (0,)))


def suffixes(g):
    total = sum(g)
    q = []
    for x in g[:-1]:
        total -= x
        q.append(total)
    return tuple(q)


def pairing(q, base):
    if q[-1] <= 0:
        return False
    levels = tuple((q[i] - int(i + 1 == len(q)), base - q[-1 - i] - int(i + 1 < len(q)))
                   for i in range(len(q)))
    return max(levels[0]) <= base - 1 and all(min(a) >= max(b) for a, b in zip(levels, levels[1:]))


def reflection(g, i):
    return sum(g[-i:]) - sum(g[:i]) if i else 0


def folded(j, n):
    return min(2 * j, n - 2 * j)


def permute(g):
    out = [0] * len(g)
    for i, x in enumerate(g):
        out[folded(i, 2 * len(g) - 1)] = x
    return tuple(out)


def update_terms(g):
    d = tuple(reflection(g, i) for i in range(len(g) // 2 + 1))
    h = tuple(max(0, -1 - x) for x in d)
    a = [0] * len(g)
    for i in range(1, len(d)):
        a[2 * i - 2] = h[i - 1] - h[i]
        a[2 * i - 1] = -a[2 * i - 2]
    e = sign(d[1])
    kick = (e, -2 * e, e) + (0,) * (len(g) - 3)
    return d, h, tuple(a), kick


def paired_update(g):
    _, _, a, kick = update_terms(g)
    return tuple(x + y + z for x, y, z in zip(permute(g), a, kick))


def sq(values):
    return sum(a * a for a in values)


def leading_output(q, base):
    radii = sorted((max(x, base - x) for x in q), reverse=True)
    return (radii[0],) + tuple(a + b - base for a, b in zip(radii, radii[1:]))


def exact_work(g, out):
    a, y1, _, y3, y4, _, y6 = g
    b, z1, _, z3, z4, _, z6 = out
    z = y6 - a
    p1, p2, p3 = tuple(max(0, reflection(g, i) + 1) for i in range(1, 4))
    assert (b - a, b + z1 - a - y6, z3 - y1, z6 - y3, z4 - y4) == (
        min(z, 1), -sign(z), -p1 + p2, 0, p2 - p3)
    f = -2 * z - 1 if z < 0 else max(1, z)
    work = 3 * (b - a) + z1 - y1 + z3 - y3 + z6 - y6
    assert (work, work - (z4 - y4)) == (p2 - f, p3 - f)


def cycles(mapping):
    complete = set()
    for start in mapping:
        path = []
        index = {}
        cursor = start
        while cursor not in complete and cursor not in index:
            index[cursor] = len(path)
            path.append(cursor)
            cursor = mapping[cursor]
        if cursor in index:
            yield tuple(path[index[cursor]:])
        complete.update(path)


def explicit_constant(m, r):
    # K >= sqrt(m)*(3*(2*r+2)*r*(r-1) + r*(16*r+6)/2) + m*r*(2*r+2).
    return m * (3 * (2 * r + 2) * r * (r - 1) + r * (8 * r + 3) + r * (2 * r + 2))


def radius_numerator(q, base):
    m = len(q) + 1
    return max(abs(m * x - (m - i) * base) for i, x in enumerate(q, 1))


def exhaustive(n, base):
    k = n // 2
    m = k + 1
    mapping = {}
    positive = Counter()
    good = set()
    inside = set()
    for seq in combinations_with_replacement(range(base), k):
        q = tuple(reversed(seq))
        word = whole(canonical(q), base)
        assert word == columns(q, base)
        out = differences(word)
        mapping[q] = out
        if q[0]:
            assert out[0]  # For an odd word the middle digit has the outstanding borrow.
        if q[-1]:
            positive[out] += 1
            assert max(abs(a - b) for a, b in zip(out, leading_output(q, base))) <= 1
        g = gaps(q, base)
        if (m + 1) * min(g) >= base + 1:
            assert pairing(q, base)
        if not pairing(q, base):
            continue
        good.add(q)
        d, h, a, _ = update_terms(g)
        if not any(h[2:]):
            inside.add(q)
        moved = permute(g)
        assert gaps(out, base) == paired_update(g)
        integral_work = 2 * sum(x * y for x, y in zip(moved, a)) + sq(a)
        reflected_work = 2 * sum((h[i] - h[i - 1]) * (d[i] + h[i] - d[i - 1] - h[i - 1])
                                 for i in range(1, len(d)))
        assert integral_work == reflected_work <= 0
        if n == 13:
            exact_work(g, gaps(out, base))
    assert max(positive.values(), default=0) <= 6 ** k
    hist = Counter()
    in_chamber = 0
    max_first = 0
    for path in cycles(mapping):
        period = len(path)
        hist[period] += 1
        for i, q in enumerate(path):
            assert (i % period == 0) == (q == path[0])
        max_first = max(max_first, max(q[0] for q in path))
        if all(q in inside for q in path):
            in_chamber += 1
        if n == 13 and any(q not in good for q in path):
            assert 16 * max(radius_numerator(q, base) for q in path) > base - 7
        if n not in (7, 11, 13):
            continue
        phase = max(range(period), key=lambda i: sq(gaps(path[i], base)))
        window = [path[(phase - i) % period] for i in range(1, 2 * k + 1)]
        if all(q in good for q in window):
            h = max(max(update_terms(gaps(q, base))[1][2:]) for q in window)
            peak = max(radius_numerator(q, base) for q in path)
            assert peak <= m * explicit_constant(m, k) * (1 + h)
    print(f"n={n} B={base}: {len(mapping)} words by differences, {len(good)} paired, "
          f"max positive-inner fiber {max(positive.values(), default=0)}, cycles {sorted(hist.items())}, "
          f"in C {in_chamber}, max first q/B {Fraction(max_first, base)}")


def entrance_charging(n, base, head, horizon):
    """Test the accounting on exact weights and literal whole-word histories."""
    point = Counter()
    joint = Counter()
    bad_original = omitted = 0
    autonomous_marks = {}
    count = base ** n - base
    histories = {}
    for word in combinations_with_replacement(range(base), n):
        q0 = differences(word)
        if not q0[0]:
            continue
        weight = factorial(n) // prod(factorial(x) for x in Counter(word).values())
        if q0 not in histories:
            index = {}
            phases = []
            q = q0
            while q not in index:
                index[q] = len(phases)
                phases.append(q)
                q = physical(q, base)
            first = index[q]
            period = len(phases) - first
            entry = phases[head] if head < len(phases) else phases[first + ((head - first) % period)]
            covered = all(pairing(q, base) for q in phases[min(head, first):])
            bad = period == 1
            if covered:
                g0 = gaps(entry, base)
                g = g0
                seen = {}
                for t in range(horizon):
                    if g in seen:
                        assert bad == (t - seen[g] == 1)
                        break
                    seen[g] = t
                    g = paired_update(g)
                else:
                    raise AssertionError("control horizon too short")
            histories[q0] = (entry, covered, bad)
        entry, covered, bad = histories[q0]
        point[entry] += weight
        bad_original += bad * weight
        if covered:
            joint[entry] += weight
            autonomous_marks[entry] = bad
        else:
            omitted += weight
    assert sum(point.values()) == count
    on_bad = sum(weight for entry, weight in joint.items() if autonomous_marks[entry])
    assert on_bad <= bad_original <= omitted + on_bad
    assert sum(joint.values()) + omitted == count
    print(f"entrance n={n} B={base} T={head}: covered {count-omitted}/{count}, "
          f"fixed-period mass {bad_original}/{count}, covered and fixed {on_bad}/{count}")


def boundary(base):
    r = base // 8
    defective = 0
    for cuts in combinations_with_replacement(range(r + 1), 6):
        excess = tuple(b - a for a, b in zip((0,) + cuts, cuts + (r,)))
        g = tuple(r + x for x in excess)
        q = suffixes(g)
        actual = gaps(physical(q, base), base)
        formatted = [base - 1]
        for j in range(6):
            formatted.extend(sorted((q[j] - int(j == 5), base - q[5 - j] - int(j < 5)), reverse=True))
        sorted_whole = sorted(whole(canonical(q), base), reverse=True)
        assert sorted(formatted, reverse=True) == sorted_whole
        if formatted != sorted_whole:
            defective += 1
            assert g[:4] == (r,) * 4
            assert actual[:2] == paired_update(g)[:2]
            assert all(abs(x - y) <= 1 for x, y in zip(formatted, sorted_whole))
        else:
            assert actual == paired_update(g)
    print(f"thirteen-digit boundary control B={base}: {defective} defective pairings")


def orbit_classes(n):
    k = n // 2
    labels = [-1] * (k + 1)
    labels[0] = 1
    sizes = {}
    for start in range(1, k + 1):
        if labels[start] >= 0:
            continue
        size = 0
        cursor = start
        while labels[cursor] < 0:
            labels[cursor] = start
            size += 1
            cursor = folded(cursor, n)
        assert cursor == start
        sizes[start] = size
    return labels, sizes


def prime_power_controls():
    p = 1069
    assert all(p % i for i in range(2, isqrt(p) + 1))
    labels, lens = orbit_classes(p)
    assert lens == {1: 178, 3: 178, 7: 178}
    assert min(i for i in range(1, p) if pow(2, i, p) in (1, p - 1)) == 178
    assert pow(2, 178, p) == p - 1 and pow(2, 356, p * p) == 1 + 417 * p
    rows = {}
    for i in range(2, 268):
        rows[i] = tuple(reflection(tuple(int(label == group) for label in labels), i) for group in (3, 7))
        assert rows[i] != (0, 0)
    assert [rows[i] for i in (2, 3, 31, 4, 34)] == [(1, 0), (2, 0), (-3, 0), (1, 1), (-1, -1)]
    n = p * p
    labels, sizes = orbit_classes(n)
    assert sizes == {1: 178 * p, 3: 178 * p, 7: 178 * p, p: 178, 3 * p: 178, 7 * p: 178}
    row = Counter(labels[-534:])
    row.subtract(labels[:534])
    assert not any(row.values())
    h = n // (2 * p)
    m = n // 2 + 1
    for a, t in ((h + 1, 1), (2 * h + 2, 2)):
        base = m * a + h * t
        g = (a,) + tuple(a + (i % p == 0) * t for i in range(1, m))
        q = suffixes(g)
        assert sum(g) == base and pairing(q, base)
        assert physical(q, base) == q
        print(f"n=1069^2 exact whole-word fixed point B={base}, peak {Fraction(radius_numerator(q, base), m * base)}")
    print("n=1069: three 178-orbits, 266 nonzero secondary rows; its square forces row 534")


def root_update(x):
    z = x[-1]
    e = sign(z)
    h = max(0, -1 - z)
    return (z - 3 * e + 2 * h, x[0] + 2 * h) + tuple(a - e + h for a in x[1:-1]), h


def root_history(source, cap=500000):
    point = source
    path = []
    index = {}
    loss = 0
    while point not in index and len(path) < cap:
        index[point] = (len(path), loss)
        path.append(point)
        point, h = root_update(point)
        loss += h
    if point not in index:
        raise AssertionError("root control exceeded cap")
    first, previous_loss = index[point]
    return path, first, loss, loss - previous_loss


def determinant(rows):
    matrix = [[Fraction(x) for x in row] for row in rows]
    total = Fraction(1)
    for j in range(len(matrix)):
        pivot = next((i for i in range(j, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            total = -total
        divisor = matrix[j][j]
        total *= divisor
        for i in range(j + 1, len(matrix)):
            factor = matrix[i][j] / divisor
            matrix[i] = [a - factor * b for a, b in zip(matrix[i], matrix[j])]
    assert total.denominator == 1
    return total.numerator


def ray_at(path, start, stop):
    before, after = path[start], path[stop]
    r = len(before)
    assert (stop - start) % r == 0
    return tuple(before[(i + start) % r] - after[(i + start) % r] for i in range(r))


def retirement(source, intervals):
    path, first, loss, recurring = root_history(source)
    assert not recurring
    rays = [ray_at(path, start, stop) for start, stop in intervals]
    rank = determinant(rays)
    assert rank and len(rays) == len(source)
    for ray, (_, stop) in zip(rays, intervals):
        assert min(ray) >= 0
        for t in range(stop):
            if ray[(len(source) - 1 - t) % len(source)]:
                assert path[t][-1] > 0
    print(f"order {len(source)} source {source}: preperiod {first}, period {len(path)-first}, "
          f"loss {loss}, insertion determinant {rank}")
    return rays


def check_lift(n, x, base, anchor):
    labels, _ = orbit_classes(n)
    j = 1
    slots = []
    while j not in slots:
        slots.append(j)
        j = folded(j, n)
    assert len(x) == len(slots)
    g = [None] * len(labels)
    g[0] = anchor
    for j, relative in zip(slots, x):
        g[j] = anchor + relative
    complement = [j for j, value in enumerate(g) if value is None]
    assert len(complement) == 1
    g[complement[0]] = base - sum(a for a in g if a is not None)
    g = tuple(g)
    visited = set()
    peak_per_phase = []
    loss = 0
    while g not in visited:
        visited.add(g)
        q = suffixes(g)
        assert pairing(q, base) and not any(update_terms(g)[1][2:])
        out = gaps(physical(q, base), base)
        trial, h = root_update(x)
        assert out == paired_update(g)
        assert tuple(out[j] - out[0] for j in slots) == trial
        assert sum(out[j] for j in [0] + slots) - sum(g[j] for j in [0] + slots) == h
        peak_per_phase.append((g, q[0]))
        loss += h
        g, x = out, trial
    first = next(i for i, (previous, _) in enumerate(peak_per_phase) if previous == g)
    return first, len(visited) - first, loss, min(q1 for _, q1 in peak_per_phase[first:])


def cycle_from(q, base):
    path = []
    at = {}
    while q not in at:
        at[q] = len(path)
        path.append(q)
        q = physical(q, base)
    return at[q], path[at[q]:]


def world_controls():
    # The known 23-digit exterior fixed family satisfies the breach alternative.
    for a in (5, 11, 101):
        base = 71 * a
        centered = (-a, -2*a, -3*a, -4*a+12, a-18, -18, -a-6, -2*a+6, -3*a+18, 2*a-12, a)
        q = tuple(((12-i) * base + x) // 12 for i, x in enumerate(centered, 1))
        assert tuple(12 * x - (12-i)*base for i, x in enumerate(q, 1)) == centered
        assert physical(q, base) == q and pairing(q, base)
        h = max(update_terms(gaps(q, base))[1][2:])
        assert h == max(2, (a-5)//2)
        print(f"23-digit fixed barrier B={base}: H={h}, absolute peak {Fraction(max(map(abs, centered)), 12)}")

    # A centered class-0 trace verifies the exact period and positive base direction.
    tails = []
    for base in (840, 847, 1400):
        x0 = (49, 7, -49, 21, 21, 14)
        q = tuple(((7-i) * base + x) // 7 for i, x in enumerate(x0, 1))
        first, path = cycle_from(q, base)
        assert first == 32 and len(path) == 177 and len(set(path)) == 177
        assert path[3] != path[0] and path[59] != path[0]
        shifted = [tuple(7 * x - (7-i) * base for i, x in enumerate(current, 1)) for current in path]
        tails.append(shifted)
        print(f"13-digit class-0 control B={base}: preperiod {first}, period {len(path)}, "
              f"peak numerator {max(max(map(abs, s)) for s in shifted)}")
    assert tails[0] == tails[1] == tails[2]

    for c in (101, 102, 10001, 10002):
        odd = c % 2
        centered = (c, c, c+14, c, -(5*c + (35 if odd else 28))//2,
                    -(5*c - (7 if odd else 0))//2, c)
        base = 48 * c + 100
        base += (-c - base) % 7
        g = tuple((base + x) // 7 for x in centered)
        assert tuple(7 * x - base for x in g) == centered
        first, path = cycle_from(suffixes(g), base)
        assert first == 0 and len(path) == (4 if odd else 8)
        peak = 0
        for i, q in enumerate(path):
            state = gaps(q, base)
            assert pairing(q, base) and paired_update(state) == gaps(path[(i+1) % len(path)], base)
            peak = max(peak, base - 7 * min(state))
        assert 2 * peak - 5 * c == (35 if odd else 42)
        print(f"13-digit even-period F control B={base}: period {len(path)}, peak R {peak}")

    # P4's endpoint describes a limiting law, as a finite exceptional orbit illustrates.
    for seq in combinations_with_replacement(range(24), 4):
        q = tuple(reversed(seq))
        if q[0] == 22:
            first, path = cycle_from(q, 24)
            if first == 0:
                assert Fraction(q[0], 24) > Fraction(5, 6)
                print(f"9-digit finite physical tail B=24: cyclic q={q}, period {len(path)}")
                break
    else:
        raise AssertionError("missing finite physical-tail control")


def strict_nine():
    source = (16, 7, 17)
    intervals = ((0, 3), (3, 9), (12, 21))
    rays = retirement(source, intervals)
    assert rays == [(4, 4, 4), (5, 0, 6), (5, 0, 0)]
    for residue in range(5):
        base = 120000 + residue
        virtual_base = 6 * base
        a = Fraction(11, 60)
        anchor = (virtual_base - a * base) // 5
        results = []
        for counts in ((0, 0, 0), (20, 0, 20), (20, 20, 20)):
            x = tuple(s + sum(c * ray[j] for c, ray in zip(counts, rays)) for j, s in enumerate(source))
            first, period, loss, terminal_q1 = check_lift(9, x, virtual_base, anchor + residue)
            assert first == 28 + sum(c * (b - a) for c, (a, b) in zip(counts, intervals))
            assert period == 2 and loss == 5 + counts[1] + 3 * counts[2]
            g = (anchor + residue, anchor + residue + x[0], anchor + residue + x[1],
                 virtual_base - 4 * (anchor + residue) - sum(x), anchor + residue + x[2])
            assert virtual_base - 5 * min(g) <= base // 5
            reported = Fraction(4, 5) + (terminal_q1 - Fraction(4 * virtual_base, 5)) / base
            assert reported > Fraction(5, 6) + Fraction(1, 1000)
            results.append(reported)
        print(f"9-digit coset {residue}: three literal lifts, minimum reported terminal q1/B {min(results)}")


def main():
    for n in (7, 9, 11, 13):
        for base in (4, 8, 13) if n == 13 else (4, 8, 13, 24):
            exhaustive(n, base)
    for base in (6, 14, 20):
        exhaustive(13, base)
    for n, base in ((7, 4), (9, 5), (11, 6), (13, 7)):
        entrance_charging(n, base, 2, 200)
    for base in (8, 48, 104):
        boundary(base)
    world_controls()
    prime_power_controls()
    strict_nine()


if __name__ == "__main__":
    main()
