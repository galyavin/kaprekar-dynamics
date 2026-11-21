"""Exact finite controls for Task 217, independent of the earlier lens scripts."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm


def differences(gaps):
    return tuple(sum(gaps[i:]) for i in range(1, len(gaps)))


def gaps_of(base, q):
    return tuple(x - y for x, y in zip((base, *q), (*q, 0)))


def literal(base, q):
    """Subtract synthetic word digits column by column, from the units end."""
    borrow = 0
    digits = []
    for difference in (*(-x for x in q), 0, *reversed(q)):
        digit = difference - borrow
        borrow = int(digit < 0)
        digits.append(digit + base * borrow)
    assert borrow == 0
    ranked = sorted(digits, reverse=True)
    return tuple(digits), tuple(
        ranked[i] - ranked[-i - 1] for i in range(len(q))
    )


def levels(base, q):
    return tuple(
        (x - int(i == len(q) - 1), base - y - int(i < len(q) - 1))
        for i, (x, y) in enumerate(zip(q, reversed(q)))
    )


def paired(base, q):
    pairs = levels(base, q)
    return q[-1] > 0 and all(
        min(a) >= max(b) for a, b in zip(pairs, pairs[1:])
    )


def chamber(base, g):
    return (
        g[0] >= 1 and min(g) >= 0 and sum(g) == base
        and paired(base, differences(g))
        and all(sum(g[-i:]) - sum(g[:i]) >= -1
                for i in range(2, len(g) // 2 + 1))
    )


def folded_slots(n):
    return [0] + [min(2 * j, n - 2 * j) for j in range(1, (n + 1) // 2)]


def orbit(slots, index):
    result = [index]
    while slots[result[-1]] != index:
        result.append(slots[result[-1]])
    return result


def permute(g, slots):
    result = [0] * len(g)
    for i, destination in enumerate(slots):
        result[destination] = g[i]
    return tuple(result)


def sign(value):
    return (value > 0) - (value < 0)


def root_move(x):
    y = x[-1]
    e, h = sign(y), max(0, -1 - y)
    return (y - 3 * e + 2 * h, x[0] + 2 * h,
            *(v - e + h for v in x[1:-1]))


def nominal_move(g, slots):
    e, h = sign(g[-1] - g[0]), max(0, g[0] - g[-1] - 1)
    out = list(permute(g, slots))
    for i, change in enumerate((e - h, -2 * e + h, e + h, -h)):
        out[i] += change
    return tuple(out)


def discrepancy_numerator(base, atoms):
    k = len(atoms)
    return max(
        abs(rank * base - k * x)
        for i, x in enumerate(sorted(atoms))
        for rank in (i, i + 1)
    )


def allowance_count(base, left, right):
    return max(
        sum(a <= t for a in source) - sum(b <= t + 2 for b in target)
        for source, target in ((left, right), (right, left))
        for t in range(-2, base + 1)
    )


def radius_numerator(base, g):
    m = len(g)
    return max(abs(m * qi - (m - i) * base)
               for i, qi in enumerate(differences(g), 1))


def check_ranks_and_steps():
    states = pairs = supported = 0
    for k in range(2, 6):
        slots = folded_slots(2 * k + 1)
        support = {0, *orbit(slots, 1)}
        if k >= 3:
            support.update(orbit(slots, 3))
        for base in range(2, 10):
            for cuts in combinations(range(base + k - 1), k):
                g = (cuts[0] + 1,
                     *(j - i - 1 for i, j in zip(cuts, cuts[1:])),
                     base + k - 2 - cuts[-1])
                q = differences(g)
                digits, successor = literal(base, q)
                z = q.count(0)
                b = allowance_count(base, q, tuple(base - x for x in q))
                folded = tuple(abs(2 * x - base) for x in q)
                assert allowance_count(base, successor, folded) <= z + 1
                assert 2 * discrepancy_numerator(base, q) <= (
                    discrepancy_numerator(base, successor)
                    + base * (b + z + 1) + 6 * k
                )
                states += 1
                if paired(base, q):
                    assert b <= 1
                    assert sorted(digits) == sorted(
                        [base - 1] + [v for pair in levels(base, q) for v in pair]
                    )
                    pairs += 1
                if k >= 3 and chamber(base, g):
                    rotated = permute(g, slots)
                    actual = gaps_of(base, successor)
                    assert nominal_move(g, slots) == actual
                    assert all(actual[j] == rotated[j] for j in range(k + 1)
                               if j not in support)
                    supported += 1
    print(f"Literal rank states: {states}, paired: {pairs}, chamber updates: {supported}")


def check_side_words():
    for n, expected_chamber in ((7, 40), (9, 100), (11, 240)):
        m = (n + 1) // 2
        base = 6 * n + 7
        accepted_pair = accepted_chamber = 0
        for word in product((0, 1), repeat=n):
            left = total = 0
            pass_pair = pass_chamber = True
            auxiliary = []
            for i, side in enumerate(word, 1):
                auxiliary.append(3 * i if side == 0 else base - 3 * i)
                left += int(side == 0)
                total += 1
                a, b = left // 2, (total - left) // 2
                pass_pair &= abs(a - b) <= 1
                pass_chamber &= not any(a < j <= b for j in range(2, m - 1))
            auxiliary.sort()
            q = tuple(base - auxiliary[2 * i - 1] for i in range(1, m))
            assert pass_pair == paired(base, q)
            assert (pass_pair and pass_chamber) == chamber(base, gaps_of(base, q))
            accepted_pair += pass_pair
            accepted_chamber += pass_pair and pass_chamber
        assert Fraction(accepted_pair, 2**n) == Fraction(5, 8)**((n - 1)//4)
        assert accepted_chamber == expected_chamber
        print(f"Folded side words at {n}: paired {accepted_pair}, chamber {accepted_chamber}")


def check_fixed_families_and_siblings():
    for n in (15, 63, 255, 1023):
        m = (n + 1) // 2
        v = (n - 3) // 6
        g = tuple(4 + 2 * int(j > 0 and j % 3 == 0) for j in range(m))
        base = sum(g)
        q = differences(g)
        assert paired(base, q) and literal(base, q)[1] == q
        assert Fraction(q[0], base) - Fraction(m - 1, m) == Fraction(v, m*(2*m+v))
        s = n // 3
        unpaired = (base // 2,) * s + (0,) * ((s - 1) // 2)
        assert not paired(base, unpaired) and literal(base, unpaired)[1] == unpaired
    for n in (9, 15, 17, 31, 33, 63, 65):
        m = (n + 1) // 2
        root = {0, *orbit(folded_slots(n), 1)}
        g = tuple(10 * m + 5 * int(j not in root) for j in range(m))
        base = sum(g)
        q = differences(g)
        assert chamber(base, g)
        y = base - 1 - q[-2]
        assert q[1] > y > q[2]
        sibling = tuple(sorted((*q[:-2], y, q[-1]), reverse=True))
        first, first_gap = literal(base, q)
        second, second_gap = literal(base, sibling)
        assert sorted(first) == sorted(second) and first_gap == second_gap == q
        assert not paired(base, sibling)
    print("Literal fixed-gap pairs: 4. Literal reflected siblings: 7.")


def check_small_time_geometry():
    hexagon = ((-1, -1), (0, -1), (1, 0), (1, 1), (0, 1), (-1, 0))
    twice_area = sum(a*d-b*c for (a, b), (c, d) in
                     zip(hexagon, (*hexagon[1:], hexagon[0])))
    assert twice_area == 6
    assert Fraction(twice_area, 2) / 3 * Fraction(40, 9) * Fraction(3, 2)**2 == 10

    counts = []
    for epsilon in (Fraction(1, 1024), Fraction(1, 2**16), Fraction(1, 2**24)):
        box = list(product(
            (Fraction(1, 2) - 2*epsilon, Fraction(1, 2) - epsilon),
            (epsilon, Fraction(5, 4)*epsilon),
            (epsilon, Fraction(5, 4)*epsilon),
        ))
        previous_max = 0
        count = 0
        while epsilon * 128 * (1 << count) <= 1:
            assert min(v for _, v, _ in box) > previous_max
            previous_max = max(v for _, v, _ in box)
            for a, v, w in box:
                b, c = a - v, 1 - 2*a + v + w
                assert 0 < c < b < a < Fraction(1, 2) and 2*a + b > 1 and w > 0
                assert tuple(sorted((a, abs(a + b - 1), abs(b + c - 1)),
                                    reverse=True)) == (a, a-w, 1-a-b)
                x, y, z = (1+a)/2, (1+c)/2, (1-b)/2
                assert min(1-x, x-y, y-z, z) > Fraction(1, 8)
            count += 1
            box = [(a, v+w, v) for a, v, w in box]
        counts.append(count)
        assert counts[-1] > 0
    print(f"Five-digit leading coefficient: 10. Seven-digit disjoint shells: {counts}.")


def history(x, move, cap):
    positions = {}
    trajectory = []
    while x not in positions and len(trajectory) < cap:
        positions[x] = len(trajectory)
        trajectory.append(x)
        x = move(x)
    return trajectory, positions.get(x)


def check_root_and_guard():
    cycles = closures = surviving = 0
    for r in range(3, 7):
        for x in product(range(-1, 2), repeat=r):
            trajectory = [x]
            for _ in range(2*r):
                trajectory.append(root_move(trajectory[-1]))
            reads = [state[-1] for state in trajectory]
            for t in range(r + 1):
                values = reads[t:t+r]
                losses = [max(0, -1 - a) for a in values]
                signs = list(map(sign, values))
                assert reads[t+r] == (
                    reads[t] - 3*signs[0] + 2*sum(losses[:2])
                    + sum(h-e for h, e in zip(losses[2:], signs[2:]))
                )
            states, first = history(x, root_move, 100)
            if first is None:
                continue
            cycle = states[first:]
            if min(state[-1] for state in cycle) < -1:
                continue
            d = gcd(r, len(cycle))
            assert all(sum(sign(state[-1]) for state in cycle[i::d]) == 0
                       for i in range(d))
            cycles += 1

    slots = folded_slots(9)
    root = orbit(slots, 1)
    support = {0, *root, *orbit(slots, 3)}

    def check_entry(base, g):
        x = tuple(g[i]-g[0] for i in root)
        states, first = history(x, root_move, 120)
        if first is None or min(z[-1] for z in states[first:]) < -1:
            return False, None
        nominal = g
        passes = True
        certificate = first + lcm(len(states) - first, len(root))
        for _ in range(certificate):
            passes &= chamber(base, nominal)
            nominal = nominal_move(nominal, slots)

        literal_g = g
        actual_history = {}
        while literal_g not in actual_history and chamber(base, literal_g):
            actual_history[literal_g] = len(actual_history)
            literal_g = gaps_of(base, literal(base, differences(literal_g))[1])
        actual_pass = literal_g in actual_history
        assert actual_pass == passes
        if passes:
            assert actual_history[literal_g] == first
            span = g
            predicted = 0
            for _ in root:
                predicted = max(predicted, radius_numerator(base, span))
                span = permute(span, slots)
            peak = max(radius_numerator(base, state)
                       for state, t in actual_history.items() if t >= first)
            assert abs(peak - predicted) <= len(g) * sum(g[j] for j in support)
        return True, passes

    for base in (9, 13, 17):
        for cuts in combinations(range(1, base), 4):
            g = tuple(b-a for a, b in zip((0, *cuts), (*cuts, base)))
            if chamber(base, g):
                closed, passes = check_entry(base, g)
                closures += closed
                surviving += bool(passes)

    assert check_entry(7, (2, 1, 1, 2, 1)) == (True, False)
    assert check_entry(58, (12, 11, 11, 13, 11)) == (True, True)
    for s in range(5):
        g = (10+s, 10+s, 10+s, 10+10*s, 10+s)
        base = sum(g)
        assert all(g[j] == g[0] for j in root) and nominal_move(g, slots) == g
        pairs = levels(base, differences(g))
        assert pairs[1][1] - pairs[2][0] == 9 - 8*s
        assert chamber(base, g) == (s < 2)
    print(f"Loss-free root cycle visits: {cycles}. Closed physical entries: {closures}, "
          f"guard-surviving: {surviving}. Translation controls passed.")


if __name__ == "__main__":
    check_ranks_and_steps()
    check_side_words()
    check_fixed_families_and_siblings()
    check_small_time_geometry()
    check_root_and_guard()
