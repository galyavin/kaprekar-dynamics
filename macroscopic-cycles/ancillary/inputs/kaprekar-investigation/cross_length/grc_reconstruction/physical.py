"""Exact physical checks through whole-word subtraction, with no repo imports."""

from collections import Counter
from fractions import Fraction


def gaps(base, differences):
    return tuple(a - b for a, b in zip((base, *differences), (*differences, 0)))


def differences(state):
    return tuple(sum(state[i:]) for i in range(1, len(state)))


def numeral(digits, base):
    result = 0
    for digit in digits:
        result = base * result + digit
    return result


def literal(radii, base):
    """Subtract the two integers made from a valid representative word."""
    forward = (*radii, *((0,) * (len(radii) + 1)))
    value = numeral(forward, base) - numeral(reversed(forward), base)
    unpacked = []
    for _ in forward:
        value, digit = divmod(value, base)
        unpacked.append(digit)
    assert value == 0
    ordered = sorted(unpacked)
    return tuple(ordered[~i] - ordered[i] for i in range(len(radii)))


def margins(state):
    return tuple(sum(state[-i:]) - sum(state[:i]) for i in range(1, len(state) // 2 + 1))


def paired(state, secondary=True):
    """Folded permutation plus exact borrow corrections from the paired formula."""
    n = 2 * len(state) - 1
    result = list(state)
    for i in range(1, len(state)):
        result[min(2 * i, n - 2 * i)] = state[i]
    d = margins(state)
    sign = (d[0] > 0) - (d[0] < 0)
    for i, factor in enumerate((1, -2, 1)):
        result[i] += factor * sign
    for wall, margin in enumerate(d, 1):
        if wall == 1 or secondary:
            correction = max(0, -1 - margin)
            for index, sign in zip(range(2 * wall - 2, 2 * wall + 2), (-1, 1, 1, -1)):
                if index < 2 * (len(state) // 2):
                    result[index] += sign * correction
    return tuple(result)


def certify_primitive(state, base):
    assert (len(state) + 1) * min(state) > base + 1
    exact = gaps(base, literal(differences(state), base))
    assert exact == paired(state)
    return exact


def corner_state(anchor, low, u, v):
    return tuple(anchor + x for x in (0, 0, 0, low + u + v + 4, 1, low + u, low + v, low, -2))


def corner_offsets():
    """Track the 244 phases, separating the permuted unit coefficients."""
    n = 17
    positions = list(range(9))
    offset = corner_state(0, 0, 0, 0)
    coefficients = ((0, 0, 0), (0, 0, 0), (0, 0, 0), (1, 1, 1), (0, 0, 0),
                    (1, 1, 0), (1, 0, 1), (1, 0, 0), (0, 0, 0))
    table = []
    for phase in range(244):
        assert max(map(abs, offset)) <= 15, (phase, offset)
        forms = tuple((*row, constant) for row, constant in zip(
            zip(*(margins(tuple(coefficients[i][axis] for i in positions)) for axis in range(3))),
            margins(offset)))
        assert forms[0][:3] == (0, 0, 0)
        for wall, form in enumerate(forms[1:], 2):
            c, u, v, intercept = form
            assert min(c, u, v) >= 0
            if (wall, phase % 4) in ((2, 0), (4, 0)):
                assert (c, u, v) == (wall // 2, 0, 0)
            else:
                assert u + v >= (2 if wall == 4 else 1)
                assert abs(intercept) <= 30 * wall
        table.append(forms)
        offset = paired(offset, secondary=False)
        moved = positions.copy()
        for i in range(1, len(positions)):
            moved[min(2 * i, n - 2 * i)] = positions[i]
        positions = moved
    assert positions == list(range(9))
    assert tuple(a - b for a, b in zip(offset, corner_state(0, 0, 0, 0))) == (
        12, 12, 12, -15, 12, -15, -15, -15, 12)
    assert min((table[t][1][-1], t) for t in range(0, 244, 4)) == (-29, 228)
    assert min((table[t][3][-1], t) for t in range(0, 244, 4)) == (-57, 240)
    return table


def corner():
    forms = corner_offsets()
    base = 175379
    counts = Counter()
    simultaneous = []
    steps = 0
    for residue in range(27):
        low = base // 350
        low += (28 + residue - low) % 27
        v = base // 96
        u = base // 95
        u += (5 * (base - 4 * low - 2 * v - 3) - u) % 9
        anchor, remainder = divmod(base - 4 * low - 2 * u - 2 * v - 3, 9)
        assert remainder == 0
        assert Fraction(1, 500) < Fraction(low, base) < Fraction(1, 250)
        assert all(Fraction(1, 100) < Fraction(x, base) < Fraction(1, 90) for x in (u, v))
        current = corner_state(anchor, low, u, v)
        copies = (low - 28) // 27 + 1
        stopped = 1 + residue
        assert low - copies * 27 == stopped
        target_phase = None
        for phase, rows in enumerate(forms):
            failures = tuple(i for i, (a, b, c, d) in enumerate(rows[1:], 2)
                             if a * stopped + b * u + c * v + d < -1)
            if failures:
                target_phase = phase
                break
        assert target_phase is not None
        for step in range(copies * 244 + target_phase + 1):
            block, phase = divmod(step, 244)
            if phase == 0:
                assert current == corner_state(anchor + 12 * block, low - 27 * block, u, v)
            values = margins(current)
            assert values == tuple(a * (low - 27 * block) + b * u + c * v + d
                                   for a, b, c, d in forms[phase])
            observed = tuple(i for i, m in enumerate(values[1:], 2) if m < -1)
            successor = certify_primitive(current, base)
            if observed:
                assert block == copies and phase == target_phase and observed == failures
                assert successor != paired(current, secondary=False)
                counts[observed] += 1
                if len(observed) > 1:
                    simultaneous.append((stopped, phase, values, successor))
            else:
                assert successor == paired(current, secondary=False)
            current = successor
            steps += 1
        print(f"corner residue {residue:2d} copies {copies:2d} phase {target_phase:3d} walls {failures}")
    assert counts == {(2,): 19, (4,): 5, (2, 4): 3}
    assert [row[:2] for row in simultaneous] == [(8, 100), (26, 228), (27, 228)]
    print("corner fresh base", base, "literal steps", steps, "counts", dict(counts))
    print("simultaneous stopped offset, phase, margins and successors:", simultaneous)


def wall_start(n, anchor, low, high):
    k = (n - 1) // 2
    state = [anchor + high[i] for i in range(k + 1)]
    for i, shift in zip((0, 1, 2, 4, k, k - 1), (0, 0, 0, 1, -2, low)):
        state[i] = anchor + shift
    return tuple(state)


def wall():
    for n, copies, threshold, decrement in ((15, 2, 22, 21), (17, 4, 29, 27)):
        k = (n - 1) // 2
        anchor = 12000
        high = tuple(811 + 137 * i for i in range(k + 1))
        aligned = wall_start(n, anchor, threshold, high)
        state = aligned
        base = sum(state)
        margins_2 = []
        for _ in range(61 * copies):
            ds = margins(state)
            assert all(x >= 0 for x in ds[1:])
            margins_2.append(ds[1])
            state = certify_primitive(state, base)
        orbit = {k - 1}
        while len(orbit) < copies:
            orbit |= {min(2 * i, n - 2 * i) for i in orbit}
        assert state == tuple(value + (3 * copies if i in (0, 1, 2, 4, k)
                                       else -15 if i in orbit else 0)
                              for i, value in enumerate(aligned))
        assert state[k - 1] - state[0] == threshold - decrement and min(margins_2) == 0
        periods = Counter()
        unresolved = []
        for low in range(threshold - decrement, threshold):
            adjusted = list(high)
            adjusted[3] += threshold - low
            state = wall_start(n, anchor, low, adjusted)
            assert sum(state) == base
            seen = {}
            for step in range(500):
                if state in seen:
                    start = seen[state]
                    period = step - start
                    if period not in ({6, 8} if n == 15 else {8, 12}):
                        unresolved.append(low)
                        break
                    assert all(min(margins(g)[:2]) >= -1 for g, t in seen.items() if t >= start)
                    periods[period] += 1
                    break
                seen[state] = step
                assert all(x >= -1 for x in margins(state)[2:])
                state = certify_primitive(state, base)
            else:
                unresolved.append(low)
        assert (periods, unresolved) == (({6: 5, 8: 16}, []) if n == 15
                                         else ({8: 20, 12: 5}, [2, 16])), (periods, unresolved)
        print("isolated wall", n, "each weight", Fraction(1, decrement),
              "repeat counts", dict(periods), "unresolved", unresolved)


if __name__ == "__main__":
    wall()
    corner()
