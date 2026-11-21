"""Exact isolated d2 return on the order-four, period-61 drift section."""
from collections import Counter
from math import lcm

ROOT = (0, 1, 2, 4)
DRIFT_ROOT = (0, 0, 1, -2)


def destination(j, n):
    return min(2 * j, n - 2 * j)


def profile(n, anchor, relative, complement):
    k = n // 2
    root = (*ROOT, k)
    state = [anchor] * (k + 1)
    for j, x in zip(root[1:], relative):
        state[j] += x
    for j, x in zip((j for j in range(k + 1) if j not in root), complement):
        state[j] = x
    return tuple(state)


def root(x):
    return tuple(x[j] - x[0] for j in (1, 2, 4, len(x) - 1))


def margin(x, j):
    return sum(x[-j:]) - sum(x[:j])


def step(x, allow_wall):
    n = 2 * len(x) - 1
    z = x[-1] - x[0]
    e = (z > 0) - (z < 0)
    h1 = max(0, -1 - z)
    h2 = max(0, -1 - margin(x, 2)) if allow_wall else 0
    y = [0] * len(x)
    for j, value in enumerate(x):
        y[destination(j, n)] = value
    for j, kick in enumerate((e-h1, -2*e+h1, e+h1-h2, -h1+h2, h2, -h2)):
        y[j] += kick
    return tuple(y), h1, h2


def drift(n):
    k = n // 2
    orbit = [k - 1]
    while destination(orbit[-1], n) != orbit[0]:
        orbit.append(destination(orbit[-1], n))
    period = lcm(61, len(orbit))
    start = profile(n, 0, DRIFT_ROOT, [0] * (k - 4))
    x = start
    advance = []
    for t in range(period):
        successor, h1, h2 = step(x, False)
        advance.append((x, h1))
        x = successor
    assert root(x) == DRIFT_ROOT
    assert len({root(x) for x, _ in advance[:61]}) == 61
    assert sum(h for _, h in advance) == 15 * (period // 61)
    change = tuple(b - a for a, b in zip(start, x))
    assert change == tuple(3 * (period // 61) if j in (*ROOT, k)
                           else -15 if j in orbit else 0 for j in range(k + 1))
    critical = -min(margin(x, 2) for x, _ in advance[::len(orbit)])
    descent = change[0] - change[k - 1]
    bounds = tuple(max(abs(x[j]) for x, _ in advance) for j in range(k + 1))
    return orbit, period, critical, descent, bounds


def kernel(n, c, orbit, cap=4096):
    k = n // 2
    high = tuple(j for j in range(k + 1) if j not in (*ROOT, k, k - 1))
    x = profile(n, 0, DRIFT_ROOT, [c if j == k - 1 else 100000 + j
                                  for j in sorted((*high, k - 1))])
    high_position = dict(zip(high, high))
    initial = x
    trace = {}
    moves = []
    bounds = [0, 0, 0]
    while (x, tuple(high_position.values())) not in trace and len(trace) < cap:
        trace[x, tuple(high_position.values())] = len(trace)
        low = orbit[len(moves) % len(orbit)]
        assert low not in high_position.values()
        bounds[0] = max(bounds[0], *(abs(x[j]) for j in (*ROOT, k)))
        bounds[1] = max(bounds[1], abs(x[low]))
        bounds[2] = max(bounds[2], *(abs(x[j] - initial[i]) for i, j in high_position.items()))
        y, h1, h2 = step(x, low == k - 1)
        moves.append((h1, h2, low, margin(x, 2), root(x)))
        high_position = {i: destination(j, n) for i, j in high_position.items()}
        x = y
    start = trace.get((x, tuple(high_position.values())))
    if start is None:
        return None, None, tuple(bounds), None
    length = len(moves) - start
    # Full repeat includes every independent high ray in its initial slot.
    assert all(m == 0 and w == 0 for m, w, *_ in moves[start:])
    contacts = tuple((t, d, x) for t, (_, h2, _, d, x) in enumerate(moves) if h2)
    return start, length, tuple(bounds), contacts


def quotient_drift(n, c, orbit):
    """The non-closing quotient retains the low gap relative to the anchor."""
    k = n // 2
    high = tuple(j for j in range(k + 1) if j not in (*ROOT, k, k - 1))
    x = profile(n, 0, DRIFT_ROOT, [c if j == k - 1 else 100000 + j
                                  for j in sorted((*high, k - 1))])
    seen = {}
    used_h1 = used_h2 = 0
    for t in range(1024):
        phase = t % len(orbit)
        low = orbit[phase]
        state = phase, root(x), x[low]-x[0]
        if state in seen:
            time, past, h1, h2 = seen[state]
            return time, t-time, tuple(a-b for a,b in zip(x,past)), used_h1-h1, used_h2-h2
        seen[state] = t, x, used_h1, used_h2
        x, h1, h2 = step(x, low == k - 1)
        used_h1 += h1
        used_h2 += h2
    raise AssertionError('quotient repeat not found')


def main():
    expected = {
        15: Counter({6: 5, 8: 16}),
        17: Counter({8: 20, 12: 5, None: 2}),
    }
    for n in (15, 17):
        orbit, period, critical, descent, bounds = drift(n)
        print(f'{n=} {orbit=} {period=} {critical=} {descent=} drift absolute offsets={bounds}')
        counts = Counter()
        for c in range(critical - descent, critical):
            pre, period, bounds, contacts = kernel(n, c, orbit)
            counts[period] += 1
            if pre is None:
                trans = quotient_drift(n,c,orbit)
                wanted = {
                    2: (29,108,(4,4,4,4,4,-12,-9,-3,4),20,9),
                    16: (122,108,(4,4,4,-3,4,-9,4,-12,4),20,9),
                }
                assert trans == wanted[c]
                print(f'  {c=} unresolved rolling section {trans=}')
            else:
                print(f'  {c=} {pre=} {period=} {bounds=} outside_updates={len(contacts)}')
        assert counts == expected[n]
        print('  resolved zero-loss contacts by return length:', sorted(counts.items(), key=str))


if __name__ == '__main__':
    main()
