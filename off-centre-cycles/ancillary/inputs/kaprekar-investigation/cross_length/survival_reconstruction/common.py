"""Independent integer dynamics, certificate tools, and literal word controls."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import accumulate
from math import ceil, factorial, lcm


def add(a, b, scale=1):
    return tuple(x + scale*y for x, y in zip(a, b))


def rotate(x, phase=1):
    phase %= len(x)
    return x[-phase:] + x[:-phase] if phase else x


def orbit_length(n):
    return (n+1)//2


def folded(n, i):
    return min(2*i, n-2*i)


def cycles(n):
    todo = set(range(orbit_length(n)))
    answer = []
    while todo:
        start = min(todo)
        cycle = []
        i = start
        while i in todo:
            todo.remove(i)
            cycle.append(i)
            i = folded(n, i)
        assert i == start
        answer.append(tuple(cycle))
    return answer


def permute(g):
    out = [0]*len(g)
    for i, value in enumerate(g):
        out[folded(2*len(g)-1, i)] = value
    return tuple(out)


def difference_rows(g):
    return tuple(b-a for _,a,b in zip(range(len(g)//2), accumulate(g), accumulate(reversed(g))))


def loss(z):
    return max(0, -1-z)


def sign(z):
    return (z > 0)-(z < 0)


def root_kick(x, e, h):
    z = x[-1]
    return (z-3*e+2*h, x[0]+2*h, *(v-e+h for v in x[1:-1]))


def root_step(x):
    return root_kick(x, sign(x[-1]), loss(x[-1]))


def gap_step(g, corrections=None):
    out = list(permute(g))
    ds = difference_rows(g)
    e = sign(ds[0])
    out[0] += e
    out[1] -= 2*e
    out[2] += e
    if corrections is None:
        corrections = range(1, len(ds)+1)
    for i in corrections:
        h = loss(ds[i-1])
        for j, factor in enumerate((-1, 1, 1, -1), 2*i-2):
            if j < len(g)-(len(g)%2):
                out[j] += factor*h
    return tuple(out)


def gap_rows(g):
    B = sum(g)
    q = tuple(sum(g[i:]) for i in range(1, len(g)))
    low_column = tuple(B-v-(i < len(q)-1) for i, v in enumerate(reversed(q)))
    high_column = q[:-1]+(q[-1]-1,)
    pairs = tuple(zip(high_column, low_column))
    return B, pairs


def paired_slacks(g):
    B, pairs = gap_rows(g)
    return (B-1-max(pairs[0]), *(min(a)-max(b) for a, b in zip(pairs, pairs[1:])))


def sorted_word(g):
    k = len(g)-1
    left = tuple(max(1, v//2) for v in g[1:])
    right = tuple(v-s for v, s in zip(g[1:], left))
    digits = [g[0]//2]
    for spacing in (*left, *reversed(right)):
        digits.append(digits[-1]+spacing)
    assert len(digits) == 2*k+1 and digits[-1] < sum(g) and len(set(digits)) == len(digits)
    return tuple(digits)


def word_gaps(word, B):
    digits = sorted(word)
    k = len(digits)//2
    q = [B] + [digits[-1-i]-digits[i] for i in range(k)] + [0]
    return tuple(a-b for a, b in zip(q, q[1:]))


def subtract(word, B):
    asc = sorted(word)
    result = 0
    for hi, lo in zip(reversed(asc), asc):
        result = B*result+hi-lo
    digits = []
    for _ in word:
        result, digit = divmod(result, B)
        digits.append(digit)
    assert result == 0
    return tuple(reversed(digits))


def literal(g, count, check=None):
    B = sum(g)
    word = sorted_word(g)
    for i in range(count):
        assert g[-1] > 0 and min(paired_slacks(g)) >= 0
        next_word = subtract(word, B)
        next_g = word_gaps(next_word, B)
        assert next_g == gap_step(g), (i, g, next_g, gap_step(g), difference_rows(g))
        if check:
            check(i, g)
        g, word = next_g, next_word
    return g


def periodic_tail(x, step):
    seen = {}
    states = []
    while x not in seen:
        seen[x] = len(states)
        states.append(x)
        x = step(x)
    return seen[x], states[seen[x]:]


def determinant(columns):
    rows = [list(map(Q, row)) for row in zip(*columns)]
    value = Q(1)
    for i in range(len(rows)):
        p = next((p for p in range(i, len(rows)) if rows[p][i]), None)
        if p is None:
            return 0
        if p != i:
            rows[i], rows[p] = rows[p], rows[i]
            value = -value
        pivot = rows[i][i]
        value *= pivot
        for j in range(i+1, len(rows)):
            factor = rows[j][i]/pivot
            for k in range(i+1, len(rows)):
                rows[j][k] -= factor*rows[i][k]
    return value


def spread(x):
    return max(0, *x)-min(0, *x)


@dataclass
class Route:
    rays: list
    lows: list
    time: int
    loss: int
    spread: int
    terminal: tuple
    period: int
    first_retire: tuple
    lc_budget: tuple


def insertions(source, intervals, end, start=None):
    """Recover guarded affine words from reference endpoint differences."""
    r = len(source)
    points = {t for pair in intervals for t in pair} | {end}
    if start is not None:
        points.add(start)
    snapshots = {}
    first_negative = [None]*r
    x = tuple(source)
    total_loss = amplitude = 0
    cycle_head = None
    period = 0
    for t in range(end+1):
        if t in points:
            snapshots[t] = x, total_loss
        amplitude = max(amplitude, spread(x))
        if start is not None and t == start:
            cycle_head = x
        if start is not None and t > start and not period and x == cycle_head:
            period = t-start
        label = (r-1-t)%r
        if x[-1] <= 0 and first_negative[label] is None:
            first_negative[label] = t
        if t < end:
            total_loss += loss(x[-1])
            x = root_step(x)
    rays = []
    word_costs = []
    prior = 0
    limit = start if start is not None else end
    for s, t in intervals:
        assert prior <= s < t <= limit and (t-s)%r == 0
        prior = t
        ray = rotate(add(snapshots[s][0], snapshots[t][0], -1), -s)
        for label, coefficient in enumerate(ray):
            assert not coefficient or first_negative[label] is None or first_negative[label] >= t
        rays.append(ray)
        word_costs.append((t-s, snapshots[t][1]-snapshots[s][1], spread(ray)))
    det = abs(determinant(rays))
    assert det
    floors = [0]*r
    for label in range(r):
        last = max(i for i, ray in enumerate(rays) if ray[label])
        pivot = rays[last][label]
        assert pivot > 0
    for i, ray in enumerate(rays):
        for label, pivot in enumerate(ray):
            if pivot and not any(later[label] for later in rays[i+1:]):
                deficit = sum((floors[j]+1)*max(0, -rays[j][label]) for j in range(i))
                floors[i] = max(floors[i], ceil(Q(deficit, pivot)))
    budgets = tuple(base+sum((floor+1)*row[i] for floor, row in zip(floors, word_costs))
                    for i, base in enumerate((end, total_loss, amplitude)))
    # With zero floors the sharp coefficient on L is the larger of the
    # reference budget and the sum of costs per additional copy.
    sharp = tuple(max(base, sum(row[i] for row in word_costs))
                  for i, base in enumerate((end, total_loss, amplitude))) if not any(floors) else budgets
    if start is not None:
        assert period and (end-start)%period == 0
        assert snapshots[start][0] == snapshots[end][0]
        assert snapshots[start][1] == snapshots[end][1]
    print(f'r={r} det={det} floors={floors} period={period or "-"} budgets={budgets} sharp={sharp}', flush=True)
    return Route(rays, floors, *sharp, snapshots[end][0], period, tuple(first_negative), budgets)


def lift(n, source, route, overrides=None, mean_scaling=False):
    decomposition = cycles(n)
    (_,), labels, *components = decomposition
    m = orbit_length(n)
    assert len(labels) == len(source)
    assert 3 not in labels and components
    weights = [0]*m
    for part in components:
        for index in part:
            weights[index] = (overrides or {}).get(part[0], 1)
    center = Q(sum(weights), m)
    width = center-min(weights) if mean_scaling else max(weights)-min(weights)
    epsilon = Q(1, 2*m*(m+1))/width
    pi = tuple(Q(1,m)+epsilon*(weight-center) for weight in weights)
    slack = min(pi)-Q(1,m+1)
    mu = min(difference_rows(pi)[1:])
    radius = max(abs(v-Q(i,m)) for i,v in enumerate(accumulate(pi),1))
    room = min(slack, mu/m, radius/m)
    error = m+2*max(route.lc_budget[1:])+1
    D = ceil(4*error/room)
    assert error > route.spread + Q(route.loss,len(labels)+1)+1
    assert Q(error, D) <= room/4 and room > 0
    assert (m+1)*(slack*D-error)>1
    assert D*(min(pi)-Q(1,2*(m+1)))-error >= 1
    # Every factor of the split count eventually exceeds this bound.
    factor = 2*(m+1)
    eta = f'{n}!/({factor}^{m}*{D}^{m-1})'
    assert Q(factorial(n), factor**m*D**(m-1)) > 0
    print(f'n={n} order={len(labels)} margin={mu} delta={radius/2} D={D} eta={eta}', flush=True)
    return labels, components, pi, D, radius/2


def lift_state(B, x, labels, components, pi, offsets=None):
    g = [None]*len(pi)
    free = tuple(i for part in components for i in part)
    offsets = offsets or (0,)*len(free)
    g[0] = (B*pi[0]-Q(sum(x),len(x)+1)).__floor__()+offsets[0]
    for j, v in zip(labels, x):
        g[j] = g[0]+v
    for j, u in zip(free[:-1], offsets[1:]):
        g[j] = (B*pi[j]).__floor__()+u
    g[free[-1]] = B-sum(v for v in g if v is not None)
    return tuple(g)


def physical_period(g, expected):
    B = sum(g)
    word = sorted_word(g)
    first_word = subtract(word, B)
    current = first_word
    length = 0
    while True:
        current = subtract(current, B)
        length += 1
        if current == first_word:
            break
        assert length < expected
    assert length == expected
    return length


def complement_period(components):
    return lcm(*(len(part) for part in components))
