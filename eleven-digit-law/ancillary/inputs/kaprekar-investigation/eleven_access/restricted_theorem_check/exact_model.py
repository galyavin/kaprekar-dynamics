"""Exact arithmetic for independent guarded paths in ordinary difference coordinates."""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import gcd, lcm

ZERO = (0,) * 5
IDENTITY = tuple(tuple(int(i == j) for j in range(7)) for i in range(5))
BASE = (0,) * 5 + (1, 0)
ONE = (0,) * 6 + (1,)
CENTRE = tuple(Q(i, 6) for i in (5, 4, 3, 2, 1))


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def combine(a, u, b, v):
    return tuple(a * x + b * y for x, y in zip(u, v))


def total(f, q, base):
    return dot(f, (*q, base, 1))


def leading(f, q):
    return dot(f, (*q, 1))


def bezout(a, b):
    old, current = a, b
    x0, x1, y0, y1 = 1, 0, 0, 1
    while current:
        q, r = divmod(old, current)
        old, current = current, r
        x0, x1, y0, y1 = x1, x0 - q * x1, y1, y0 - q * y1
    sign = 1 if old >= 0 else -1
    return sign * old, sign * x0, sign * y0


def pivot(columns, row):
    columns = list(columns)
    j = next((j for j, col in enumerate(columns) if dot(row, col)), None)
    if j is None:
        return 0, None, tuple(columns)
    first = columns.pop(j)
    for j, other in enumerate(columns):
        a, b = dot(row, first), dot(row, other)
        if b:
            d, x, y = bezout(a, b)
            first, columns[j] = combine(x, first, y, other), combine(-b // d, first, a // d, other)
    value = dot(row, first)
    if value < 0:
        first = tuple(-x for x in first)
    return abs(value), first, tuple(columns)


def canonical(basis, *offsets):
    columns = []
    for axis in range(5):
        row = IDENTITY[axis][:5]
        value, first, basis = pivot(basis, row)
        if not value:
            continue
        offsets = tuple(combine(1, p, -(p[axis] // value), first) for p in offsets)
        columns = [combine(1, c, -(c[axis] // value), first) for c in columns]
        columns.append(first)
    return tuple(columns), offsets


@dataclass(frozen=True)
class Sheet:
    base: int
    stride: int
    point: tuple
    slope: tuple
    basis: tuple

    def reduced(self):
        basis, (point, slope) = canonical(self.basis, self.point, self.slope)
        return Sheet(self.base, self.stride, point, slope, basis)

    def push(self, word):
        return Sheet(self.base, self.stride,
                     tuple(total(f, self.point, self.base) for f in word),
                     tuple(dot(f, (*self.slope, self.stride)) for f in word),
                     tuple(tuple(dot(f, col) for f in word) for col in self.basis)).reduced()

    def exits(self, guard, decrement, drift):
        row = guard[:5]
        value, first, rest = pivot(self.basis, row)
        delta = dot(guard, (*self.slope, self.stride))
        divisor = gcd(decrement, value)
        refine = divisor // gcd(divisor, delta)
        index = decrement // divisor
        advance = tuple(refine * x for x in self.slope)
        if index > 1:
            move = -(refine * delta // divisor) * pow(value // divisor, -1, index) % index
            advance = combine(1, advance, move, first)
        stride = self.stride * refine
        assert dot(guard, (*advance, stride)) % decrement == 0
        out_slope = combine(1, advance, dot(guard, (*advance, stride)) // decrement, drift)
        kernel = rest if first is None else (*rest, tuple(index * x for x in first))
        assert all(dot(row, c) % decrement == 0 for c in kernel)
        out_basis = tuple(combine(1, c, dot(row, c) // decrement, drift) for c in kernel)
        for k in range(refine):
            base = self.base + k * self.stride
            start = combine(1, self.point, k, self.slope)
            for j in range(index):
                point = combine(1, start, j, first) if j else start
                count, remainder = divmod(total(guard, point, base), decrement)
                out = combine(1, point, 1 + count, drift)
                yield remainder, index, Sheet(base, stride, out, out_slope, out_basis).reduced()


def sort_step(state, sign):
    """Form the borrowed columns from right-to-left subtraction, then sort."""
    assert sign(state[-1]) > 0, 'this verification starts at positive height'
    digits = [combine(1, f, -(i == 4), ONE) for i, f in enumerate(state)]
    digits += [combine(1, BASE, -1, ONE)]
    digits += [combine(-1, f, 1, combine(1, BASE, -(i > 0), ONE))
               for i, f in reversed(tuple(enumerate(state)))]
    from functools import cmp_to_key

    def compare(i, j):
        score = sign(combine(1, digits[i], -1, digits[j]))
        return -score if score else i - j

    order = sorted(range(11), key=cmp_to_key(compare))
    guards = [combine(1, state[-1], -1, ONE)]
    guards.extend(combine(1, combine(1, digits[a], -1, digits[b]), -(a > b), ONE)
                  for a, b in zip(order, order[1:]))
    paired = True
    for level in range(4):
        a = max((level, 6 + level), key=order.index)
        b = min((level + 1, 7 + level), key=order.index)
        gap = combine(1, digits[a], -1, digits[b])
        allowed = sign(gap) >= 0
        paired = paired and allowed
        guards.append(gap if allowed else combine(-1, gap, -1, ONE))
    out = tuple(combine(1, digits[a], -1, digits[b]) for a, b in zip(order[:5], order[::-1]))
    return out, guards, paired, tuple(order)


def pull(f, chart):
    return tuple(sum(f[i] * chart[i][j] for i in range(5)) + (f[5] if j == 5 else 0)
                 for j in range(6))


def make_at(chart, point):
    return tuple(leading(f, point) for f in chart)


def stage(chart, point, sheet, horizon=120, skip_limit=250):
    """Return one guarded translation and its alternatives, together with its short prefix."""
    trial_guards = set()
    for shorts in range(skip_limit + 1):
        ref = make_at(chart, point)
        den = lcm(*(x.denominator for x in ref))
        lead_point = tuple(int(den * x) for x in ref) + (den,)
        exact_point = (*sheet.point, sheet.base, 1)

        def sign(f):
            value = dot(f, lead_point) or dot(f, exact_point)
            return (value > 0) - (value < 0)

        state = IDENTITY
        used = []
        first = None
        losses = 0
        for length in range(1, horizon + 1):
            state, guards, paired, _ = sort_step(state, sign)
            losses += not paired
            if first is None:
                first = state
            for guard in guards:
                form = pull(guard, chart)
                if any(form):
                    if leading(form, point) <= 0:
                        return trial_guards, None
                    trial_guards.add(form)
                else:
                    assert all(not dot(guard, c) for c in sheet.basis)
                    assert not dot(guard, (*sheet.slope, sheet.stride))
                    assert dot(guard, exact_point) >= 0
                used.append((guard, form))
            if tuple(pull(f, chart) for f in state) != chart:
                relaxed = tuple(dot(f, lead_point) for f in state)
                diffs = (*lead_point[:5], den)
                before = (den - diffs[0], *(a - b for a, b in zip(diffs[:4], diffs[1:5])), diffs[4])
                after = (den - relaxed[0], *(a - b for a, b in zip(relaxed, relaxed[1:])), relaxed[4])
                if sum(g * g for g in after) < sum(g * g for g in before):
                    break
                continue
            drift = tuple(dot(f, exact_point) - old for f, old in zip(state, sheet.point))
            if not any(drift) or tuple(dot(f, drift) for f in state) != drift:
                continue
            assert all(tuple(dot(f, c) for f in state) == c for c in sheet.basis)
            falling = {}
            for g, form in used:
                s = -dot(g, drift)
                if s <= 0:
                    continue
                timed = tuple(Q(a, s) for a in form)
                if timed in falling:
                    h, t = falling[timed]
                    if Q(dot(g, exact_point), s) >= Q(dot(h, exact_point), t):
                        continue
                falling[timed] = (g, s)
            if falling and min(leading(f, point) for f in falling) > 0:
                return trial_guards, (shorts, length, losses, chart, sheet, drift, falling)
        chart = tuple(pull(f, chart) for f in first)
        sheet = sheet.push(first)
    raise RuntimeError((chart, sheet, ref, 'search budget'))


def exit_chart(chart, drift, time):
    return tuple(combine(1, c, d, time) for c, d in zip(chart, drift))


def centered(q):
    return tuple(6 * (x - y) for x, y in zip(q, CENTRE))


def radius(q):
    x = (0, *centered(q), 0)
    return max(b - a for a, b in zip(x, x[1:]))
