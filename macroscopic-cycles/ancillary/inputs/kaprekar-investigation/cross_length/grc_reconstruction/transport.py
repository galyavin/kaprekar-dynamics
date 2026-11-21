"""Independent affine digit compiler and finite-quotient lattice transport."""

from collections import deque
from dataclasses import dataclass, replace
from fractions import Fraction as Q
from functools import cache, lru_cache
from math import gcd

from physical import literal

replay = lru_cache(maxsize=32768)(literal)


def sum_rows(*rows):
    return tuple(map(sum, zip(*rows)))


def times(scalar, row):
    return tuple(scalar * x for x in row)


def minus(left, right):
    return tuple(a - b for a, b in zip(left, right))


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


@cache
def identity(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def linear(rows, vector):
    return tuple(dot(row, vector) for row in rows)


def multiply(rows, matrix):
    return tuple(linear(tuple(zip(*matrix)), row) for row in rows)


def specialize(row, point, base):
    return dot(row, (*point, base, 1))


@cache
def borrowed_rows(height, k):
    axes = identity(k + 2)
    B, one = axes[-2:]
    if not height:
        return ((0,) * (k + 2),) * (2 * k + 1)
    left = (*axes[:height - 1], minus(axes[height - 1], one))
    middle = (minus(B, one),) * (2 * k + 1 - 2 * height)
    right = tuple(minus(minus(B, axes[i]), one) for i in reversed(range(1, height)))
    return (*left, *middle, *right, minus(B, axes[0]))


@cache
def primitive(height, order, k):
    """On a stable digit sort, return affine outputs and all branch inequalities."""
    digits = borrowed_rows(height, k)
    axes = identity(k + 2)
    guards = []
    if height:
        guards.append(minus(axes[height - 1], axes[-1]))
    if height < k:
        guards.extend(times(-1, row) for row in axes[height:k])
    for a, b in zip(order, order[1:]):
        slack = minus(digits[a], digits[b])
        if a > b:
            slack = minus(slack, axes[-1])
        guards.append(slack)
    return tuple(minus(digits[order[j]], digits[order[~j]]) for j in range(k)), tuple(guards)


def whole_word(point, base, length):
    """Compile one exact branch using the literal integer point as a witness."""
    k = len(point)
    path = []
    itinerary = []
    for _ in range(length):
        height = sum(q > 0 for q in point)
        digits = linear(borrowed_rows(height, k), (*point, base, 1))
        order = tuple(sorted(range(2 * k + 1), key=lambda j: (-digits[j], j)))
        itinerary.append((height, order))
        point = tuple(digits[order[j]] - digits[order[~j]] for j in range(k))
        path.append(point)
    composed, guards = compile_word(tuple(itinerary), k)
    return composed, guards, tuple(path)


@cache
def compile_word(itinerary, k):
    composed = identity(k + 2)
    guards = set()
    for height, order in itinerary:
        rows, inequalities = primitive(height, order, k)
        guards.update(multiply(inequalities, composed))
        composed = (*multiply(rows, composed), *identity(k + 2)[-2:])
    return composed, tuple(sorted(guards))


def support(increments, modulus):
    """BFS in a cyclic group, with an explicit lift of each residue."""
    units = identity(len(increments))
    zero = (0,) * len(increments)
    representatives = {0: zero}
    queue = deque([0])
    while queue:
        residue = queue.popleft()
        for step, unit in zip(increments, units):
            next_residue = (residue + step) % modulus
            if next_residue not in representatives:
                representatives[next_residue] = sum_rows(representatives[residue], unit)
                queue.append(next_residue)
    return representatives


def translate(columns, vector):
    return tuple(dot(row, vector) for row in zip(*columns))


@dataclass(frozen=True)
class Sheet:
    point: tuple
    columns: tuple
    base: int
    period: int
    base_shift: tuple
    macroscopic: tuple
    clock: tuple
    weight: Q = Q(1)

    @property
    def k(self):
        return len(self.point)

    def push(self, word):
        moved = linear(word, (*self.point, self.base, 1))
        directions = tuple(linear(word, (*column, 0, 0))[:self.k] for column in self.columns)
        base_shift = linear(word, (*self.base_shift, self.period, 0))[:self.k]
        return replace(self, point=moved[:self.k], columns=directions, base_shift=base_shift,
                       macroscopic=multiply(word, self.macroscopic))


def initial(base, period, leading):
    k = len(leading)
    point = tuple(int(base * x) for x in leading)
    return Sheet(point, identity(k), base, period, (0,) * k, identity(k + 2), (0,) * (k + 2))


def branch(sheet, length, expected_drift):
    point = sheet.point
    for _ in range(length):
        point = replay(point, sheet.base)
    if minus(point, sheet.point) != expected_drift:
        return None
    word, guards, path = whole_word(sheet.point, sheet.base, length)
    drift = minus(path[-1], sheet.point)
    assert drift == expected_drift
    departures = tuple(minus(row, unit) for row, unit in zip(word, identity(sheet.k + 2)))
    tests = ([(*column, 0, 0) for column in sheet.columns]
             + [(*sheet.base_shift, sheet.period, 0), (*drift, 0, 0)])
    if any(any(linear(departures, test)) for test in tests):
        return None
    return word, guards, drift


@cache
def zero_leading_guards(macroscopic, guards, reference):
    tied = []
    for row in guards:
        pullback = multiply((row,), macroscopic)[0]
        leading = dot(pullback[:-1], (*reference, 1))
        assert leading >= 0
        if not leading:
            assert not any(pullback[:-1]), (row, pullback, reference)
            tied.append(row)
    return tuple(tied)


def check_guards(sheet, guards, reference):
    """A witness gives an open original neighborhood, including structural ties."""
    assert all(specialize(row, sheet.point, sheet.base) >= 0 for row in guards)
    for row in zero_leading_guards(sheet.macroscopic, guards, reference):
        assert not any(dot(row, (*v, 0, 0)) for v in sheet.columns)
        assert dot(row, (*sheet.base_shift, sheet.period, 0)) == 0


def bounded_prefix(sheet, length, drift, maximum, reference):
    for steps in range(maximum + 1):
        compiled = branch(sheet, length, drift)
        if compiled:
            return sheet, compiled, steps
        word, guards, path = whole_word(sheet.point, sheet.base, 1)
        assert path[-1] == literal(sheet.point, sheet.base)
        check_guards(sheet, guards, reference)
        sheet = sheet.push(word)
    raise AssertionError(("return never found", drift, sheet.point))


def falling_rows(sheet, guards, drift):
    candidates = []
    for row in guards:
        speed = -dot(row, (*drift, 0, 0))
        if speed > 0:
            pulled = multiply((row,), sheet.macroscopic)[0]
            form = tuple(Q(x, speed) for x in (*pulled[:-1], 0))
            exact = Q(specialize(row, sheet.point, sheet.base), speed)
            candidates.append((exact, form, row, speed))
    return sorted(candidates)


def exit_children(sheet, compiled, length, reference):
    _, guards, drift = compiled
    check_guards(sheet, guards, reference)
    race = falling_rows(sheet, guards, drift)
    _, leading, guard, decrement = race[0]
    winning = dot(leading, (*reference, 1, 0))
    for _, competitor, _, _ in race[1:]:
        other = dot(competitor, (*reference, 1, 0))
        assert other > winning or competitor == leading, (guard, competitor, leading, other, winning)
    assert winning > 0
    increments = tuple(dot(guard, (*column, 0, 0)) for column in sheet.columns)
    group = support(increments, decrement)
    step = gcd(decrement, *increments)
    assert len(group) == decrement // step
    lift_step = group[step % decrement]
    # The generators below replace integer normal form. For z in the kernel,
    # h(z)/step is a multiple of decrement/step, so these span every such z.
    kernel = [minus(unit, times(increment // step, lift_step))
              for unit, increment in zip(identity(len(increments)), increments)]
    kernel.append(times(decrement // step, lift_step))
    assert all(dot(increments, vector) % decrement == 0 for vector in kernel)
    base_increment = dot(guard, (*sheet.base_shift, sheet.period, 0))
    assert -base_increment % decrement in group, (sheet.period, base_increment, step)
    base_adjustment = translate(sheet.columns, group[-base_increment % decrement])
    adjusted_base = sum_rows(sheet.base_shift, base_adjustment)

    def outgoing(direction, base_delta=0):
        count, residue = divmod(dot(guard, (*direction, base_delta, 0)), decrement)
        assert residue == 0
        return sum_rows(direction, times(count, drift))

    columns = tuple(dict.fromkeys(filter(any, (outgoing(translate(sheet.columns, vector))
                                              for vector in kernel))))
    if not columns:
        columns = ((0,) * sheet.k,)
    base_shift = outgoing(adjusted_base, sheet.period)
    macro = tuple(sum_rows(row, times(d, leading))
                  for row, d in zip(sheet.macroscopic, (*drift, 0, 0)))
    elapsed = sum_rows(sheet.clock, times(length, leading))
    representative_guard = specialize(guard, sheet.point, sheet.base)
    for shift, lift in sorted(group.items()):
        point = sum_rows(sheet.point, translate(sheet.columns, lift))
        value = specialize(guard, point, sheet.base)
        residue = value % decrement
        assert residue == (representative_guard + shift) % decrement
        count = 1 + value // decrement
        assert count >= 1
        assert min(specialize(row, point, sheet.base) for row in guards) >= 0
        last = sum_rows(point, times(count - 1, drift))
        after = sum_rows(last, drift)
        assert min(specialize(row, last, sheet.base) for row in guards) >= 0
        assert min(specialize(row, after, sheet.base) for row in guards) < 0
        yield residue, len(group), replace(
            sheet, point=after, columns=columns, base_shift=base_shift, macroscopic=macro,
            clock=elapsed, weight=sheet.weight / len(group))


def real_rank(columns):
    echelon = {}
    for col in columns:
        col = tuple(map(Q, col))
        for index in range(len(col)):
            if not col[index]:
                continue
            if index not in echelon:
                echelon[index] = times(1 / col[index], col)
                break
            col = minus(col, times(col[index], echelon[index]))
    return len(echelon)
