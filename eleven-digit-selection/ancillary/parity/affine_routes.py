"""Guarded affine routes computed from the centered borrowed-pair map."""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cache
from math import lcm
from primitive import CENTER, pairs, step, centered

WIDTH = Q(1, 10**6)
INITIAL_X = centered(CENTER)
HORIZON = 65


@dataclass(frozen=True)
class Form:
    """Affine function of the initial centered coordinates X."""
    row: tuple
    constant: Q = Q(0)

    def __neg__(self):
        return Form(tuple(-a for a in self.row), -self.constant)

    def __add__(self, other):
        if not isinstance(other, Form):
            return Form(self.row, self.constant + other)
        return Form(tuple(a + b for a, b in zip(self.row, other.row)),
                    self.constant + other.constant)

    __radd__ = __add__

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, scalar):
        return Form(tuple(a * scalar for a in self.row), self.constant * scalar)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * Q(1, scalar)

    def evaluate(self, x):
        return sum(a * b for a, b in zip(self.row, x)) + self.constant

    @property
    def leading(self):
        return Form(self.row)

    def __str__(self):
        return '(' + ','.join(map(str, (*self.row, self.constant))) + ')'


ZERO = Form((0,) * 5)
VARIABLES = tuple(Form(tuple(int(i == j) for i in range(5))) for j in range(5))
PARTNER = (1, 2, 3, 4, -1)


@cache
def normalized_interval(row):
    center = Q(sum(a * x for a, x in zip(row, INITIAL_X)), 180000)
    # X_i = iB - 6 sum_{j < i} G_j.
    radius = 6 * WIDTH * sum(abs(sum(row[j:])) for j in range(5))
    return center - radius, center + radius


def sign(form):
    """Eventual weak sign on the full closed normalized rectangle."""
    if not any(form.row):
        return form.constant >= 0
    lo, hi = normalized_interval(form.row)
    if lo > 0:
        return True
    if hi < 0:
        return False
    raise ArithmeticError(f'nonseparated guard {form} in [{lo}, {hi}]')


def strict_lower(form):
    lo, _ = normalized_interval(form.row)
    return lo


def compare(x):
    return tuple(a - b for a, b in pairs(x))


def choose(x):
    return sum(sign(f) << i for i, f in enumerate(compare(x)))


def guards(x, word):
    # All five differences are divisible by 6 on every physical coset.
    return tuple(g if word >> i & 1 else -g - 6
                 for i, g in enumerate(compare(x)))


def translate(word, delta):
    current = delta
    slopes = []
    for branch in word:
        # Difference of the primitive affine update on translated inputs.
        slopes.extend(s - t for s, t in zip(
            guards(current, branch), guards((0,) * 5, branch)))
        current = tuple(s - t for s, t in zip(
            step(current, branch)[0], step((0,) * 5, branch)[0]))
    return tuple(slopes), current


def macros(state):
    """Return the first substantial translation or one guarded primitive step."""
    cur, pulled, word = state, [], []
    for _ in range(HORIZON):
        choice = choose(cur)
        pulled.extend(guards(cur, choice))
        word.append(choice)
        cur = step(cur, choice)[0]
        delta = tuple(b - a for a, b in zip(state, cur))
        if any(any(v.row) for v in delta):
            continue
        v = tuple(int(f.constant) for f in delta)
        if not any(v):
            raise ArithmeticError('unexpected cycle with nonzero leading state')
        slopes, translated = translate(word, v)
        if translated != v:
            continue
        candidates = [(-slope, guard) for guard, slope in zip(pulled, slopes)
                      if slope < 0]
        if not candidates:
            raise ArithmeticError('unexpected unbounded translation')
        # Unambiguous minimum throughout the rectangle for sufficiently large B.
        k, form = min(candidates, key=lambda pair: (
            sum(a * b for a, b in zip(pair[1].row, INITIAL_X)) / pair[0],
            pair[1].constant / pair[0]))
        for divisor, guard in candidates:
            if not sign(guard / divisor - form / k):
                raise ArithmeticError('selected return cannot give first exit')
        if not any(form.row):
            continue
        return tuple(word), v, k, form
    return (choose(state),), (), 0, ZERO


def repeat(state, elapsed, block, remainder):
    word, v, divisor, numerator = block
    count = (numerator - remainder) / divisor + 1
    return tuple(x + count * shift for x, shift in zip(state, v)), elapsed + count * len(word)


def witness(odd):
    state = tuple(v + odd * d for v, d in zip(VARIABLES, PARTNER))
    root = INITIAL_X
    elapsed = ZERO
    congruences = []
    print('ODD' if odd else 'EVEN')
    while any(any(s.row) for s in state):
        block = macros(state)
        word, delta, divisor, numerator = block
        if not divisor:
            print('step', word[0], 'at', elapsed.evaluate(root))
            state, elapsed = step(state, word[0])[0], elapsed + 1
        else:
            remainder = numerator.evaluate(root) % divisor
            count = (numerator.evaluate(root) - remainder) / divisor + 1
            print('repeat', word, 'at', elapsed.evaluate(root), 'delta', delta,
                  'numerator', numerator, 'divisor', divisor, 'remainder', remainder, 'count', count)
            state, elapsed = repeat(state, elapsed, block, remainder)
            congruences.append((numerator, divisor, remainder))
            print('arrival', tuple(str(s) for s in state))
    print('constant endpoint', tuple(s.constant for s in state),
          'elapsed', elapsed, 'at center', elapsed.evaluate(root), 'clock', strict_lower(elapsed))
    return congruences, elapsed


if __name__ == '__main__':
    for odd in range(2):
        witness(odd)
