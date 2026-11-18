"""Exact guarded affine words and maximal translation excursions.

Each row is an integer linear form in (u,v,B,1). Guards mean >= 0.
Sorting ties use the digit's original index, so all cells are disjoint.
No floating point, approximate branch tests, or asymptotic exits are used.
"""
from dataclasses import dataclass
from itertools import permutations

from five_digit import valid

U, V, B, ONE = ((1, 0, 0, 0), (0, 1, 0, 0),
                (0, 0, 1, 0), (0, 0, 0, 1))
IDENTITY = (U, V, B, ONE)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def value(row, base, pair):
    return sum(a * b for a, b in zip(row, (*pair, base, 1)))


def pull(row, matrix):
    return tuple(sum(row[k] * matrix[k][j] for k in range(4)) for j in range(4))


def compose(after, before):
    return tuple(pull(row, before) for row in after)


DOMAIN = (sub(B, (0, 0, 0, 2)), sub(U, ONE), sub(sub(B, U), ONE), V, sub(U, V))
DIGITS = (U, sub(V, ONE), sub(sub(B, V), ONE), sub(B, U))


@dataclass(frozen=True)
class Word:
    labels: tuple
    matrix: tuple
    guards: tuple
    cell_translation: tuple = None

    def contains(self, base, pair):
        return all(value(g, base, pair) >= 0 for g in self.guards)

    def destination(self, base, pair):
        return tuple(value(row, base, pair) for row in self.matrix[:2])

    def then(self, following):
        return Word(self.labels + following.labels,
                    compose(following.matrix, self.matrix),
                    tuple(dict.fromkeys(self.guards + tuple(
                        pull(g, self.matrix) for g in following.guards))))

    @property
    def translation(self):
        """Return (du,dv) only for a base-independent pure translation."""
        if self.cell_translation is not None:
            return self.cell_translation
        if all(self.matrix[i][:3] == IDENTITY[i][:3] for i in range(2)):
            return self.matrix[0][3], self.matrix[1][3]
        return None

    def restrict_displacement(self, delta):
        """Refine by the exact missing datum F(q)-q=delta.

        The two added equalities are guards, so the exit calculation checks
        them too. This handles integer-width cells and equality-line drift.
        """
        guards = list(self.guards)
        for i in range(2):
            equation = sub(sub(self.matrix[i], IDENTITY[i]), (0, 0, 0, delta[i]))
            guards.extend((equation, neg(equation)))
        return Word(self.labels, self.matrix, tuple(dict.fromkeys(guards)), tuple(delta))

    def repetitions(self, base, pair):
        """Maximal legal repetitions; None means a stationary periodic word.

        At repetition j each guard is g(q)+j*g_linear(delta). The first
        negative guard gives the exact exit. A nonzero translation in the
        bounded state space always has at least one decreasing guard.
        """
        delta = self.translation
        if delta is None:
            raise ValueError('Word is not a constant translation.')
        if not self.contains(base, pair):
            return 0
        if delta == (0, 0):
            return None
        limits = []
        for g in self.guards:
            slope = g[0] * delta[0] + g[1] * delta[1]
            if slope < 0:
                limits.append(value(g, base, pair) // -slope + 1)
        if not limits:
            raise AssertionError('Nonzero translation cannot stay in a bounded domain.')
        return min(limits)


def all_branches():
    branches = {}
    for order in permutations(range(4)):
        guards = list(DOMAIN) + [sub(V, ONE)]
        for left, right in zip(order, order[1:]):
            guard = sub(DIGITS[left], DIGITS[right])
            if left > right:  # lower index wins an equality
                guard = sub(guard, ONE)
            guards.append(guard)
        matrix = (sub(sub(B, ONE), DIGITS[order[3]]),
                  sub(DIGITS[order[0]], DIGITS[order[2]]), B, ONE)
        label = ''.join(map(str, order))
        branches[label] = Word((label,), matrix, tuple(guards))
    # v=0: the other three subtraction digits are all B-1.
    low, high = sub(U, ONE), sub(B, U)
    for label, first, second, guard in (
        ('Z0', low, high, sub(low, high)),
        ('Z1', high, low, sub(sub(high, low), ONE)),
    ):
        matrix = (sub(sub(B, ONE), second), sub(sub(B, ONE), first), B, ONE)
        branches[label] = Word((label,), matrix, DOMAIN + (neg(V), guard))
    return branches


BRANCHES = all_branches()


def branch(base, pair):
    if not valid(base, pair):
        raise ValueError('Invalid base or pair.')
    u, v = pair
    if not v:
        label = 'Z0' if u - 1 >= base - u else 'Z1'
    else:
        order = sorted(range(4), key=lambda i: (-value(DIGITS[i], base, pair), i))
        label = ''.join(map(str, order))
    return BRANCHES[label]


def itinerary(base, pair, length):
    word = Word((), IDENTITY, ())
    for _ in range(length):
        step = branch(base, pair)
        word = word.then(step)
        pair = step.destination(base, pair)
    return word


def excursion(base, pair, max_word=4):
    """Try lengths 1..max_word, returning an exact jump and elapsed time.

    A one-step fallback covers every state. This is NOT a proved bounded-block
    algorithm: the number of fallback steps/revisited words is uncontrolled.
    Nor may it skip over a target and claim first arrival without a capture test.
    """
    word = Word((), IDENTITY, ())
    q = pair
    for _ in range(max_word):
        step = branch(base, q)
        word = word.then(step)
        q = step.destination(base, q)
        if word.translation is None:
            refined = word.restrict_displacement((q[0] - pair[0], q[1] - pair[1]))
        else:
            refined = word
        if refined.translation is not None:
            repeats = refined.repetitions(base, pair)
            if repeats is None:
                return {'destination': pair, 'steps': 0, 'word': word.labels,
                        'repetitions': None, 'periodic': True}
            if repeats > 1:
                du, dv = refined.translation
                return {'destination': (pair[0] + repeats * du, pair[1] + repeats * dv),
                        'steps': len(word.labels) * repeats, 'word': word.labels,
                        'repetitions': repeats, 'periodic': False}
    step = branch(base, pair)
    return {'destination': step.destination(base, pair), 'steps': 1,
            'word': step.labels, 'repetitions': 1, 'periodic': False}


def first_capture(base, pair, jump, targets):
    """First target at offset 0..steps, including every skipped word prefix.

    targets maps pairs to cycle labels. An affine prefix of repetition j is
    prefix(q) + j * linear(prefix)(delta); solve two integer equations.
    """
    if pair in targets:
        return 0, targets[pair]
    if jump['periodic']:
        raise AssertionError('Unknown periodic state: classical cycle bridge failed.')
    repeats = jump['repetitions']
    length = len(jump['word'])
    delta = tuple((jump['destination'][i] - pair[i]) // repeats for i in range(2))
    prefix, best = IDENTITY, None
    for offset, label in enumerate(jump['word']):
        start = tuple(value(row, base, pair) for row in prefix[:2])
        slope = tuple(row[0] * delta[0] + row[1] * delta[1] for row in prefix[:2])
        for target, cycle in targets.items():
            candidate = None
            possible = True
            for i in range(2):
                difference = target[i] - start[i]
                if slope[i] == 0:
                    possible &= difference == 0
                elif difference % slope[i]:
                    possible = False
                else:
                    index = difference // slope[i]
                    if candidate is not None and candidate != index:
                        possible = False
                    candidate = index
            if candidate is None:
                candidate = 0
            if possible and 0 <= candidate < repeats:
                time = candidate * length + offset
                if best is None or time < best[0]:
                    best = time, cycle
        prefix = compose(BRANCHES[label].matrix, prefix)
    if jump['destination'] in targets:
        capture = jump['steps'], targets[jump['destination']]
        if best is None or capture[0] < best[0]:
            best = capture
    return best


def accelerated_arrival(base, pair, cache=None, targets=None, max_word=4):
    """Exact pair-clock arrival; safe even when a jump passes a cycle state.

    Termination follows from the finite functional graph and classical cycle
    completeness, not from an unproved uniform bound on the block count.
    """
    from five_digit import classical_cycles
    if cache is None:
        cache = {}
    if targets is None:
        targets = {p: label for label, c in classical_cycles(base).items() for p in c}
    q, path, seen = pair, [], set()
    while q not in cache:
        if q in targets:
            cache[q] = 0, targets[q]
            break
        if q in seen:
            raise AssertionError('Unknown terminal cycle.')
        seen.add(q)
        jump = excursion(base, q, max_word)
        capture = first_capture(base, q, jump, targets)
        if capture is not None:
            cache[q] = capture
            break
        path.append((q, jump['steps']))
        q = jump['destination']
    time, label = cache[q]
    for q, elapsed in reversed(path):
        time += elapsed
        cache[q] = time, label
    return cache[pair]
