"""Rotate every fifth long return to its central-sheet phase for comparison."""
from collections import Counter, defaultdict
from fractions import Fraction as Q

import exclude
from affine_routes import VARIABLES, ZERO
from primitive import step

A = (Q(1, 6), Q(55, 6), -5, Q(-55, 6), Q(95, 6))
H = tuple(a / Q(2) for a in A)
SHEET = (A, H, ZERO.row, A, H)
PREFIX = (Q(5, 36), Q(-55, 36), Q(5, 6), Q(55, 36), Q(-95, 36))
PREDICTED_SPEEDS = {(Q(18, 5), Q(12, 5)),
                    (Q(33, 10), Q(21, 10)),
                    (Q(16, 5), Q(2))}

original = exclude.explore
counts = Counter()
fibers = defaultdict(set)
rotations = defaultdict(set)
phase0_words = defaultdict(set)
prefix_constants = []


def explore(state, elapsed, x, basis, depth=0, history=()):
    if depth != 4:
        return original(state, elapsed, x, basis, depth, history)
    for _ in range(exclude.MAX_SHORT):
        word, delta, modulus, _ = exclude.macro(state)
        if modulus:
            break
        state, elapsed = step(state, word[0])[0], elapsed + 1
    else:
        raise ArithmeticError('missing sheet return')
    assert elapsed.row == PREFIX
    prefix_constants.append(elapsed.constant)

    found = set()
    for j, choice in enumerate(word):
        if tuple(s.row for s in state) == SHEET:
            a, b, c, d, e = state
            f = exclude.plain((c, d - a, e - b))
            assert all(z % 6 == 0 for z in f)
            v = tuple(Q(-z, len(word)) for z in delta[:2])
            assert delta[2] == 0 and delta[3:] == delta[:2]
            assert v in PREDICTED_SPEEDS
            found.add(v)
            fword = word[j:] + word[:j]
            rotations[v].add((j, len(word)))
            phase0_words[v].add(fword)
            fibers[v].add(f)
        state = step(state, choice)[0]
        delta = tuple(s - t for s, t in zip(
            step(delta, choice)[0], step((0,) * 5, choice)[0]))
    assert len(found) == 1, (found, word)
    counts[found.pop()] += 1


if __name__ == '__main__':
    exclude.explore = explore
    basis = tuple(tuple(6 * int(i == j) for j in range(5)) for i in range(5))
    explore(VARIABLES, ZERO, (0,) * 5, basis)
    assert set(counts) == PREDICTED_SPEEDS
    assert sum(counts.values()) == 6243
    print('prefix time form', PREFIX, 'offset range',
          (min(prefix_constants), max(prefix_constants)))
    for v in sorted(counts):
        print('sheet speed', v, 'entry cosets', counts[v],
              'phase rotations and lengths', sorted(rotations[v]))
        print('phase-zero fibers', sorted(fibers[v]))
        print('phase-zero words', sorted(phase0_words[v]))
