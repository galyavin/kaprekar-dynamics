"""Integer centered dynamics, derived from the borrowed digits."""
from collections import Counter

CENTER = (26293, 27557, 28481, 40447, 29865, 27357)


def centered(gaps):
    base = sum(gaps)
    return tuple(6 * sum(gaps[i:]) - (6 - i) * base for i in range(1, 6))


def pairs(x):
    a, b, c, d, e = x
    return ((a, -e - 6), (b, -d - 6), (c, -c - 6),
            (d, -b - 6), (e - 6, -a))


def step(x, word=None):
    levels = pairs(x)
    if word is None:
        word = sum((left >= right) << i
                   for i, (left, right) in enumerate(levels))
    ranked = [-6]
    for i, pair in enumerate(levels):
        ranked.extend(pair if word >> i & 1 else reversed(pair))
    return tuple(ranked[i] - ranked[10 - i] for i in range(5)), word


def trace(x):
    states, words, seen = [], [], {}
    while x not in seen:
        seen[x] = len(states)
        states.append(x)
        x, word = step(x)
        words.append(word)
    return states, words, seen[x]


def pieces(states, words, stop, minimum=3, max_word=60):
    t = 0
    while t < stop:
        for length in range(1, min(max_word, (stop - t - 1) // minimum) + 1):
            word = words[t:t + length]
            shift = tuple(y - x for x, y in zip(states[t], states[t + length]))
            end = t + length
            while (end + length < stop and words[end:end + length] == word
                   and all(y - x == v for x, y, v in zip(
                       states[end], states[end + length], shift))):
                end += length
            if end >= t + minimum * length:
                yield t, end, tuple(word), shift, states[t]
                t = end
                break
        else:
            yield t, t + 1, (words[t],), (), states[t]
            t += 1


if __name__ == '__main__':
    for residue in range(2):
        g = tuple(v + (residue if i == 4 else 0)
                  for i, v in enumerate(CENTER))
        states, words, end = trace(centered(g))
        print(g, 'start', states[0], 'end', end, 'period', len(states) - end)
        for block in pieces(states, words, end + 1):
            if block[3]:
                print(block)
        print('counts', Counter(words))
