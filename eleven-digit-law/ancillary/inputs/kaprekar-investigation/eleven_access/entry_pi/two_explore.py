from collections import Counter
from random import Random


def induced(x, y, m):
    if x > y:
        return (x - y + m, x, y), "A"
    if y > x + m:
        return (x + m, y - x, y - x - m), "B"
    return (x + m - y, y, x), "C"


rng = Random(98)
final_words = Counter()
block_words = Counter()
ratios = []
for _ in range(1000):
    start = p = tuple(rng.randint(1, 1000) for _ in range(3))
    word = ""
    for i in range(200):
        p, w = induced(*p)
        word += w
    ratios.append((max(p), max(start)))
    final_words[word[-12:]] += 1
    for j in range(40, len(word) - 6):
        block_words[word[j:j+6]] += 1
print("Final max larger:", sum(final > initial for final, initial in ratios))
print("Final max unchanged:", sum(final == initial for final, initial in ratios))
print("Most common suffixes:", final_words.most_common(15))
print("Most common blocks:", block_words.most_common(20))
