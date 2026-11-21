"""Independent exhaustive zero-loss homogeneous circuit and drift checks."""

from collections import Counter, defaultdict
from functools import cache
from itertools import product
from math import gcd


def negate(row):
    return tuple(-x for x in row)


def primitive(row):
    divisor = gcd(*row)
    return tuple(x // divisor for x in row) if divisor else row


def canonical(rows):
    return tuple(sorted({primitive(row) for row in rows}))


@cache
def eliminate(rows):
    """Decide exact rational feasibility of rows (constant last), each >= 0."""
    if any(row[-1] < 0 and not any(row[:-1]) for row in rows):
        return False
    rows = tuple(row for row in rows if any(row[:-1]))
    if not rows:
        return True
    dim = len(rows[0]) - 1
    index = min(range(dim), key=lambda i: (
        sum(row[i] > 0 for row in rows) * sum(row[i] < 0 for row in rows), i
    ))
    pos = [row for row in rows if row[index] > 0]
    neg = [row for row in rows if row[index] < 0]
    joined = [row for row in rows if row[index] == 0]
    for a in pos:
        for b in neg:
            joined.append(tuple(a[index] * y - b[index] * x for x, y in zip(a, b)))
    return eliminate(canonical(row[:index] + row[index + 1:] for row in joined))


@cache
def feasible(dim, weak, strict=()):
    # Homogeneity lets us rescale every finite collection of strict margins to >= 1.
    if any(not any(row) for row in strict):
        return False
    rows = [row + (0,) for row in weak]
    rows.extend(row + (-1,) for row in strict)
    rows.extend(tuple(int(i == j) for i in range(dim)) + (-1,) for j in range(dim))
    return eliminate(canonical(rows))


def rank_form(phase, weights):
    """A sum-zero linear form on slots, expressed in successive level widths."""
    assert sum(weights) == 0
    return tuple(sum(w for level, w in zip(phase, weights) if level > i)
                 for i in range(max(phase)))


@cache
def reflections(phase):
    return tuple(rank_form(phase, (-1,) * i + (0,) * (7 - 2 * i) + (1,) * i)
                 for i in range(1, 4))


def regular_step(s):
    return tuple(s[i] for i in (0, 6, 1, 5, 2, 4, 3))


def edges(s):
    """Every zero-loss successor preserving the anchor, with its sign constraints."""
    if s[6] < s[0]:
        return ()
    d1, d2, d3 = reflections(s)
    choices = []
    for keep_left, keep_right in product((True, False), repeat=2):
        weak = {d1}
        out = [s[0], s[6]]
        for i, keep, left, right in ((1, keep_left, d1, d2), (2, keep_right, d2, d3)):
            weak.update((left, right) if keep else (negate(left), negate(right)))
            pair = (s[i], s[6 - i])
            out.extend(pair if keep else reversed(pair))
        out.append(s[3])
        weak = tuple(sorted(weak))
        if feasible(max(s), weak):
            choices.append((tuple(out), weak))
    return tuple(sorted(set(choices)))


def groups():
    words = defaultdict(list)
    for count in range(1, 8):
        for word in product(range(count), repeat=7):
            multiplicities = tuple(word.count(i) for i in range(count))
            if all(multiplicities):
                words[multiplicities].append(word)
    return words


def circuits(words):
    """Enumerate simple circuits and retain every jointly feasible guard set."""
    graph = {word: edges(word) for word in words}
    active = set(words)
    while True:
        outgoing = {source for source in active if any(dest in active for dest, _ in graph[source])}
        incoming = {dest for source in active for dest, _ in graph[source] if dest in active}
        smaller = outgoing & incoming
        if smaller == active:
            break
        active = smaller
    result = defaultdict(set)

    def visit(path, weak):
        for dest, new in graph[path[-1]]:
            if dest not in active or dest < path[0]:
                continue
            guards = tuple(sorted(set(weak).union(new)))
            if not feasible(max(dest), guards):
                continue
            if dest == path[0]:
                result[path].add(guards)
            elif dest not in path:
                visit(path + (dest,), guards)

    for start in sorted(active):
        visit((start,), ())
    return len(active), result


def sign_patterns(dim, rows, alternatives):
    """Split a circuit's feasible guard sets into strict signs and ties."""
    def visit(signs, weak, strict):
        if len(signs) == len(rows):
            yield signs, weak, strict
            return
        row = rows[len(signs)]
        other = negate(row)
        for sign in (-1, 0, 1):
            w = weak if sign else tuple(sorted(set(weak).union((row, other))))
            s = tuple(sorted(set(strict).union((row if sign == 1 else other,)))) if sign else strict
            if feasible(dim, w, s):
                yield from visit(signs + (sign,), w, s)

    for weak in alternatives:
        yield from visit((), weak, ())


def classify(path, alternatives):
    dim = max(path[0])
    rows = tuple(row for phase in path for row in reflections(phase))
    if all(regular_step(phase) == path[(i + 1) % len(path)] for i, phase in enumerate(path)):
        assert all(not feasible(dim, weak, (negate(row),))
                   for weak in alternatives for row in set(rows))
        return "regular", 0

    kinds = set()
    all_signs = set()
    for signs, weak, strict in sign_patterns(dim, rows, alternatives):
        all_signs.add(signs)
        ds = [signs[i:i + 3] for i in range(0, len(signs), 3)]
        # The top-rank identity needs separation both before and after each update.
        def top_slots(phase, group):
            return min(phase[i] for i in group) > max(phase[i] for i in set(range(7)) - group)

        if all(sign[1] < 0 for sign in ds) and all(top_slots(s, {0, 1, 3, 6}) for s in path):
            kind = "top four"
        elif all(sign[1] < 0 for sign in ds) and all(top_slots(s, {0, 1, 3, 4, 6}) for s in path):
            kind = "top five"
        elif all(sign[2] < 0 for sign in ds):
            kind = "W3"
        else:
            # Ask directly whether some tied-end kick choices violate the energy bound.
            works = tuple(rank_form(path[(i + 1) % len(path)], (1, -2, 1, 0, 0, 0, 0))
                          for i in range(len(path)))
            signs_of_kicks = [(-1, 1) if sign[0] == 0 else (sign[0],) for sign in ds]
            for kicks in product(*signs_of_kicks):
                bad = tuple(1 + 2 * sum(e * row[j] for e, row in zip(kicks, works)) for j in range(dim))
                if feasible(dim, weak, tuple(sorted(set(strict).union((bad,))))):
                    kind = "exception"
                    break
            else:
                kind = "energy"
        kinds.add(kind)
    assert len(kinds) == 1, (path, kinds)
    return kinds.pop(), len(all_signs)


def check_177(path, kind, alternatives):
    """Task 093's separate strict-sign classification for leading periods 1 and 3."""
    if not max(path[0]):
        return "origin"
    if all(s[6] > s[0] for s in path):
        return "strict anchor"
    if len(path) == 1:
        assert kind in ("top four", "top five", "W3", "energy")
        phase = path[0]
        d1, d2, d3 = reflections(phase)
        assert not any(d1)
        assert all(not feasible(max(phase), tuple(sorted(set(w).union((d2, d3))))) for w in alternatives)
        assert max(phase[:4]) > 0
        return "tied fixed"
    highest = max(path[0])
    assert all(s.count(highest) == 1 for s in path)
    loc = tuple(s.index(highest) for s in path)
    assert any(loc[i:] + loc[:i] == (6, 1, 3) for i in range(3))
    for phase in path:
        d1, d2, _ = reflections(phase)
        i = phase.index(highest)
        if i == 6:
            assert phase[6] > phase[0]
        else:
            assert max(phase[:4]) > 0
        if i == 1:
            assert all(not feasible(highest, tuple(sorted(set(w).union((d2,))))) for w in alternatives)
    return "tied triple"


def main():
    initial = active = patterns = 0
    totals = Counter()
    by_kind = Counter()
    task_093 = Counter()
    singular = []
    for multiplicities, words in sorted(groups().items()):
        initial += len(words)
        retained, feasible_circuits = circuits(words)
        active += retained
        for path, alternatives in sorted(feasible_circuits.items()):
            totals[len(path)] += 1
            kind, number = classify(path, alternatives) if max(path[0]) else ("origin", 0)
            by_kind[kind] += 1
            patterns += number
            if kind not in ("regular", "origin"):
                singular.append((multiplicities, path, kind, number))
            if len(path) in (1, 3):
                task_093[check_177(path, kind, alternatives)] += 1

    print(f"rank words {initial}; locally feasible graph survivors {active}")
    print("circuits by period:", sorted(totals.items()))
    print("circuits by certificate:", sorted(by_kind.items()), "; singular sign patterns", patterns)
    print("177-dividing leading circuits:", sorted(task_093.items()))
    for multiplicities, path, kind, number in singular:
        text = " ".join("".join(map(str, phase)) for phase in path)
        print(f"{multiplicities} {kind} {number}: {text}")
        if kind == "exception":
            assert multiplicities == (2, 5) and set(path) == {(1, 1, 0, 1, 0, 1, 1), (1, 1, 1, 1, 0, 0, 1)}
    assert totals == {1: 7, 2: 8, 3: 6, 4: 1, 5: 3, 6: 789, 7: 1, 8: 2}
    assert task_093 == {"origin": 1, "strict anchor": 5, "tied fixed": 5, "tied triple": 2}
    assert by_kind == {"origin": 1, "regular": 795, "top four": 8, "top five": 1,
                       "W3": 2, "energy": 9, "exception": 1}
    assert patterns == 27


if __name__ == "__main__":
    main()
