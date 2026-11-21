"""Exact finite checks for Task 226, using only the Python standard library."""

from collections import Counter
from itertools import product
from random import Random


def sign(value):
    return (value > 0) - (value < 0)


def step(root):
    read = root[-1]
    e = sign(read)
    loss = max(-1 - read, 0)
    successor = (
        read - 3 * e + 2 * loss,
        root[0] + 2 * loss,
        *(value - e + loss for value in root[1:-1]),
    )
    return successor, loss


def initial_from_reads(order, reads):
    assert len(reads) >= order and min(reads) >= -1
    initial = []
    for slot, first_read in enumerate(range(order - 1, -1, -1)):
        corrections = sum(
            -sign(reads[t]) for t in range(first_read) if slot + t > 0
        )
        initial.append(reads[first_read] - corrections)
    return tuple(initial)


def replay_cycle(order, word):
    prefix = tuple(word[t % len(word)] for t in range(order))
    root = initial_from_reads(order, prefix)
    states = []
    current = root
    for read in word:
        assert current[-1] == read, (order, word, current)
        states.append(current)
        current, loss = step(current)
        assert loss == 0
    assert current == root
    return states, len(set(states))


def negate(state):
    u, e = state
    return -u, -e


def graph():
    choices = (-1, 0, 1)
    states = tuple(product(range(-3, 4), choices))
    edges = {}
    for u, e in states:
        allowed_signs = {
            -1: tuple(f for f in choices if u <= f),
            0: (u,) if u in choices else (),
            1: tuple(f for f in choices if u >= f),
        }
        edges[u, e] = tuple(
            (u - 2 * e, f) for f in allowed_signs[e] if (u - 2 * e, f) in states
        )

    def can_return(source):
        seen = set()
        pending = list(edges[source])
        while pending:
            target = pending.pop()
            if target == source:
                return True
            if target not in seen:
                seen.add(target)
                pending.extend(edges[target])
        return False

    live = {s for s in states if can_return(s)}
    assert len(states) == 21 and len(live) == 13
    assert all({negate(t) for t in edges[s]} == set(edges[negate(s)]) for s in states)
    print(f"Graph has {len(states)} nominal states, {len(live)} on closed walks.")
    return {s: tuple(t for t in edges[s] if t in live) for s in sorted(live)}


def expected_count(order):
    if order % 2:
        return 2**order + 1
    return (2 ** (order // 2) + 1) ** 2


def walk_counts(edges, limit):
    distributions = {source: Counter({source: 1}) for source in edges}
    for order in range(1, limit + 1):
        for source, current in distributions.items():
            extended = Counter()
            for state, frequency in current.items():
                for target in edges[state]:
                    extended[target] += frequency
            distributions[source] = extended
        if order >= 3:
            count = sum(counts[negate(source)] for source, counts in distributions.items())
            assert count == expected_count(order), (order, count)
    print(f"Transfer counts agree with the formula at orders 3 through {limit}.")


def graph_words(edges, order):
    def extend(path, remaining, target):
        if remaining == 0:
            if path[-1] == target:
                yield path[:-1]
            return
        for state in edges[path[-1]]:
            yield from extend((*path, state), remaining - 1, target)

    for source in edges:
        for path in extend((source,), order, negate(source)):
            states = path + tuple(map(negate, path))
            following = states[1:] + states[:1]
            yield tuple(1 + u - f if e == 1 else e for (u, e), (_, f) in zip(states, following))


def ternary_words(order):
    for first in product((-1, 0, 1), repeat=order):
        signs = first + tuple(-e for e in first)
        period = len(signs)
        writes = tuple(
            3 * signs[t] + sum(signs[(t + j) % period] for j in range(2, order))
            for t in range(period)
        )
        reads = tuple(writes[t] - 1 if e == 1 else e for t, e in enumerate(signs))
        if all(
            sign(reads[t]) == signs[t] and reads[(t + order) % period] == reads[t] - writes[t]
            for t in range(period)
        ):
            yield reads


def verify_entry(states):
    target = states[0]
    cyclic_predecessor = states[-1]
    cyclic_read = cyclic_predecessor[-1]
    # The candidates come from the interiors of the two affine branches and the two boundary reads.
    found = {}
    for read in {target[0] + 3, 1 - target[0], -1, 0}:
        e = sign(read)
        loss = max(-1 - read, 0)
        predecessor = (
            target[1] - 2 * loss,
            *(value + e - loss for value in target[2:]),
            read,
        )
        successor, actual_loss = step(predecessor)
        if successor == target:
            assert actual_loss == 0
            found[read] = predecessor
    expected_reads = {0, 3} if target[0] == 0 else {-1, 5} if target[0] == 2 else {cyclic_read}
    assert set(found) == expected_reads
    assert found[cyclic_read] == cyclic_predecessor
    for read in expected_reads - {cyclic_read}:
        assert found[read] not in states
    return len(found) - 1


def verify_words(edges):
    words_checked = 0
    alternative_entries = 0
    for order in range(3, 14):
        words = set()
        phase_periods = Counter()
        for word in graph_words(edges, order):
            assert word not in words
            words.add(word)
            states, least_period = replay_cycle(order, word)
            phase_periods[least_period] += 1
            alternative_entries += verify_entry(states)
            if 5 in word and 0 not in word:
                assert order % 2
        assert len(words) == expected_count(order)
        if order <= 9:
            assert words == set(ternary_words(order))
        assert all(phases % period == 0 for period, phases in phase_periods.items())
        counts = " ".join(
            f"{period}:{phases // period}" for period, phases in sorted(phase_periods.items())
        )
        print(f"r={order:2d} marked words={len(words):4d}   cycles by least period {counts}")
        words_checked += len(words)
    print(f"Replayed {words_checked} marked words and {alternative_entries} alternative predecessors.")
    print("Exhaustive ternary sign enumeration agrees at orders 3 through 9.")


def lc3_core(order):
    if order % 2:
        pairs = (order - 7) // 2
        return (
            (-1, 3, 3) + (-1, 3) * pairs + (-1, 1, 1, -1),
            (1, -1, -1) + (3, -1) * pairs + (3, -1, -1, 5),
        )
    pairs = (order - 10) // 2
    return (
        (-1, 3, 3) + (-1, 3) * pairs + (-1, 0, 3, -1, -1, 3, 1),
        (1, -1, -1) + (3, -1) * pairs + (2, 0, -1, 1, 1, -1, -1),
    )


def verify_six_blocks():
    u_block = (-1, 3, 3, -1, 3, -1)
    v_block = (1, -1, -1, 3, -1, 1)
    checked = 0
    for core_order in (7, *range(9, 32)):
        u, v = lc3_core(core_order)
        for repeats in range(7):
            order = core_order + 6 * repeats
            word = u_block * repeats + u + v_block * repeats + v
            states, period = replay_cycle(order, word)
            assert len(word) == period == 2 * order
            assert word.count(5 if order % 2 else 2) == 1
            if (core_order, repeats) == (9, 1):
                assert states[0] == (2, 2, 1, 0, 3, 0, 3, 2, -1, 0, 3, 0, 3, 2, -1)
            checked += 1
    print(f"Six-block family: {checked} words at orders 7 and 9 through 67.")
    print("LC3(9) plus one six-block agrees with the registered order-15 lock at phase zero.")


def verify_other_periods():
    controls = (
        (4, (-1, -1, 0, 1, 1, 0)),
        (6, (-1, 3, -1, 3, 0)),
        (8, (-1, -1, -1, 5, -1, 5, 0, 1, 1, -1, 1, -1, 1, 0)),
        (9, (-1, 3, -1, 3, 0, -1, 3, 0)),
        (10, (-1, 3, -1, 3, -1, 3, -1, 3, 0)),
        (10, (-1, 3, 0)),
        (11, (-1, 3, -1, 3, -1, 3, -1, 3, 0, 0)),
        (11, (-1, -1, 3, 0, -1, 3, 0, -1, 3, 0, -1, 2, 2)),
    )
    for order, word in controls:
        states, least_period = replay_cycle(order, word)
        assert len(word) == least_period and (2 * order) % least_period
        print(f"Outside-period control r={order:2d} period={least_period:2d} root={states[0]}")


def verify_label_drift():
    rng = Random(226)
    checked = 0
    for order in range(3, 18):
        for _ in range(20):
            initial = tuple(rng.randint(-20, 50) for _ in range(order))
            current = initial
            signs = [0] * order
            losses = [0] * order
            for t in range(order * rng.randint(1, 10)):
                signs[t % order] += sign(current[-1])
                current, loss = step(current)
                losses[t % order] += loss
            kick = sum(signs) - sum(losses)
            for label, (old, new) in enumerate(zip(reversed(initial), reversed(current))):
                successor = (label + 1) % order
                assert new - old == -kick - 2 * signs[label] + losses[label] + signs[successor] + losses[successor]
                checked += 1
    print(f"Label drift identity agrees on {checked} exact label displacements.")


if __name__ == "__main__":
    graph_edges = graph()
    walk_counts(graph_edges, 80)
    verify_words(graph_edges)
    verify_six_blocks()
    verify_other_periods()
    verify_label_drift()
