from investigate import transition


def stable(base):
    return [(u, v) for u in range(3, base, 2) for v in range(1, u, 2)]


def cycles(base):
    unseen = set(stable(base))
    result = []
    while unseen:
        start = min(unseen)
        pair, path = start, []
        while pair not in path:
            path.append(pair)
            unseen.remove(pair)
            pair = transition(base, pair)
        assert pair == start
        result.append(path)
    return result
