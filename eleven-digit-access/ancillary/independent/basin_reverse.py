"""Reconstruct the inverse closure of the four-cycle at five bases."""

from collections import Counter
from pathlib import Path
from time import monotonic

from basin_core import F, Fint, gaps, preds


def cycle(m):
    return [
        (4*m + 1, 3*m, 3*m - 1, 2*m, m + 2),
        (4*m, 3*m, 2*m - 2, m, m - 1),
        (4*m + 1, 3*m, 3*m, 2*m + 1, m + 2),
        (4*m, 3*m, 2*m - 2, m + 1, m - 1),
    ]


def basin(m):
    base = 5*m
    frontier = cycle(m)
    for source, target in zip(frontier, frontier[1:] + frontier[:1]):
        assert F(source, base) == Fint(source, base) == target
    depth = dict.fromkeys(frontier, 0)
    layer = 0
    while frontier:
        layer += 1
        next_frontier = []
        for target in frontier:
            for source in preds(target, base):
                if source not in depth:
                    depth[source] = layer
                    next_frontier.append(source)
        frontier = next_frontier
    return depth


def frozen_depths():
    states = {}
    for line in Path(__file__).with_name("basin_m20.txt").read_text().splitlines():
        if line.startswith("S "):
            depth, *source = map(int, line.split()[1:])
            assert len(source) == 5 and tuple(source) not in states
            states[tuple(source)] = depth
    return states


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    for m in (20, 21, 50, 1000, 20001):
        start = monotonic()
        depths = basin(m)
        layers = Counter(depths.values())
        layer_sizes = [layers[i] for i in range(max(layers) + 2)]
        minimum_gaps = Counter(min(gaps(source, 5*m)) for source in depths)
        weight = sum((gap + 1)*count for gap, count in minimum_gaps.items())
        zero_inner = sum(source[-1] == 0 for source in depths)
        assert len(depths) == 235
        assert layer_sizes == [4, 27, 31, 38, 57, 36, 25, 17, 0]
        assert minimum_gaps == {0: 141, 1: 48, 2: 30, 3: 16}
        assert weight == 391 and zero_inner == 45
        if m == 20:
            assert depths == frozen_depths()
            print("m=20 matches the exhaustive list state for state and depth for depth.")
        print(f"m={m} B={5*m}: {len(depths)} states, layers {layer_sizes}, "
              f"minimum gap counts {dict(sorted(minimum_gaps.items()))}, "
              f"weighted sum {weight}, {zero_inner} zero-innermost states "
              f"({monotonic() - start:.2f} s)", flush=True)
