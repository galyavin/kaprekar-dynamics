"""Replay the retirement structure and the one-active obstruction exactly."""

import argparse
import json
from itertools import accumulate
from math import prod
from pathlib import Path


def step(state):
    endpoint = state[-1]
    sign = (endpoint > 0) - (endpoint < 0)
    loss = max(0, -1 - endpoint)
    return (
        endpoint - 3 * sign + 2 * loss,
        state[0] + 2 * loss,
        *(value - sign + loss for value in state[1:-1]),
    ), loss


def verify_retirement_structure(certificate, description="archived"):
    order = certificate["order"]
    intervals = certificate["insertions"]
    assert len(intervals) == order
    assert all(0 <= start < end for start, end in intervals)
    assert all(left[1] <= right[0] for left, right in zip(intervals, intervals[1:]))
    assert intervals[-1][1] <= certificate["cycle"][0]
    checkpoints = {time for interval in intervals for time in interval}
    state = tuple(certificate["source"])
    records = {}
    first_nonpositive = {}
    for time in range(certificate["cycle"][1]):
        if time in checkpoints:
            records[time] = state
        if state[-1] <= 0:
            label = (order - 1 - time) % order
            first_nonpositive.setdefault(label, time)
        state, _ = step(state)

    retirement_order = sorted(first_nonpositive, key=first_nonpositive.get)
    assert len(retirement_order) == order
    diagonal = []
    for index, (start, end) in enumerate(intervals):
        assert (end - start) % order == 0
        displacement = tuple(
            records[start][(label + start) % order] - records[end][(label + end) % order]
            for label in retirement_order
        )
        assert all(first_nonpositive[label] <= start for label in retirement_order[:index])
        assert all(first_nonpositive[label] >= end for label in retirement_order[index:])
        assert all(value == 0 for value in displacement[:index])
        assert displacement[index] > 0
        diagonal.append(displacement[index])
    last_start, last_end = intervals[-1]
    determinant = prod(diagonal)
    print(f"{description} order={order} triangular determinant={determinant} "
          f"last projected word={last_end - last_start} last drift={diagonal[-1]}")
    return [first_nonpositive[label] for label in retirement_order], diagonal


def last_active_turn(state):
    assert state[-1] == 0
    # Give the active endpoint the formal offset zero and force its positive branch.
    moved = (-3, state[0], *(value - 1 for value in state[1:-1]))
    reads = []
    for _ in range(len(state) - 1):
        reads.append(moved[-1])
        moved, _ = step(moved)
    return (*moved[:-1], 0), moved[-1], tuple(reads)


def ordinary_cycle(state):
    loss = 0
    seen = {}
    for time in range(500):
        if state in seen:
            entry, loss_at_entry = seen[state]
            return entry, time - entry, loss - loss_at_entry
        seen[state] = (time, loss)
        state, step_loss = step(state)
        loss += step_loss
    raise AssertionError("Obstruction replay did not reach its stated cycle")


def verify_obstruction():
    projected = (
        (-1, 1, 4, 4, 0),
        (2, 3, 2, 3, 0),
        (0, -2, 0, -2, 0),
        (0, 2, -2, 4, 0),
    )
    turns = [last_active_turn(state) for state in projected]
    assert tuple(row[0] for row in turns) == projected[1:] + projected[:1]
    offsets = tuple(accumulate((row[1] for row in turns), initial=0))
    assert offsets == (0, -1, -5, -5, -6)
    drift = -offsets[-1]
    print(f"order=5 obstructed projection period={len(projected) * 5} drift={drift}")

    exits = {}
    for phase, state in enumerate(projected):
        previous = [offset + (drift if j >= phase else 0)
                    for j, offset in enumerate(offsets[:-1])]
        guard = min(previous) - offsets[phase]
        for endpoint in range(1 - guard, 1):
            exits[phase, endpoint] = ordinary_cycle((*state[:-1], endpoint))
    assert exits == {
        (0, 0): (39, 100, 36),
        (1, 0): (9, 100, 36),
        (2, -3): (0, 37, 12),
        (2, -2): (104, 100, 36),
        (2, -1): (237, 100, 36),
        (2, 0): (11, 100, 36),
    }
    for (phase, endpoint), (entry, period, loss) in exits.items():
        print(f"  exit phase={phase} endpoint={endpoint} entry={entry} "
              f"ordinary period={period} cycle loss={loss}")
    full_rank_source = {
        "order": 5,
        "source": (111, 136, 48, 16, 1),
        "insertions": ((0, 5), (11, 31), (37, 92), (254, 264), (303, 323)),
        "cycle": (362, 462),
    }
    times, drifts = verify_retirement_structure(full_rank_source, "obstructed")
    assert times == [5, 31, 92, 264, 323]
    assert drifts == [6, 14, 28, 6, 6] and prod(drifts) == 84672
    assert ordinary_cycle(full_rank_source["source"]) == (362, 100, 36)
    state = full_rank_source["source"]
    for _ in range(323):
        state, _ = step(state)
    assert state == projected[0]


def verify_terminal(order):
    if order % 2:
        pairs = (order - 7) // 2
        start = (2, 2, 1, 0, *(3, 0) * pairs, 3, 2, -1)
        endpoints = (-1, 3, 3, *(-1, 3) * pairs, -1, 1, 1, -1,
                     1, -1, -1, *(3, -1) * pairs, 3, -1, -1, 5)
    else:
        pairs = (order - 10) // 2
        start = (2, 2, -1, 0, 3, 0, 0, *(3, 0) * pairs, 3, 2, -1)
        endpoints = (-1, 3, 3, *(-1, 3) * pairs, -1, 0, 3, -1, -1, 3, 1,
                     1, -1, -1, *(3, -1) * pairs, 2, 0, -1, 1, 1, -1, -1)
    state = start
    seen = set()
    for endpoint in endpoints:
        assert state not in seen
        seen.add(state)
        assert state[-1] == endpoint
        state, loss = step(state)
        assert loss == 0
    assert state == start and len(seen) == 2 * order
    print(f"order={order} constructed least terminal period={len(seen)} loss=0")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    paths = (args.archives / "lc_classification/root-certificates.jsonl",
             args.archives / "lc2_order12/certificate.jsonl")
    certificates = [json.loads(line) for path in paths for line in path.read_text().splitlines()]
    assert sorted(row["order"] for row in certificates) == list(range(3, 14))
    for certificate in sorted(certificates, key=lambda row: row["order"]):
        verify_retirement_structure(certificate)
    verify_obstruction()
    for order in (7, 9, 10, 12, 13):
        verify_terminal(order)


if __name__ == "__main__":
    main()
