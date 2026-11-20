"""Reconstruct ordinary roots from a staged projected itinerary."""
import json
import math
import sys
from pathlib import Path

R = 12
LAST = 1 << (R - 1)
ALL = (1 << R) - 1
INF = 10**30


def projected(state, mask):
    return tuple(0 if mask >> i & 1 else v for i, v in enumerate(state))


def forced_step(state, mask):
    z = state[-1]
    if mask & LAST:
        e, h = 1, 0
    else:
        e, h = (z > 0) - (z < 0), max(0, -1-z)
    return ((z - 3*e + 2*h, state[0] + 2*h, *(v - e + h for v in state[1:-1])),
            ((mask << 1) & ALL) | (mask >> (R - 1)), h)


def read_stage(line):
    left, right = line.split("|")
    mask, row = left.split(":")
    return int(mask), tuple(map(int, row.split())), *map(int, right.split())


def itinerary(filename):
    lines = Path(filename).read_text().splitlines()
    return [read_stage(line) for line in lines if "|" in line]


def build(stages):
    source = [None] * R
    state, mask = (0,) * R, ALL
    reads = [INF] * R
    clock = 0
    intervals = []
    derived = []
    for expected_mask, expected_state, s, p, q, value in stages:
        assert mask == expected_mask and projected(state, mask) == expected_state
        period = math.lcm(p, R)
        assert s <= q < s + period
        for j in range(q):
            if mask & LAST:
                label = (R - 1 - clock - j) % R
                reads[label] = min(reads[label], state[-1])
            state, mask, _ = forced_step(state, mask)
        at_q, mask_q = state, mask
        assert mask_q & LAST
        block_reads = [INF] * R
        for j in range(period):
            if mask & LAST:
                label = (R - 1 - clock - q - j) % R
                block_reads[label] = min(block_reads[label], state[-1])
            state, mask, _ = forced_step(state, mask)
        assert mask == mask_q and projected(at_q, mask) == projected(state, mask)
        delta = tuple(b-a for a,b in zip(at_q,state))
        retiring = (R - 1 - clock - q) % R
        assert delta[-1] < 0
        assert value + block_reads[retiring] - state[-1] > 0
        step_size = -delta[-1]
        threshold = 1-reads[retiring]-value+at_q[-1]
        repetitions = max(1, (threshold+step_size-1)//step_size)
        initial = value - at_q[-1] - repetitions*delta[-1]
        assert initial + reads[retiring] >= 1
        assert source[retiring] is None
        source[retiring] = initial
        for i, d in enumerate(delta):
            label = (i - clock - q) % R
            if mask >> i & 1:
                reads[label] = min(reads[label], block_reads[label] + min(0,(repetitions-1)*d))
        state = tuple(v+repetitions*d for v,d in zip(at_q,delta))
        assert source[retiring] + state[-1] == value
        state = (*state[:-1], value)
        mask ^= LAST
        intervals.append([clock + q, clock + q + period])
        clock += q + repetitions*period
        derived.append({"label":retiring,"period":period,"repeats":repetitions,"value":value,
                        "phase":clock,"drift":delta[-1]})
    assert mask == 0 and all(x is not None for x in source)
    print("retirement state", state, "at", clock, file=sys.stderr)
    seen = {}
    loss = time = 0
    while state not in seen:
        seen[state] = (time, loss)
        state, _, h = forced_step(state, 0)
        time += 1
        loss += h
        assert time <= 100000
    first, first_loss = seen[state]
    assert first_loss == loss
    print("stage data", json.dumps(derived, separators=(",", ":")), file=sys.stderr)
    return {"order":R,"source":source,"insertions":intervals,"cycle":[clock+first,clock+time]}


if __name__ == "__main__":
    print(json.dumps(build(itinerary(sys.argv[1])), separators=(",", ":")))
