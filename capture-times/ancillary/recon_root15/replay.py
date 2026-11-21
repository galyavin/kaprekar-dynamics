#!/usr/bin/env python3
"""Reconstruct the four finite order-15 witnesses and the length screen."""

import argparse
from collections import deque
from hashlib import sha256
import json
from math import isqrt, prod
from pathlib import Path


ORDER = 15
INPUT_HASHES = {
    "order15-alternate-certificates.jsonl":
        "aa48f0d2b31757e23e457d6fe518b992113b57a554da52494c8a01050156c9ee",
    "order15-alternate-anatomy.jsonl":
        "094bb7389bd2095d641eca9ced8f2bee89fd947411a090b1a26aede40fc34a29",
}


def sign_and_loss(z):
    return (z > 0) - (z < 0), max(0, -1 - z)


def root_step(state, sign, loss):
    """P9 and P10 slot order, with slot 15 rotating into slot 1."""
    common = loss - sign
    return (
        state[-1] - 3 * sign + 2 * loss,
        state[0] + 2 * loss,
        *(value + common for value in state[1:-1]),
    )


def cycle_at(candidate):
    """Derive the least loss-free return of the archived phase candidate."""
    state = candidate
    states = []
    reads = []
    seen = set()
    while state not in seen:
        assert len(states) < 2 * ORDER, "No return within the claimed period"
        sign, loss = sign_and_loss(state[-1])
        assert loss == 0, "Terminal candidate has borrowing loss"
        states.append(state)
        reads.append(state[-1])
        seen.add(state)
        state = root_step(state, sign, loss)
    assert state == candidate
    return states, reads


def in_read_labels(state, time):
    # At time zero the rightmost slot reads first. Each step rotates once.
    return tuple(state[(time - 1 - label) % ORDER] for label in range(ORDER))


def anchored_span(state):
    return max(0, max(state)) - min(0, min(state))


def exact_determinant(rows):
    """Fraction-free Bareiss elimination independent of the diagonal check."""
    rows = [row.copy() for row in rows]
    sign = previous_pivot = 1
    for col in range(len(rows) - 1):
        pivot = next((j for j in range(col, len(rows)) if rows[j][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            rows[col], rows[pivot] = rows[pivot], rows[col]
            sign = -sign
        current_pivot = rows[col][col]
        for j in range(col + 1, len(rows)):
            for k in range(col + 1, len(rows)):
                numerator = current_pivot * rows[j][k] - rows[j][col] * rows[col][k]
                rows[j][k], remainder = divmod(numerator, previous_pivot)
                assert remainder == 0
            rows[j][col] = 0
        previous_pivot = current_pivot
    return sign * rows[-1][-1]


def replay_reference(source, intervals, arrival, through, terminal):
    """Capture interval endpoints, literal projected words, and first reads."""
    checkpoint_times = {time for pair in intervals for time in pair}
    checkpoints = {}
    retirements = {}
    retire_order = []
    words = []
    active_word = None
    stage = 0
    cumulative_loss = 0
    reference_norm = reference_span = 0
    first_terminal_visit = None
    terminal_states = set(terminal)
    state = source

    # Cross-check the slot update against the registered order-15 read
    # recurrence using the last 15 read triples and a rolling sum of h-e.
    recent = deque()
    recent_kicks = 0
    for time in range(through + 1):
        z = state[-1]
        sign, loss = sign_and_loss(z)
        label = time % ORDER
        if len(recent) == ORDER:
            previous, previous_sign, previous_loss = recent.popleft()
            _, next_sign, next_loss = recent[0]
            remaining = (
                recent_kicks - (previous_loss - previous_sign)
                - (next_loss - next_sign)
            )
            assert z == previous - 3 * previous_sign + 2 * (
                previous_loss + next_loss
            ) + remaining, time
            recent_kicks -= previous_loss - previous_sign
        recent.append((z, sign, loss))
        recent_kicks += loss - sign

        reference_norm = max(reference_norm, max(state), -min(state))
        reference_span = max(reference_span, anchored_span(state))
        if first_terminal_visit is None and state in terminal_states:
            first_terminal_visit = time
        if time == arrival or time == through:
            assert state == terminal[0], time
        if time in checkpoint_times:
            checkpoints[time] = (state, cumulative_loss)

        if stage < ORDER:
            start, stop = intervals[stage]
            if time == stop:
                words.append(active_word)
                active_word = None
                stage += 1
            if stage < ORDER and time == intervals[stage][0]:
                active_word = {
                    "retired": retire_order.copy(),
                    "hasher": sha256(),
                    "head": [],
                    "length": 0,
                    "min_active_read": None,
                }
            if active_word is not None:
                if label in retirements:
                    token = z
                else:
                    assert z > 0, (stage + 1, time, label)
                    token = "+"
                    old_min = active_word["min_active_read"]
                    active_word["min_active_read"] = (
                        z if old_min is None else min(old_min, z)
                    )
                if len(active_word["head"]) < 30:
                    active_word["head"].append(token)
                if active_word["length"]:
                    active_word["hasher"].update(b",")
                active_word["hasher"].update(str(token).encode("ascii"))
                active_word["length"] += 1

        if label not in retirements and z <= 0:
            retirements[label] = (time, z, state)
            retire_order.append(label)
        if time < through:
            if time >= arrival:
                assert loss == 0, time
            cumulative_loss += loss
            state = root_step(state, sign, loss)

    assert first_terminal_visit is not None
    return {
        "checkpoints": checkpoints,
        "retirements": retirements,
        "retire_order": retire_order,
        "words": words,
        "reference_norm": reference_norm,
        "reference_span": reference_span,
        "reference_loss": cumulative_loss,
        "first_terminal_visit": first_terminal_visit,
    }


def check_stages(intervals, reference, expected):
    stages = []
    previous_retirement = -1
    for i, (start, stop) in enumerate(intervals):
        label = reference["retire_order"][i]
        exit_time, exit_value, exit_state = reference["retirements"][label]
        assert previous_retirement < start < stop <= exit_time
        previous_retirement = exit_time
        word = reference["words"][i]
        assert word["length"] == stop - start
        retired = reference["retire_order"][:i]
        assert retired == word["retired"]

        before, before_loss = reference["checkpoints"][start]
        after, after_loss = reference["checkpoints"][stop]
        before = in_read_labels(before, start)
        after = in_read_labels(after, stop)
        ray = [left - right for left, right in zip(before, after)]
        assert all(ray[j] == 0 for j in retired)
        assert min(ray) >= 0 and ray[label] > 0

        stage = {
            "label": label,
            "exit_time": exit_time,
            "exit_value": exit_value,
            "exit_state": list(exit_state),
            "word_interval": [start, stop],
            "word_sha256": word["hasher"].hexdigest(),
            "word_head": word["head"],
            "projection_labels": retired,
            "projection_start": [before[j] for j in retired],
            "word_loss": after_loss - before_loss,
            "drift": ray[label],
            "ray": ray,
        }
        assert stage == expected[i], (i + 1, stage, expected[i])
        stage["min_active_read"] = word["min_active_read"]
        stages.append(stage)
    return stages


def check_witness(cert, anatomy):
    source = tuple(cert["source"])
    intervals = cert["insertions"]
    arrival, through = cert["cycle"]
    assert cert["order"] == anatomy["order"] == len(source) == ORDER
    assert len(anatomy["stages"]) == len(intervals) == ORDER
    assert anatomy["source"] == list(source)
    assert anatomy["cycle_replay"] == [arrival, through]
    previous_stop = 0
    for start, stop in intervals:
        assert previous_stop <= start < stop <= arrival
        assert (stop - start) % ORDER == 0
        previous_stop = stop

    terminal, reads = cycle_at(tuple(anatomy["cycle_first_state"]))
    assert len(terminal) == through - arrival == anatomy["least_period"] == 30
    assert reads == anatomy["terminal_word"]
    canonical = min(reads[i:] + reads[:i] for i in range(len(reads)))
    assert canonical == anatomy["terminal_canonical_word"]
    assert all(sign_and_loss(reads[i])[0] + sign_and_loss(reads[i + ORDER])[0] == 0
               for i in range(ORDER))

    reference = replay_reference(source, intervals, arrival, through, terminal)
    for field in ("reference_norm", "reference_span", "reference_loss"):
        assert reference[field] == anatomy[field]
    labels = reference["retire_order"]
    assert labels == anatomy["retirement_labels"] == cert["retirement_labels"]
    stages = check_stages(intervals, reference, anatomy["stages"])
    rays = [stage["ray"] for stage in stages]
    matrix = [[ray[label] for ray in rays] for label in labels]
    determinant = abs(exact_determinant(matrix))
    assert determinant == prod(stage["drift"] for stage in stages)
    assert determinant == anatomy["determinant"] and anatomy["rank"] == ORDER
    assert anatomy["coefficient_floors"] == [0] * ORDER

    norm = reference["reference_norm"] + sum(max(ray) for ray in rays)
    bounds = {
        "steps_through_one_period": through + sum(stop - start for start, stop in intervals),
        "lifetime_norm": norm,
        "lifetime_displacement_from_source": 2 * norm,
        "lifetime_span": reference["reference_span"] + sum(anchored_span(ray) for ray in rays),
        "lifetime_loss": reference["reference_loss"] + sum(stage["word_loss"] for stage in stages),
        "radius_box_growth_numerator": 1,
        "radius_box_growth_denominator": (2 * norm)**ORDER,
        "radius_box_min_L": 2 * norm,
    }
    assert bounds == anatomy["bounds"]
    return {
        "certified_phase_time": arrival,
        "first_cycle_entry": reference["first_terminal_visit"],
        "cycle_first_state": terminal[0],
        "least_period": len(terminal),
        "terminal_word": reads,
        "reference_steps_checked": through,
        "reference_norm": reference["reference_norm"],
        "reference_span": reference["reference_span"],
        "reference_loss": reference["reference_loss"],
        "retirement_labels": labels,
        "stages": stages,
        "rank": ORDER,
        "absolute_determinant": determinant,
        "coefficient_floors": [0] * ORDER,
        "bounds": bounds,
    }


def divisors(value):
    return {
        factor for d in range(1, isqrt(value) + 1) if value % d == 0
        for factor in (d, value // d)
    }


def folded_orbits(n):
    k = n // 2
    slot_orbit = [0] * (k + 1)
    orbits = {}
    for start in range(1, k + 1):
        if slot_orbit[start] != 0:
            continue
        orbit = []
        j = start
        while not slot_orbit[j]:
            slot_orbit[j] = start
            orbit.append(j)
            j = min(2 * j, n - 2 * j)
        assert j == start
        orbits[start] = orbit
    slot_orbit[0] = 1
    return slot_orbit, orbits


def strict_length_screen():
    candidates = set.union(*(divisors(2**ORDER + sign) for sign in (-1, 1)))
    results = []
    excluded = []
    for n in sorted(n for n in candidates if n >= 7):
        k = n // 2
        r = next(i for i in range(1, ORDER + 1) if pow(2, i, n) in (1, n - 1))
        slot_orbit, orbits = folded_orbits(n)
        assert r == len(orbits[1])
        if r != ORDER:
            excluded.append({"length": n, "folded_order": r})
            continue

        if n == 99:
            # Integral primal certificate, indexed by least positive slot.
            offsets = {1: 0, 3: 1, 5: 1, 9: 1, 11: 0, 15: 2, 33: 3}
            assert offsets.keys() == orbits.keys()
        else:
            offsets = {j: int(j != 1) for j in orbits}
        weights = [offsets[j] for j in slot_orbit]
        assert weights[0] == 0
        assert all(weights[j] == weights[min(2 * j, n - 2 * j)]
                   for j in range(1, k + 1))

        literal = walk = 0
        margins = []
        for i in range(1, (k + 1) // 2 + 1):
            literal += weights[k + 1 - i] - weights[i - 1]
            # The quotient walk on a prime, and the same alternating
            # folded-orbit walk evaluated on slots at a composite.
            odd = 2 * i - 1
            even = 2 * i - 2
            walk += weights[min(odd, n - odd)] - weights[min(even, n - even)]
            assert literal == walk
            if i == 1:
                assert literal == 0
            else:
                assert literal > 0, (n, i, literal)
                margins.append(literal)

        # Shift toward the constant profile far enough to certify P9's
        # sufficient pairing margin (m + 1) min(pi) > 1 exactly.
        m = k + 1
        shift = 1 + sum(weights)
        normalizer = m * shift + sum(weights)
        pairing_numerator = (m + 1) * shift - normalizer
        assert min(weights) == 0 and pairing_numerator == 1
        suffix = radius_numerator = 0
        for i in range(k, 0, -1):
            suffix += weights[i] + shift
            radius_numerator = max(radius_numerator, abs(m * suffix - (m - i) * normalizer))
        assert radius_numerator > 0
        results.append({
            "length": n,
            "folded_order": r,
            "folded_return_sign": 1 if pow(2, r, n) == 1 else -1,
            "orbit_count": len(orbits),
            "root_orbit": orbits[1],
            "secondary_count": len(margins),
            "min_integer_secondary_slack": min(margins),
            "max_integer_secondary_slack": max(margins),
            "profile_common_shift": shift,
            "profile_normalizer": normalizer,
            "profile_pairing_margin_numerator": pairing_numerator,
            "profile_radius_numerator": radius_numerator,
            "profile_radius_denominator": m * normalizer,
            "offsets": offsets if n == 99 else "zero on the root and one off it",
        })
    return {"strict_candidates": results, "lower_order_divisors": excluded}


def frozen_rows(input_directory, name):
    data = (input_directory / name).read_bytes()
    assert sha256(data).hexdigest() == INPUT_HASHES[name], name
    return [json.loads(line) for line in data.splitlines()]


def main():
    if not __debug__:
        raise SystemExit("Run without -O so that all exact checks execute.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_directory", type=Path)
    args = parser.parse_args()
    certs, anatomies = [frozen_rows(args.input_directory, name) for name in INPUT_HASHES]
    assert len(certs) == len(anatomies) == 4
    witness_results = [
        check_witness(cert, anatomy)
        for cert, anatomy in zip(certs, anatomies)
    ]
    assert all(witness["terminal_word"] == witness_results[0]["terminal_word"]
               for witness in witness_results)
    output = {
        "input_hashes": INPUT_HASHES,
        "witnesses": witness_results,
        "length_screen": strict_length_screen(),
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
