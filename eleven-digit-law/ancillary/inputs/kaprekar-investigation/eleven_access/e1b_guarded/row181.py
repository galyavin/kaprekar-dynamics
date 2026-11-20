"""Replay the obstruction on all cyclic phases of the remaining unpaired row."""

from gzip import open as gzip_open
import json
from pathlib import Path
import sys

from anchor import constrained_counts, gap_order_profiles, majorization_guards
from covectors import add_fixed_constraints
from guarded import ExactCones, I, branch_word, dot, family, load_rows
from guarded import normal_basis, null_basis, phases_for
from repeated import replay_ray


def main(path):
    row = load_rows()[181]
    phases, tangents = phases_for(row)
    _, _, _, points = family(row)
    period = len(phases)
    assert period == 3 and row["dimension"] == 1
    with gzip_open(path, "rt") as source:
        evidence = json.load(source)["181"]
    repeats = evidence["repeated_words"]
    assert all(ray["type"] == "contracting" for ray in repeats)
    labels = tuple(tuple(ray["labels"]) for ray in repeats)
    assert labels == ((0, 1, 0), (1, 1, 0), (1, 1, 1))

    cones = ExactCones()
    for phase in range(period):
        span = normal_basis(tangents[phase])
        for word, ray in zip(labels, repeats):
            rotated = word[phase:] + word[:phase]
            W, guards = branch_word(phases, phase, rotated)
            if phase == 0:
                replay_ray(cones, W, guards, ray)
            # Shift the strict contracting eigenray by the same phase prefix.
            span = add_fixed_constraints(span, W)
        assert not null_basis(span)
    print("Row 181 contracting return words (exact labels):", labels)
    print("Row 181 upper dimensions on cyclic phases (exact):", [0, 0, 0])

    for word, coordinate in zip(labels, (0, 0, 4)):
        W, _ = branch_word(phases, 0, word)
        assert W[coordinate] == I[coordinate]
        assert all(v[coordinate] == 0 for v in tangents[0])
    assert points[0][0] * 8 == 7 and points[0][4] * 4 == 1
    print("Fixed phase-zero hyperplanes on the three convergent repetitions (exact):")
    print("q1=7/8, q1=7/8, q5=1/4")

    summary = []
    for parameter, profile in gap_order_profiles(tangents, points[:period]):
        constraints = []
        for q, T in zip(profile, tangents):
            strict, weak = [], []
            for normal in majorization_guards(q):
                (strict if all(dot(normal, t) == 0 for t in T) else weak).append(normal)
            constraints.append((tuple(strict), tuple(weak)))
        counts, exhausted = constrained_counts(cones, phases, 0, constraints, 15)
        assert exhausted and (counts[3], counts[9], counts[15]) == (8, 70, 310)
        summary.append((str(parameter), counts[3], counts[9], counts[15]))
    assert len(summary) == 3
    print("Row 181 strict prefixes at primitive depths 3, 9, 15 (exact):", summary)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
