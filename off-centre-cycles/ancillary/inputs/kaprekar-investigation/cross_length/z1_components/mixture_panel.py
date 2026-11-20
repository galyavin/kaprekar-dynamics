"""Measure last-active component/residue arrivals on fresh uniform order-8 cubes."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random
import subprocess
import sys

N = 8192
SCALES = (100_000, 1_000_000)
GRID = 1_000
ALPHA = Fraction(1, 800)


def cumulative(n, k, a, b):
    """Exact probability of at most k successes under Bernoulli(a/b)."""
    if k < 0:
        return Fraction(0)
    if k >= n:
        return Fraction(1)
    if k > n // 2:
        return 1 - cumulative(n, n - k - 1, b - a, b)
    if a == 0:
        return Fraction(1)
    if a == b:
        return Fraction(0)
    adverse = b - a
    term = total = adverse**n
    for j in range(1, k + 1):
        term = term * (n - j + 1) * a // (j * adverse)
        total += term
    return Fraction(total, b**n)


def confidence_interval(count):
    """Exact simultaneous binomial interval rounded outward to 1/GRID."""
    def last_true(predicate):
        left, right = 0, GRID + 1
        while right - left > 1:
            center = (left + right) // 2
            if predicate(center):
                left = center
            else:
                right = center
        return left

    lower = last_true(lambda a: cumulative(N, N - count, GRID - a, GRID) <= ALPHA)
    upper_flip = last_true(lambda a: cumulative(N, count, GRID - a, GRID) <= ALPHA)
    return [str(Fraction(lower, GRID)), str(1 - Fraction(upper_flip, GRID))]


def panel(binary, components, scale):
    rng = random.Random(176080000 + scale)
    draws = [f"{N} 1000000"] + [
        " ".join(str(rng.randrange(scale)) for _ in range(8)) for _ in range(N)
    ]
    data = subprocess.run([binary, "mix", "8"], text=True, check=True, capture_output=True,
                          input="\n".join(draws) + "\n").stdout
    known = {
        (c["period"], c["descent"], tuple(c["canonical_section"])): c for c in components
    }
    arrival = Counter()
    hits = Counter()
    phases = {c["name"]: Counter() for c in components}
    others = Counter()
    for row in data.splitlines():
        mark, _, _, _, _, period, descent, residue, *section = map(int, row.split())
        assert mark in (0, 1)  # All sampled terminal marks resolved under the cap.
        key = period, descent, tuple(section)
        if key in known:
            component = known[key]
            name = component["name"]
            assert mark == (residue in component["good_residues"])
            phases[name][residue] += 1
        elif period:
            name = "other_verified_return"
            others[key] += 1
        else:
            name = "no_verified_return"
        arrival[name] += 1
        hits[name] += mark
    assert sum(arrival.values()) == N
    counts = {
        **{f"arrivals_{key}": arrival[key] for key in
           ("period1712", "period368", "other_verified_return", "no_verified_return")},
        **{f"zero_loss_from_{key}": hits[key] for key in ("period1712", "period368")},
        "zero_loss_from_remaining": hits["other_verified_return"] + hits["no_verified_return"],
        "zero_loss_total": sum(hits.values()),
        "possible_arrivals_period1712": arrival["period1712"] + arrival["no_verified_return"],
        "possible_arrivals_period368": arrival["period368"] + arrival["no_verified_return"],
    }
    assert len(counts) == 10

    return {
        "L": scale, "draws": N, "seed": 176080000 + scale,
        "exact_sample_counts": counts,
        "simultaneous_95_percent_intervals": {
            key: confidence_interval(count) for key, count in counts.items()
        },
        "named_component_residue_counts": {
            key: dict(sorted(counts.items())) for key, counts in phases.items()
        },
        "other_components": [
            {"period": key[0], "descent": key[1], "canonical_section": key[2], "arrivals": count}
            for key, count in others.most_common()
        ],
        "raw_classifications_sha256": sha256(data.encode()).hexdigest(),
    }


def main(binary, weights_path):
    assert 2 * 10 * len(SCALES) * ALPHA == Fraction(1, 20)
    assert cumulative(5, 2, 1, 2) == sum(comb(5, k) for k in range(3)) / Fraction(32)
    weights = json.loads(weights_path.read_text())
    panels = [panel(binary, weights["components"], scale) for scale in SCALES]
    print(json.dumps(panels, indent=2))


if __name__ == "__main__":
    main(str(Path(sys.argv[1]).resolve()), Path(sys.argv[2]))
