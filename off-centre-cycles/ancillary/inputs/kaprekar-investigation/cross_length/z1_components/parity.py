"""Generate the frozen draw streams and test exact simultaneous binomial tails."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path
import random
import subprocess
import sys

SAMPLES = 512
SCALES = (20_000, 200_000)
ALPHA = Fraction(1, 1440)
RULES = {
    **{r: ("above", Fraction(1, 10)) for r in (3, 4, 5, 7, 8)},
    **{r: ("above", Fraction(1, 50)) for r in (9, 11, 13)},
    6: ("below", Fraction(13, 200)),
    **{r: ("below", Fraction(1, 40)) for r in (10, 12)},
    **{r: ("below", Fraction(2, 25)) for r in range(14, 21)},
}


def tails(n, threshold):
    """Cumulative integer binomial weights with common denominator b**n."""
    a, b = threshold.as_integer_ratio()
    cumulative = [0]
    for k in range(n + 1):
        cumulative.append(cumulative[-1] + comb(n, k) * a**k * (b - a)**(n - k))
    assert cumulative[-1] == b**n
    return cumulative, b**n


def row(binary, order, scale_index):
    scale = SCALES[scale_index - 1]
    seed = 17620261120 + 100 * order + scale_index
    generator = random.Random(seed)
    inputs = [str(SAMPLES) + " 1000000"]
    for _ in range(SAMPLES):
        inputs.append(" ".join(str(generator.randrange(scale)) for _ in range(order)))
    panel = subprocess.run(
        [binary, "panel", str(order)], input="\n".join(inputs) + "\n",
        capture_output=True, text=True, check=True,
    )
    records = [tuple(map(int, line.split())) for line in panel.stdout.splitlines()]
    assert len(records) == SAMPLES
    marks = Counter(record[0] for record in records)
    hits, misses, unknown = marks[1], marks[0], marks[-1]
    assert SAMPLES == hits + misses + unknown
    direction, threshold = RULES[order]
    cumulative, denominator = tails(SAMPLES, threshold)
    above = Fraction(denominator - cumulative[hits], denominator) <= ALPHA
    below = Fraction(cumulative[hits + unknown + 1], denominator) <= ALPHA
    assert not (above and below)
    marked = ("above" if above else "below") if above or below else "unresolved"
    decision = "confirmed" if marked == direction else marked
    if marked != direction and marked != "unresolved":
        decision = "contradicted"
    need_above = min(k for k in range(SAMPLES + 1)
                     if Fraction(denominator - cumulative[k], denominator) <= ALPHA)
    need_below = max(k for k in range(SAMPLES + 1)
                     if Fraction(cumulative[k + 1], denominator) <= ALPHA)
    return (order, scale, seed, hits, misses, unknown, direction, str(threshold), need_above,
            need_below, decision, max(record[1] for record in records),
            max(record[2] for record in records), max(record[4] for record in records),
            sha256(panel.stdout.encode()).hexdigest())


def main():
    binary = str(Path(sys.argv[1]).resolve())
    print("order\tscale\tseed\tzero\tpositive\tunknown\tprediction\tthreshold\t"
          "min_zero_for_above\tmax_possible_zero_for_below\tdecision\t"
          "max_explicit_updates\tmax_clock\tmax_projected_search\tmarks_sha256", flush=True)
    for order in range(3, 21):
        for scale_index in range(1, 3):
            print(*row(binary, order, scale_index), sep="\t", flush=True)


if __name__ == "__main__":
    main()
