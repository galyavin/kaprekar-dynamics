"""Restrict the established eleven-digit periodic catalogue to strict pairing."""
from collections import Counter
from fractions import Fraction
from gzip import open as gzip_open
import json
from pathlib import Path

from lp import family, reduce_rows


def leading(q):
    digits = sorted((1, *q, *(1 - a for a in q)), reverse=True)
    return tuple(digits[i] - digits[-1 - i] for i in range(5))


def paired(q):
    levels = tuple((q[i], 1 - q[4 - i]) for i in range(5))
    return (
        0 < q[4] and max(levels[0]) < 1
        and all(min(high) > max(low) for high, low in zip(levels, levels[1:]))
    )


def signs(q):
    return "".join("+" if a + b > 1 else "-" if a + b < 1 else "0"
                   for a, b in zip(q[:3], q[::-1][:3]))


def output():
    with gzip_open(Path("evidence/energy.json.gz"), "rt") as source:
        established = json.load(source)
    rows = []
    for row in established:
        if not row["in_ball"]:
            continue
        faces = tuple((middle, tuple(map(tuple, groups)), edge)
                      for middle, groups, edge in row["faces"])
        _, equations, strict = family(faces, 5)
        _, pivots = reduce_rows([*equations, (1,) * 6], 6)
        dim = 6 - len(pivots)
        start = tuple(map(Fraction, row["q_witness"]))
        phases = [start]
        for _ in range(len(faces) - 1):
            phases.append(leading(phases[-1]))
        assert leading(phases[-1]) == start
        all_paired = all(map(paired, phases))
        word = tuple(map(signs, phases)) if all_paired else ()
        if all_paired:
            assert dim <= 5
        rows.append({
            "dimension": dim,
            "period": row["witness_period"],
            "paired": all_paired,
            "word": word,
            "witness": list(map(str, start)),
            "faces": row["faces"],
        })
    assert len(rows) == 191
    print("All in-ball cells by (paired, dimension, period)")
    print(sorted(Counter((r["paired"], r["dimension"], r["period"]) for r in rows).items()))
    print("Paired words by (dimension, period, word)")
    for row in sorted(rows, key=lambda r: (-r["dimension"], r["period"], r["word"])):
        if row["paired"]:
            print(row["dimension"], row["period"], " ".join(row["word"]), " ".join(row["witness"]))
    Path("catalogue.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    output()
