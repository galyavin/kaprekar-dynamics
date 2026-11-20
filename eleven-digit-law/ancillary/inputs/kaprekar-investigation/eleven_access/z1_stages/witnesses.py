"""Certify positive original mass for low-rank histories on open cells."""
from argparse import ArgumentParser
from ast import literal_eval
from fractions import Fraction as Q
from math import factorial, lcm, prod
from pathlib import Path
from random import Random
import re
from subprocess import PIPE, Popen

from cover import ID, START, centered_row, reply, within_interlevel, xchart


GRID = 2**38
CONTROLS = (
    (48000, (5, 4, 3, 2, 2, 1, 0), (213133136167, 190311175378, 161941236035, 112373615996, 15532002865)),
    (48002, (5, 4, 3, 2, 2, 1, 1, 0), (256853830325, 142659989201, 102881701575, 56341069472, 48365276043)),
    (48000, (5, 4, 3, 2, 1, 0), (198026272578, 156123904658, 101767147319, 66916358067, 32315112159)),
    (48000, (5, 4, 3, 2, 1, 1, 0), (215574821786, 138892895160, 109803020123, 91024528535, 13999268143)),
    (48000, (5, 4, 3, 2, 1, 1, 1, 0), (265311890134, 243417445419, 172454829105, 156380545639, 143711481260)),
)


def gaps(q):
    return tuple(Q(a - b, GRID) for a, b in zip((GRID, *q), (*q, 0)))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def entrance(q):
    signs, pairing = within_interlevel()
    for n in range(501):
        current = gaps(q)
        if min(current) > 0 and all(dot(row, current) > 0 for row in (*signs, *pairing)):
            return n, q, current
        labels = sorted((GRID, *q, *(GRID - v for v in q)))
        q = tuple(labels[10 - i] - labels[i] for i in range(5))
    raise ValueError("control missed the strict chamber")


def at(chart, x):
    return tuple((dot(chart[i:i + 5], x) + chart[i + 5]) / chart[i + 6] for i in range(0, 35, 7))


def restrict_box(rows, center, radius):
    for row in rows:
        gain = dot(row, center)
        norm = sum(abs(a - row[-1]) for a in row[:-1])
        assert gain > 0
        if norm:
            radius = min(radius, gain / (2 * norm))
    return radius


def one(reader, base, center, rng, prefix=()):
    x = xchart(center)
    within, inter = within_interlevel()
    radius = restrict_box((*within, *inter), center, min(center) / 10)
    root = base % 6
    sheet = (root, 6, 5, *(0,) * 5, *(((i - 5) * root) % 6 for i in range(5)), *ID)
    chart = START
    ranks = [5]
    lane = Q(1)
    history = []
    while len(ranks) < 10:
        point = at(chart, x)
        scale = lcm(*(c.denominator for c in point))
        phase, details, exits, guards = reply(reader, (scale, *(int(scale * c) for c in point)), chart, sheet, 500)
        if phase != "stage":
            return None
        radius = restrict_box(guards, center, radius)
        ordered = sorted(exits, key=lambda option: dot(option[0][:-1], x) + option[0][-1])
        lead, stage, chart, children = ordered[0]
        differences = [centered_row(tuple(a - b for a, b in zip(competitor[0], lead)))
                       for competitor in ordered[1:]]
        if any(dot(row, center) == 0 for row in differences):
            return None
        radius = restrict_box(differences, center, radius)
        if not children:
            assert len(history) >= len(prefix)
            assert all(c == 0 for i, c in enumerate(chart) if i % 7 != 6)
            history.append((details, stage, None, None))
            ranks.append(0)
            mass = factorial(11) * (2 * radius) ** 5 * prod(c / 2 for c in center) * lane
            return tuple(ranks), radius, lane, mass, history
        candidates = [(child, rem, index) for child, rem, index in children if base % child[1] == child[0]]
        if len(history) < len(prefix):
            incoming, outgoing, required, index, descriptor = prefix[len(history)]
            assert sheet[:2] == incoming and (*details, *stage) == descriptor
            candidates = [(child, rem, stride) for child, rem, stride in candidates
                          if child[:2] == outgoing and rem == required and stride == index]
        assert candidates
        sheet, remainder, index = rng.choice(candidates)
        history.append((details, stage, remainder, sheet))
        ranks.append(sheet[2])
        lane /= index
    return None


def show(answer):
    ranks, radius, weight, mass, history = answer
    print("ranks", ranks, "half-width", radius, "lane weight", weight, "original mass bound", mass)
    for record in history:
        print("stage", *record)


def check_third_keep(reader, rng):
    here = Path(__file__).resolve().parent
    evidence = (here / "third-keep-cell.out").read_text()
    vertices_text = evidence.split(" vertices ")[1].split(" chart ")[0]
    rows = [Q(int(a), int(b)) for a, b in re.findall(r"Fraction\((-?\d+), (\d+)\)", vertices_text)]
    assert len(rows) > 0 and len(rows) % 6 == 0
    center = tuple(sum(rows[i::6]) / (len(rows) // 6) for i in range(6))
    prefix = literal_eval((here / "third-keep-prefix.txt").read_text())
    answer = one(reader, 48006, center, rng, prefix)
    assert answer is not None and answer[0] == (5, 4, 3, 2, 2, 2, 1, 1, 0)
    print("base residue 6 modulo 12; original gap centre", *(str(x) for x in center))
    show(answer)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("reader", type=Path, help="compiled exact chart reader in a private snapshot")
    parser.add_argument("--third-keep", action="store_true")
    args = parser.parse_args()
    rng = Random(137097)
    with Popen([str(args.reader.resolve())], stdin=PIPE, stdout=PIPE, text=True) as reader:
        if args.third_keep:
            check_third_keep(reader, rng)
        else:
            for base, target, original in CONTROLS:
                n, regular, center = entrance(original)
                for trial in range(250):
                    answer = one(reader, base, center, rng)
                    if answer is None or answer[0] != target:
                        continue
                    print("base", base, "entrance", n, "grid", GRID, "q", regular, "trial", trial)
                    show(answer)
                    break
                else:
                    raise RuntimeError(f"target history not found for {base}")
        reader.stdin.close()
        assert reader.wait() == 0


if __name__ == "__main__":
    main()
