"""Certify an uncompressed source for the index-11 tangent torsion."""
from fractions import Fraction
import json
from math import lcm
from pathlib import Path
import re
import subprocess
import tempfile

from ingress import (ACCESS, ID, START, Poly, build, centered_row, leading_chart,
                     line_data, race_domain, reply, within_interlevel, xchart)
from pi_geometry import dirichlet_mass


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def small_points(vertices):
    """Queries near the interior barycentre, with controlled denominators."""
    middle = tuple(sum(v[i] for v in vertices) / len(vertices) for i in range(6))
    for power in (24, 29, 34, 37):
        denominator = 2**power
        first = [Fraction(round(coordinate * denominator), denominator)
                 for coordinate in middle[:5]]
        yield (*first, 1 - sum(first))


def represent(chart, gaps):
    source = (*xchart(gaps), 1)
    point = tuple(dot(row, source) for row in leading_chart(chart))
    denominator = lcm(*(x.denominator for x in point))
    return (denominator, *(int(x * denominator) for x in point))


def reconstruct(reader, vertices):
    desc, inter = within_interlevel()
    source = Poly.root().restrict(sorted(inter) + sorted(desc))
    chart = START
    sheet = (0, 6, 5, *(0,) * 10, *ID)
    results = []
    target_history = (4, 3, 2, 2, 1)
    for remainder, next_rank in zip((6, 0, 36, 18, 18), target_history):
        for gaps in small_points(vertices):
            point = represent(chart, gaps)
            if max(map(abs, point)) >= 2**63:
                continue
            kind, details, exits, guards = reply(reader, point, chart, sheet, 500)
            if kind == "stage" and all(dot(row, v) >= 0 for row in guards for v in vertices):
                break
        else:
            raise AssertionError("no guard chart containing the candidate vertices")
        deadlines = tuple(row for row, *_ in exits)
        supported = []
        for winner, descriptor, target, children in exits:
            cuts = [centered_row(tuple(a - b for a, b in zip(rival, winner)))
                    for rival in deadlines if rival != winner]
            if all(dot(row, vertex) >= 0 for row in cuts for vertex in vertices):
                supported.append((winner, descriptor, target, children))
        assert len(supported) == 1
        winner, descriptor, target, children = supported[0]
        source = race_domain(source.restrict(sorted(guards)), deadlines, winner)
        assert source is not None
        candidates = [(child, index) for child, rem, index in children
                      if rem == remainder and child[:3] == (0, 6, next_rank)]
        assert len(candidates) == 1
        sheet, index = candidates[0]
        chart = target
        results.append(dict(rank=next_rank, shorts=details[0], word_length=details[2],
                            unpaired_phases=details[3], remainder=remainder, index=index,
                            decrement=descriptor[0], guard_base_coefficient=descriptor[-2]))
    index = 1
    for edge in results:
        index *= edge["index"]
    period, _, _, _ = line_data(chart, sheet)
    assert period == 66
    assert results[2]["guard_base_coefficient"] == -1 and results[2]["decrement"] == 66
    assert results[2]["index"] == 11 and index == 5544
    source_mass = dirichlet_mass(source)
    mass = source_mass / index
    assert mass > Fraction(418, 10**15)
    return dict(edges=results, source_vertices=len(source.vertices), source_mass=str(source_mass),
                lane_index=index, direct_source_mass=str(mass), centered_period=period,
                arrival=sheet)


def main():
    archived = (ACCESS / "class_compat" / "period66-witness.out").read_text()
    vertex_data = archived.split(" vertices ", 1)[1].split(" chart ", 1)[0]
    flat = [Fraction(int(a), int(b)) for a, b in re.findall(r"Fraction\((-?\d+), (\d+)\)", vertex_data)]
    assert len(flat) % 6 == 0
    vertices = [tuple(flat[i:i+6]) for i in range(0, len(flat), 6)]
    with tempfile.TemporaryDirectory(prefix="kaprekar-task147-witness-") as tmp:
        binary = build(Path(tmp), paired_fast=False)
        with subprocess.Popen([str(binary)], stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, text=True) as reader:
            result = reconstruct(reader, vertices)
            reader.stdin.close()
            assert reader.wait() == 0
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
