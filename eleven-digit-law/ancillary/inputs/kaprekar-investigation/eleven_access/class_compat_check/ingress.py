"""Read the Task 137 physical chart oracle and independently classify its inherited exits.

The polytope and plane-coordinate routines are Task 137 infrastructure.
The worklist, arithmetic inventory, physical normal key and terminal race
coverage are reconstructed here. No Task 143/144 program is executed.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
import fcntl
import json
from math import gcd, lcm
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = ACCESS.parents[1]
sys.path.insert(0, str(ACCESS / "z1_stages"))
from cover import (START, ID, centered_row, exact_point, integer_vertices, padded_form, project_plane,
                   record_stage_probes, reply)
from pi_geometry import Poly, within_interlevel
from geometry import xchart


@lru_cache(maxsize=16384)
def guard_cut(poly, guard):
    return poly.cut(guard)


def invert(row):
    return tuple(-value for value in row)


def sections(poly, guards):
    """Disjoint interiors split by the first failed primitive comparison."""
    for guard in sorted(guards):
        if min(sum(a * b for a, b in zip(guard, vertex))
               for vertex in integer_vertices(poly.vertices)) >= 0:
            continue
        outside = guard_cut(poly, invert(guard))
        accepted = guard_cut(poly, guard)
        assert accepted is not None
        if outside is not None:
            yield outside, False
        poly = accepted
    yield poly, True


@lru_cache(maxsize=8192)
def race_domain(domain, deadlines, winner):
    for rival in deadlines:
        if rival != winner:
            domain = guard_cut(domain, centered_row(tuple(a - b for a, b in zip(rival, winner))))
            if domain is None:
                break
    return domain


def leading_chart(chart):
    return tuple(tuple(Fraction(value, chart[i + 6]) for value in chart[i:i + 6])
                 for i in range(0, 35, 7))


def zero_chart(chart):
    return all(not value for i, value in enumerate(chart) if i % 7 != 6)


def line_data(chart, sheet):
    _, stride, rank, *coordinates = sheet
    assert rank == 1
    slope = coordinates[:5]
    point = coordinates[5:10]
    column = coordinates[10:]
    index = gcd(*column)
    direction = tuple(entry // index for entry in column)
    axis = next(i for i, entry in enumerate(direction) if entry)
    assert direction[axis] > 0 and index % 6 == 0
    parameter_slope = Fraction(slope[axis], direction[axis])
    assert parameter_slope.denominator == 1
    assert all(s == int(parameter_slope) * u for s, u in zip(slope, direction))
    reference = leading_chart(chart)
    for u, expression in zip(direction, reference):
        assert all(direction[axis] * a == u * b for a, b in zip(expression, reference[axis]))
    # A physical fixed-base saturated line has tangent step 6*u. Reduce by its
    # first nonzero coordinate. This is an equality test in the integral quotient.
    quotient = tuple(x - (point[axis] // (6 * direction[axis])) * 6 * u
                     for x, u in zip(point, direction))
    centered_period = stride * index // gcd(index, int(parameter_slope))
    return centered_period, (direction, quotient), reference[axis], axis


def projected_interval(poly, function):
    *coefficients, constant = function
    values = [sum(x * y for x, y in zip(coefficients, xchart(gaps))) + constant
              for gaps in poly.vertices]
    return min(values), max(values)


def union(intervals):
    current = []
    for lo, hi in sorted(intervals):
        assert lo < hi
        if current and lo <= current[-1][1]:
            current[-1] = current[-1][0], max(hi, current[-1][1])
        else:
            current.append((lo, hi))
    return current


def line_chart(direction, axis, interval):
    lo, hi = interval
    return tuple(item
                 for u in direction
                 for item in padded_form((-Fraction(u, direction[axis]) * (hi-lo) / 6,
                                           0, 0, 0, 0,
                                           Fraction(u, direction[axis]) * (lo + (hi-lo) / 6))))


def inspect(reader, poly, chart, sheet, rng):
    for attempt in range(20):
        point = exact_point(poly, chart, rng, attempt)
        if max(map(abs, point)) < 2**63:
            outcome = reply(reader, point, chart, sheet, 500)
            if outcome[0] != "ambiguous":
                assert outcome[0] in ("stage", "terminal"), outcome[0:2]
                return outcome
    raise RuntimeError(f"no exact witness in the interior of {sheet=}")


def split_permanence(direction, axis, intervals):
    """Give separate canonical outer intervals and widened inner half-rays."""
    slope = tuple(Fraction(x, direction[axis]) for x in direction)
    steps = [b-a for a, b in zip((0, *slope), (*slope, 0))]
    left = -1 / (7 * max(-a for a in steps))
    right = 1 / (7 * max(steps))
    regions = defaultdict(list)
    for lo, hi in intervals:
        if lo < left:
            regions["outer"].append((lo, min(hi, left)))
        if hi > right:
            regions["outer"].append((max(lo, right), hi))
        if max(lo, left) < min(hi, 0):
            regions["inner-negative"] = [(left, Fraction(0))]
        if max(lo, 0) < min(hi, right):
            regions["inner-positive"] = [(Fraction(0), right)]
    return regions.values()


def certify_prime_returns(reader, root, sections_by_key, canonical):
    """Cover the inherited intervals or the canonical permanence overcover."""
    rng = random.Random(148151)
    simplex = Poly.root()
    tally = Counter()
    section_count = 0
    for (direction, anchor), raw_intervals in sorted(sections_by_key.items()):
        axis = next(i for i, value in enumerate(direction) if value)
        tangent = tuple(6 * value for value in direction)
        sheet = (root, 6, 1, 0, 0, 0, 0, 0, *anchor, *tangent)
        regions = split_permanence(direction, axis, raw_intervals) if canonical else (raw_intervals,)
        section_count += len(regions)
        for interval in (part for group in regions for part in union(group)):
            assert not interval[0] < 0 < interval[1]
            chart = line_chart(direction, axis, interval)
            pending = [simplex]
            while pending:
                poly = pending.pop()
                kind, details, options, guards = inspect(reader, poly, chart, sheet, rng)
                assert kind == "stage" and details[3:] == (0, 1), (root, anchor, interval, kind, details)
                for part, accepted in sections(poly, guards):
                    if not accepted:
                        pending.append(part)
                        continue
                    deadlines = tuple(row for row, *_ in options)
                    for winner, descriptor, target, children in options:
                        if race_domain(part, deadlines, winner) is None:
                            continue
                        decrement, *_, coefficient, _ = descriptor
                        assert coefficient == 0 and zero_chart(target) and not children
                        tally[decrement] += 1
    return section_count, dict(tally)


def audit_root(reader, root, canonical, restrict_prime_index):
    desc, inter = within_interlevel()
    chamber = Poly.root().restrict(sorted(inter) + sorted(desc))
    assert len(chamber.vertices) == 77
    start = (root, 6, 5, *(0,) * 5, *((i * root) % 6 for i in range(1, 6)), *ID)
    work = [(chamber, START, start, (5,), False)]
    rng = random.Random(137049)
    covered = defaultdict(list)
    seen = set()
    simplex = Poly.root()
    reads = Counter()
    arrivals = Counter()
    prime_history = Counter()
    prime_sections = defaultdict(list)
    nreads = 0
    while work:
        poly, chart, sheet, ranks, compressed = work.pop()
        if ranks[-1] == 2 and not compressed:
            poly, chart, sheet, faces = project_plane(poly, chart, sheet, simplex, False)
            key = sheet, ranks, chart
            region = covered[key]
            vertices = integer_vertices(poly.vertices)
            if any(all(all(sum(a * b for a, b in zip(face, v)) >= 0 for v in vertices)
                       for face in old) for old in region):
                continue
            region.append(faces)
            compressed = True
        record = poly, chart, sheet, ranks
        if record in seen:
            continue
        seen.add(record)
        reads[ranks[-1], sheet[1]] += 1
        kind, details, options, guards = inspect(reader, poly, chart, sheet, rng)
        nreads += 1
        if nreads % 5000 == 0:
            print(f"root {root}: {nreads} physical chart reads", file=sys.stderr, flush=True)
        for piece, accepted in sections(poly, guards):
            if not accepted:
                work.append((piece, chart, sheet, ranks, compressed))
                continue
            if kind == "terminal":
                assert zero_chart(chart)
                continue
            assert details[-1] == ranks[-1]
            deadlines = tuple(row for row, *_ in options)
            for winner, _, target, children in options:
                domain = race_domain(piece, deadlines, winner)
                if domain is None:
                    continue
                by_child = {child: index for child, _, index in children}
                for child, index in by_child.items():
                    rank = child[2]
                    path = (*ranks, rank)
                    if zero_chart(target):
                        continue
                    assert rank >= 1
                    if restrict_prime_index and rank == 2 and (
                        (ranks[-1] == 3 and index % 11 != 0) or len(path) > 5
                    ):
                        continue
                    if rank > 1:
                        work.append((domain, target, child, path, False))
                        continue
                    assert child[0] % 6 == root
                    assert all((a - i * root) % 6 == 0 for i, a in enumerate(child[8:13], 1))
                    period, key, parameter, _ = line_data(target, child)
                    arrivals[child[1], period] += 1
                    if period == 66:
                        prime_history["".join(map(str, path))] += 1
                        prime_sections[key].append(projected_interval(domain, parameter))
    print(f"root {root}: closed after {nreads} reads, checking terminal returns",
          file=sys.stderr, flush=True)
    section_count, terminal = certify_prime_returns(reader, root, prime_sections, canonical)
    result = {
        "root": root,
        "read_counts": [[*key, count] for key, count in sorted(reads.items())],
        "arrival_counts": [[*key, count] for key, count in sorted(arrivals.items())],
        "period_66_history": dict(prime_history),
        "period_66_normal_lines": len(prime_sections),
        "period_66_terminal_returns": terminal,
    }
    if canonical:
        result["period_66_sections"] = section_count
    return result


def run_root(binary, root, canonical, restrict_prime_index):
    with subprocess.Popen([str(binary)], stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE, text=True) as reader:
        result = audit_root(reader, root, canonical, restrict_prime_index)
        reader.stdin.close()
        assert reader.wait() == 0
    return result


def build(snapshot, paired_fast):
    for source in (ACCESS / "z1_stages" / "chart.cpp",
                   ACCESS / "stage_graph" / "lottery.cpp",
                   ACCESS / "stage_graph" / "graph.cpp"):
        shutil.copyfile(source, snapshot / source.name)
    record_stage_probes(snapshot / "graph.cpp", paired_fast=paired_fast)
    subprocess.run(["c++", "-O3", "-std=c++17", "-Wno-return-type", "chart.cpp",
                    "-o", "chart"], cwd=snapshot, check=True)
    return snapshot / "chart"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roots", nargs="+", type=int, choices=range(6), default=list(range(6)))
    parser.add_argument("--canonical", action="store_true", help="also use the canonical permanence overcover")
    parser.add_argument("--index-11-at-rank-three", action="store_true",
                        help="inspect only paths with one rank-two keep and a multiple-of-11 third index")
    args = parser.parse_args()
    with (ROOT / ".task-monitor" / "compute.lock").open() as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix="kaprekar-task147-") as tmp:
            binary = build(Path(tmp), paired_fast=True)
            with ProcessPoolExecutor(max_workers=len(args.roots)) as pool:
                jobs = [pool.submit(run_root, binary, root, args.canonical, args.index_11_at_rank_three)
                        for root in args.roots]
                payload = [job.result() for job in jobs]
            print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
