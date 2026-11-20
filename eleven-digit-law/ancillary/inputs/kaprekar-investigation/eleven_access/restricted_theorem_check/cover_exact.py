"""Independent exact chart coverage.  The reader certifies whole paths, the driver covers domains."""
import argparse
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import fcntl
import random
import shutil
import subprocess
import sys
import tempfile

from exact_model import CENTRE, IDENTITY, Sheet, ZERO, bezout, combine, dot, exit_chart, leading, radius

HORIZON = 120
SHORT_CAP = 120
CENTER = (0, 0, 0, 0, 0, 1)
ROOT = Path(__file__).resolve().parents[3]


def normalized(row):
    den = lcm(*(Q(x).denominator for x in row))
    nums = tuple(int(den * x) for x in row)
    common = gcd(*nums)
    return tuple(x // common for x in nums) if common else nums


@dataclass(frozen=True)
class Poly:
    dimension: int
    vertices: tuple
    faces: tuple

    def clip(self, row):
        numerators = [dot(row, v) for v in integer_points(self.vertices)]
        if all(v >= 0 for v in numerators):
            return self
        kept = [p for p, v in enumerate(numerators) if v >= 0]
        if all(numerators[p] == 0 for p in kept):
            return None
        points = integer_points(self.vertices)
        values = [Q(n, p[-1]) for n, p in zip(numerators, points)]
        masks = [0] * len(values)
        facet_masks = [0] * len(self.faces)
        for j, face in enumerate(self.faces):
            for p in face:
                masks[p] |= 1 << j
                facet_masks[j] |= 1 << p
        vertices = {self.vertices[p]: masks[p] | ((values[p] == 0) << len(self.faces)) for p in kept}
        for i in kept:
            if values[i] == 0:
                continue
            for j, below in enumerate(values):
                if below >= 0:
                    continue
                pair = (1 << i) | (1 << j)
                common = shared = masks[i] & masks[j]
                if bin(shared).count('1') < self.dimension - 1:
                    continue
                face = (1 << len(values)) - 1
                while shared:
                    bit = shared & -shared
                    face &= facet_masks[bit.bit_length() - 1]
                    shared ^= bit
                if face != pair:
                    continue
                point = combine(-below/(values[i]-below), self.vertices[i],
                                values[i]/(values[i]-below), self.vertices[j])
                vertices[point] = common | (1 << len(self.faces))
        points = tuple(sorted(vertices))
        faces = tuple(sorted({tuple(p for p, point in enumerate(points) if vertices[point] & (1 << j))
                             for j in range(len(self.faces) + 1)} - {()}))
        assert len(points) >= self.dimension + 1
        return Poly(self.dimension, points, faces)

    def interior(self, randomizer, trial):
        if trial:
            chosen = [randomizer.randrange(1, 31) for _ in self.vertices]
        else:
            chosen = (1,) * len(self.vertices)
        total = sum(chosen)
        return tuple(dot(values, chosen) / total for values in zip(*self.vertices))


@lru_cache(maxsize=8192)
def integer_points(points):
    result = []
    for point in points:
        den = lcm(*(p.denominator for p in point))
        result.append(tuple(int(p*den) for p in point) + (den,))
    return tuple(result)


def chamber():
    vertices = tuple(tuple(Q(j < k) for j in range(5)) for k in range(6))
    poly = Poly(5, vertices, tuple(tuple(j for j in range(6) if j != k) for k in range(6)))
    f = tuple(form[:6] for form in IDENTITY)
    guards = set()
    for _ in range(5):
        for i in range(5):
            guards.add(normalized(combine(1, combine(1, f[i], 1, f[4-i]), -1, CENTER)))
            if i < 4:
                guards.add(normalized(combine(-1, combine(1, f[4-i], 1, f[i+1]), 1, CENTER)))
        f = (f[0], combine(1, f[0], -1, f[4]), combine(1, f[1], -1, f[4]),
             combine(1, f[1], -1, f[3]), combine(1, f[2], -1, f[3]))
    assert f == tuple(form[:6] for form in IDENTITY)
    for guard in sorted(guards):
        poly = poly.clip(guard)
    assert len(poly.vertices) == 77
    energy = [sum((a-b-Q(1, 6))**2 for a, b in zip((1, *x), (*x, 0))) for x in poly.vertices]
    assert max(energy) == Q(1, 30) and energy.count(Q(1, 30)) == 1
    return poly


def turn(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def hull(points):
    points = sorted(set(points))
    def boundary(walk):
        chain = []
        for point in walk:
            while len(chain) >= 2 and turn(*chain[-2:], point) <= 0:
                chain.pop()
            chain.append(point)
        return chain[:-1]
    result = boundary(points) + boundary(reversed(points))
    return tuple((*p, Q(0), Q(0), Q(0)) for p in result)


def plane(poly, chart, sheet):
    direction_one, direction_two = sheet.basis
    i, j, determinant = next((i, j, direction_one[i]*direction_two[j] - direction_one[j]*direction_two[i])
                             for i, j in combinations(range(5), 2)
                             if direction_one[i]*direction_two[j] != direction_one[j]*direction_two[i])
    one = tuple(Q(a*direction_two[j] - b*direction_one[j], determinant)
                for a, b in zip(direction_one, direction_two))
    two = tuple(Q(b*direction_one[i] - a*direction_two[i], determinant)
                for a, b in zip(direction_one, direction_two))
    new_chart = []
    for k in range(5):
        constant = combine(1, chart[k], -1, combine(one[k], chart[i], two[k], chart[j]))
        assert not any(constant[:5])
        new_chart.append((one[k], two[k], 0, 0, 0, constant[5]))
    points = hull(tuple((leading(chart[i], p), leading(chart[j], p)) for p in poly.vertices))
    assert len(points) >= 3
    faces = tuple(tuple(sorted((i, (i+1) % len(points)))) for i in range(len(points)))
    return Poly(2, points, faces), tuple(new_chart)


class Reader:
    def __init__(self, exe):
        self.proc = subprocess.Popen([exe], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     text=True)

    def close(self):
        self.proc.stdin.close()
        assert self.proc.wait() == 0

    def __call__(self, chart, point, sheet, budget=None):
        den = lcm(*(Q(v).denominator for row in chart for v in row))
        qden = lcm(*(Q(v).denominator for v in point))
        args = (HORIZON, SHORT_CAP if budget is None else budget, len(sheet.basis), den,
                *(int(den * x) for f in chart for x in f),
                *(int(qden * x) for x in point), qden,
                *sheet.point, sheet.base, 1, *sheet.slope, sheet.stride,
                *(x for c in sheet.basis for x in c))
        self.proc.stdin.write(' '.join(map(str, args)) + '\n')
        self.proc.stdin.flush()
        result = self.proc.stdout.readline().split()
        assert result, ('reader aborted', self.proc.poll(), chart, point, sheet)
        kind, *status = result
        if kind == 'stage':
            drift = self.row()
            nums = tuple(self.row() for _ in range(5))
            den, = self.row()
            chart = tuple(tuple(Q(a, den) for a in f) for f in nums)
            *point, base, unit = self.row()
            *slope, stride = self.row()
            assert unit == 1
            rank, = self.row()
            sheet = Sheet(base, stride, tuple(point), tuple(slope), tuple(self.row() for _ in range(rank))).reduced()
            races, = self.row()
            choices = {}
            for _ in range(races):
                *num, den = self.row()
                guard = self.row()
                step, = self.row()
                choices[tuple(Q(a, den) for a in num)] = (guard, step)
            read = (*map(int, status), chart, sheet, drift, choices)
        else:
            read = None
        count, = self.row()
        guards = tuple(self.row() for _ in range(count))
        return kind, read, guards

    def row(self):
        row = self.proc.stdout.readline()
        assert row, 'reader aborted mid-response'
        return tuple(map(int, row.split()))


def partition(poly, guards):
    for guard in guards:
        part = poly.clip(guard)
        assert part is not None
        if part is not poly:
            remainder = poly.clip(tuple(-a for a in guard))
            assert remainder is not None
            yield remainder, False
        poly = part
    yield poly, True


@lru_cache(maxsize=16384)
def winning_cell(domain, race, winner):
    for other in race:
        if other != winner:
            domain = domain.clip(normalized(combine(1, other, -1, winner)))
            if domain is None:
                break
    return domain


def central(chart):
    return all(f == (0, 0, 0, 0, 0, c) for f, c in zip(chart, CENTRE))


class Ranks:
    def __init__(self, reader, root, final_rank):
        self.reader = reader
        self.root = root
        self.final_rank = final_rank
        self.stats = Counter()
        self.random = random.Random(498 + root)
        self.seen = set()
        self.lines = Lines(reader, root)

    def run(self):
        root = self.root
        initial = Sheet(root, 6, ZERO, ZERO, tuple(f[:5] for f in IDENTITY))
        queue = [(chamber(), tuple(f[:6] for f in IDENTITY), initial, (5,), False)]
        while queue:
            domain, chart, sheet, history, compressed = queue.pop()
            rank = history[-1]
            if rank == 2 and not compressed:
                domain, chart = plane(domain, chart, sheet)
                compressed = True
            key = (domain, chart, sheet, history)
            if key in self.seen:
                continue
            self.seen.add(key)
            for trial in range(35):
                point = domain.interior(self.random, trial)
                kind, read, guards = self.reader(chart, point, sheet)
                if kind != 'wall':
                    break
            assert kind == 'stage', (kind, history, domain, chart, point, sheet)
            for piece, match in partition(domain, guards):
                if match:
                    break
                queue.append((piece, chart, sheet, history, compressed))
            self.stats[('reads', *history)] += 1
            shorts, length, losses, shortchart, shortsheet, drift, choices = read
            self.stats[('max shorts', *history)] = max(self.stats[('max shorts', *history)], shorts)
            self.stats[('period', rank, length, bool(losses))] += 1
            for time, (guard, decrement) in choices.items():
                active = winning_cell(piece, tuple(choices), time)
                if active is None:
                    continue
                self.stats[('exits', *history)] += 1
                out = exit_chart(shortchart, drift, time)
                children = {child for _, _, child in shortsheet.exits(guard, decrement, drift)}
                (next_rank,) = {len(child.basis) for child in children}
                if rank >= 3:
                    assert next_rank == rank - 1
                else:
                    assert len(history) <= 6 and (len(history) < 6 or next_rank == 1), history
                path = (*history, next_rank)
                if next_rank < self.final_rank:
                    if next_rank == 1:
                        for child in children:
                            self.lines.add_sheet(active, out, child)
                            self.stats[('arrivals', *path)] += 1
                    continue
                for child in children:
                    queue.append((active, out, child, path, False))
            reads = len(self.seen)
            if reads % 2500 == 0:
                print('ranks progress', root, reads, len(queue), flush=True)
        print('ranks closed', root, sorted(self.stats.items()), flush=True)
        if self.final_rank == 2:
            self.lines.run()


def normal_form(sheet, chart):
    (column,) = sheet.basis
    size = gcd(*column)
    direction = tuple(x // size for x in column)
    axis = next(i for i, x in enumerate(direction) if x)
    if direction[axis] < 0:
        direction = tuple(-a for a in direction)
    coefficients = []
    n = 0
    for x in direction:
        n, a, b = bezout(n, x)
        coefficients = [a*c for c in coefficients] + [b]
    assert n == 1 and dot(coefficients, direction) == 1
    def normal(v):
        return combine(1, v, -dot(coefficients, v), direction)
    slope = tuple(Q(a, sheet.stride) for a in normal(sheet.slope))
    stride = lcm(6, *(a.denominator for a in slope))
    assert sheet.stride % stride == 0
    base = sheet.base % stride
    offset = combine(1, normal(sheet.point), base-sheet.base, slope)
    out_slope = tuple(stride*a for a in slope)
    assert all(a.denominator == 1 for a in (*offset, *out_slope))
    rebuilt = Sheet(base, stride, tuple(map(int, offset)), tuple(map(int, out_slope)), (direction,))
    relative = tuple(combine(1, row, -centre, CENTER) for row, centre in zip(chart, CENTRE))
    pivot = relative[axis]
    ratio = tuple(Q(x, direction[axis]) for x in direction)
    assert all(combine(1, row, -a, pivot) == CENTER_ZERO for row, a in zip(relative, ratio))
    new_chart = tuple((a, 0, 0, 0, 0, c) for a, c in zip(ratio, CENTRE))
    return rebuilt, new_chart, pivot


CENTER_ZERO = (0,) * 6


def ranges(lo, hi, covered):
    for a, b in covered:
        if b <= lo:
            continue
        if a >= hi:
            break
        if a > lo:
            yield lo, a
        lo = max(b, lo)
        if lo >= hi:
            return
    if lo < hi:
        yield lo, hi


def union(old, extra):
    combined = []
    for a, b in sorted([*old, extra]):
        if combined and a <= combined[-1][1]:
            combined[-1] = combined[-1][0], max(b, combined[-1][1])
        else:
            combined.append((a, b))
    return combined


def interval_test(lo, hi, guard):
    assert not any(guard[1:5])
    a, b = guard[0], guard[5]
    if a == 0:
        return (lo, hi) if b >= 0 else None
    x = Q(-b, a)
    if a > 0:
        lo = max(lo, x)
    else:
        hi = min(hi, x)
    return (lo, hi) if lo < hi else None


class Lines:
    def __init__(self, reader, root):
        self.reader = reader
        self.root = root
        self.cover = defaultdict(list)
        self.edges = defaultdict(set)
        self.work = []
        self.minimum = None
        self.stats = Counter()

    def add(self, sheet, out_chart, lo, hi):
        sheet, chart, pivot = normal_form(sheet, out_chart)
        lo, hi = sorted((pivot[0]*lo+pivot[5], pivot[0]*hi+pivot[5]))
        return self.insert(sheet, chart, lo, hi)

    def add_sheet(self, poly, out_chart, sheet):
        sheet, chart, pivot = normal_form(sheet, out_chart)
        values = [leading(pivot, x) for x in poly.vertices]
        self.insert(sheet, chart, min(values), max(values))

    def insert(self, sheet, chart, lo, hi):
        assert lo < hi
        positive = Q(1, 7*radius(tuple(leading(f, (1, 0, 0, 0, 0)) for f in chart)))
        negative = -Q(1, 7*radius(tuple(leading(f, (-1, 0, 0, 0, 0)) for f in chart)))
        fragments = [(0, lo, min(hi, negative)), (0, max(lo, positive), hi)]
        if max(lo, negative) < min(hi, 0):
            fragments.append((-1, negative, Q(0)))
        if max(lo, 0) < min(hi, positive):
            fragments.append((1, Q(0), positive))
        keys = set()
        for code, a, b in fragments:
            if a >= b:
                continue
            key = (sheet, code)
            keys.add(key)
            previous = self.cover[key]
            for interior in ranges(a, b, previous):
                self.work.append((key, chart, interior))
            self.cover[key] = union(previous, (a, b))
        return keys

    def run(self):
        print('line ingress', self.root, len(self.cover), len(self.work), flush=True)
        while self.work:
            key, chart, (lo, hi) = self.work.pop()
            sheet, code = key
            for n in range(2, 35):
                q = lo+(hi-lo)/n
                kind, read, guards = self.reader(chart, (q, 0, 0, 0, 0), sheet, budget=250)
                if kind != 'wall':
                    break
            assert kind == 'stage', (kind, key, lo, hi, q)
            for guard in guards:
                rejected = interval_test(lo, hi, tuple(-a for a in guard))
                interval = interval_test(lo, hi, guard)
                assert interval
                if rejected and interval != (lo, hi):
                    self.work.append((key, chart, rejected))
                lo, hi = interval
            shorts, length, failures, short_chart, short_sheet, drift, choices = read
            assert not failures
            self.stats[('period', length)] += 1
            self.stats['reads'] += 1
            self.stats['shorts'] = max(self.stats['shorts'], shorts)
            for time, (guard, step) in choices.items():
                active = lo, hi
                for other in choices:
                    if other != time:
                        active = interval_test(*active, combine(1, other, -1, time))
                        if active is None:
                            break
                if active is None:
                    continue
                if code:
                    assert not time[5]
                    at = (sum(active)/2, 0, 0, 0, 0)
                    r = radius(tuple(leading(f, at) for f in chart))
                    ratio = length * leading(time, at) / r
                    assert ratio > 0
                    if self.minimum is None or ratio < self.minimum:
                        self.minimum = ratio
                target = exit_chart(short_chart, drift, time)
                if central(target):
                    continue
                for child in {child for _, _, child in short_sheet.exits(guard, step, drift)}:
                    for child_key in self.add(child, target, *active):
                        assert not code or child_key[1]
                        self.edges[key].add(child_key)
            if self.stats['reads'] % 5000 == 0:
                print('line progress', self.root, self.stats['reads'], len(self.cover), len(self.work), flush=True)
        remaining = {key: 0 for key in self.cover if not key[1]}
        for key in remaining:
            for target in self.edges[key]:
                if target in remaining:
                    remaining[target] += 1
        queue = deque(key for key, degree in remaining.items() if not degree)
        visited = 0
        while queue:
            key = queue.popleft()
            visited += 1
            for target in self.edges[key]:
                if target in remaining:
                    remaining[target] -= 1
                    if not remaining[target]:
                        queue.append(target)
        count = len(remaining) - visited
        print('line closed', self.root, len(self.cover), sum(map(len, self.edges.values())),
              'outer remainder', count, 'minimum', self.minimum, self.stats, flush=True)
        assert count == 0 and self.minimum == Q(5, 48)


def audit(executable, root, rank):
    reader = Reader(executable)
    try:
        Ranks(reader, root, rank).run()
    finally:
        reader.close()


def run(args):
    if args.executable:
        with ProcessPoolExecutor(max_workers=len(args.roots)) as pool:
            futures = [pool.submit(audit, args.executable, root, args.rank) for root in args.roots]
            for future in as_completed(futures):
                future.result()
        return
    with tempfile.TemporaryDirectory(prefix='kaprekar-task145-') as temporary, open(ROOT / '.task-monitor/compute.lock', 'rb') as lock:
        folder = Path(temporary)
        for source in ('reader.cpp', 'cover_exact.py', 'exact_model.py'):
            shutil.copy2(Path(__file__).with_name(source), folder)
        fcntl.flock(lock, fcntl.LOCK_EX)
        reader = folder / 'reader'
        subprocess.run(['clang++', '-std=c++17', '-O3', '-ftrapv', str(folder/'reader.cpp'), '-o', str(reader)], check=True)
        rank = str(args.rank)
        subprocess.run([sys.executable, '-u', '-B', str(folder/'cover_exact.py'),
                        '--executable', str(reader), '--rank', rank,
                        '--roots', *map(str, args.roots)], cwd=folder, check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--roots', type=int, nargs='+', default=list(range(6)))
    parser.add_argument('--rank', type=int, choices=(2, 3), default=3)
    parser.add_argument('--executable')
    run(parser.parse_args())
