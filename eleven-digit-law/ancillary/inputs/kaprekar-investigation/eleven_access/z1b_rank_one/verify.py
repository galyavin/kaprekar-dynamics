"""Replay rank-one ingress and verify a physical overcover in one dimension."""
import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction
from functools import lru_cache
from io import StringIO
from itertools import chain
from math import gcd, lcm
from pathlib import Path
import fcntl
import pickle
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = ACCESS.parents[1]
sys.path.insert(0, str(ACCESS / 'z1_stages'))
import cover

RANK_ONE = 1
UNSAFE, NEGATIVE, POSITIVE = 0, -1, 1


@lru_cache(maxsize=512)
def normal_projection(column):
    """Return the primitive tangent and an integral retraction onto it."""
    common = gcd(*column)
    u = tuple(x // common for x in column)
    common = 0
    functional = []
    for a in u:
        common, old, new = cover.extended_gcd(common, a)
        functional = [old * x for x in functional] + [new]
    assert common == 1 and sum(a * b for a, b in zip(u, functional)) == 1
    axis = next(j for j, a in enumerate(u) if a)
    assert u[axis] > 0
    return u, tuple(functional), axis


def project_normal(vector, tangent, functional):
    coordinate = sum(a * b for a, b in zip(vector, functional))
    return tuple(x - coordinate * u for x, u in zip(vector, tangent))


@lru_cache(maxsize=4096)
def canonical(sheet):
    """Enlarge to all physical lifts having the same exact affine normal mark."""
    base, period, rank, *affine = sheet
    slope, point, column = affine[:5], affine[5:10], tuple(affine[10:])
    beta_point = tuple(((j - 5) * base) % 6 for j in range(5))
    assert all((x - b) % 6 == 0 for x, b in zip(point, beta_point))
    assert period % 6 == 0 and all(x % 6 == 0 for x in slope)
    if rank == RANK_ONE:
        tangent, functional, axis = normal_projection(column)
    else:
        assert rank == 0 and not column
        tangent = functional = (0,) * 5
        axis = None
    normal_slope = project_normal(slope, tangent, functional)
    a = tuple(Fraction(x, period) for x in normal_slope)
    normal_period = 6 * lcm(*(x.denominator for x in a))
    assert period % normal_period == 0
    normal_base = base % normal_period
    shift = base - normal_base
    relative = project_normal(tuple(x - b for x, b in zip(point, beta_point)), tangent, functional)
    outpoint = tuple(b + n - shift * x for b, n, x in zip(beta_point, relative, a))
    assert all(x.denominator == 1 and (x - b) % 6 == 0 for x, b in zip(outpoint, beta_point))
    new_sheet = (normal_base, normal_period, rank, *(int(normal_period * x) for x in a),
                 *map(int, outpoint), *(6 * u for u in tangent if rank))
    if rank:
        direction = tuple(Fraction(x, tangent[axis]) for x in tangent)
        intercept = tuple(x - a[axis] * u for x, u in zip(a, direction))
    else:
        intercept, direction = a, (0,) * 5
    return new_sheet, intercept, direction, axis


@lru_cache(maxsize=64)
def ingress_interval(poly, chart, axis):
    row = chart[axis * 7:(axis + 1) * 7]
    values = [Fraction(sum(a * x for a, x in zip(row[:5], cover.xchart(g))) + row[5], row[6])
              for g in poly.vertices]
    return min(values), max(values)


@lru_cache(maxsize=512)
def check_ingress_chart(chart, intercept, direction, axis):
    pivot = tuple(Fraction(x, chart[axis * 7 + 6]) for x in chart[axis * 7:axis * 7 + 6])
    for k, (a, v) in enumerate(zip(intercept, direction)):
        row = chart[k * 7:(k + 1) * 7]
        assert tuple(Fraction(x, row[6]) for x in row[:6]) == tuple(v * x for x in pivot[:5]) + (a + v * pivot[5],)


def split_regions(intercept, direction, interval):
    """Only homogeneous permanent half-rays are widened."""
    lo, hi = interval
    if any(intercept) or not any(direction):
        yield UNSAFE, interval
        return
    path = (0, *direction, 0)
    derivative = tuple(b - a for a, b in zip(path, path[1:]))
    inner_lo = Fraction(-1, 7 * max(-a for a in derivative))
    inner_hi = Fraction(1, 7 * max(derivative))
    for sign, left, right in ((NEGATIVE, inner_lo, Fraction(0)), (POSITIVE, Fraction(0), inner_hi)):
        if max(lo, left) < min(hi, right):
            yield sign, (left, right)
    if lo < inner_lo:
        yield UNSAFE, (lo, min(hi, inner_lo))
    if hi > inner_hi:
        yield UNSAFE, (max(lo, inner_hi), hi)


class Sections:
    def __init__(self):
        self.coverage = {}
        self.pending = deque()
        self.edges = {}
        self.boundary_edges = {}
        self.stats = Counter()
        self.contraction = None

    def add(self, sheet, interval):
        sheet, intercept, direction, axis = canonical(sheet)
        for region, limits in split_regions(intercept, direction, interval):
            key = sheet, region
            previous = self.coverage.get(key, ())
            uncovered = tuple(cover.remaining_intervals(limits, previous))
            if uncovered:
                self.coverage[key] = cover.add_interval(limits, previous)
                self.pending.extend((key, segment) for segment in uncovered)
            yield key

    def arrival(self, poly, chart, sheet, ranks):
        self.stats['arrivals', ranks] += 1
        _, intercept, direction, axis = canonical(sheet)
        check_ingress_chart(chart, intercept, direction, axis)
        interval = ingress_interval(poly, chart, axis)
        self.stats['ingress', 'nonhomogeneous' if any(intercept) else 'homogeneous'] += 1
        self.stats['ingress', 'all'] += 1
        list(self.add(sheet, interval))

    def child(self, source, sheet, interval):
        for target in self.add(sheet, interval):
            if source[1] != UNSAFE:
                assert target[1] != UNSAFE, (source, target, interval)
            self.edges.setdefault(source, set()).add(target)
            yield target


class IngressQueue(deque):
    def append(self, item):
        poly, chart, sheet, ranks, _, _ = item
        if sheet[2] == RANK_ONE:
            self.sections.arrival(poly, chart, sheet, ranks)
        else:
            super().append(item)


def path_point(intercept, direction, t):
    values = tuple(a + t * v for a, v in zip(intercept, direction))
    den = lcm(*(p.denominator for p in values))
    point = (den, *(int(den * p) for p in values))
    assert max(map(abs, point)) < 2**63
    return point


def clip(interval, row):
    """Clip by a*g+b >= 0 in the synthetic first gap g."""
    if interval is None:
        return None
    lo, hi = interval
    a, b = row
    if a == 0:
        return interval if b >= 0 else None
    cut = -Fraction(b, a)
    if a > 0:
        lo = max(lo, cut)
    else:
        hi = min(hi, cut)
    return (lo, hi) if lo < hi else None


def restrict(interval, rows):
    for row in rows:
        interval = clip(interval, row)
    return interval


def gap_row(form):
    """Use x0=1-6g; the other source coordinates have zero coefficient."""
    a, *rest, b = form
    assert not any(rest)
    return -6 * a, a + b


def output_interval(out, axis, interval):
    if axis is None:
        return Fraction(0), Fraction(1)
    row = out[axis * 7:(axis + 1) * 7]
    a, b = gap_row(row[:6])
    values = [(a * t + b) / row[6] for t in interval]
    return min(values), max(values)


def record_boundary(sections, source, target, bounds, winning, duration, out, axis):
    """Retain zero-duration boundary transitions between nonpermanent sections."""
    if source[1] != UNSAFE or target[1] != UNSAFE or axis is None:
        return
    slope, intercept = duration
    if not slope:
        return
    zero = -intercept / slope
    if zero not in winning:
        return
    left, right = bounds
    entry = left + zero * (right - left)
    side = 1 if zero == winning[0] else -1
    row = out[axis * 7:(axis + 1) * 7]
    a, b = gap_row(row[:6])
    assert a
    exit_at_zero = (a * zero + b) / row[6]
    exit_side = side if a > 0 else -side
    start = source, entry, side
    end = target, exit_at_zero, exit_side
    sections.boundary_edges.setdefault(start, set()).add(end)


def certify_read(reader, sections, key, bounds, args, unseen):
    sheet, region = key
    sheet, intercept, direction, axis = canonical(sheet)
    lo, hi = bounds
    chart = cover.new_line_chart(intercept, direction, bounds)
    for fraction in chain((Fraction(1, 2),), (Fraction(1, n) for n in range(3, 28))):
        point = path_point(intercept, direction, lo + fraction * (hi - lo))
        kind, details, exits, guards = cover.reply(reader, point, chart, sheet, args.shorts)
        if kind != 'ambiguous':
            break
    inside = (Fraction(0), Fraction(1))
    for row in sorted(guards):
        assert len(set(row[1:])) == 1
        form = row[0] - row[1], row[1]
        rest = clip(inside, tuple(-x for x in form))
        if rest:
            if lo == hi:
                raise AssertionError('a proper wall on a fixed point')
            unseen.append((lo + rest[0] * (hi - lo), lo + rest[1] * (hi - lo)))
        inside = clip(inside, form)
        assert inside is not None
    if kind not in {'terminal', 'stage'}:
        raise RuntimeError(f'{kind=} {details=} {point=} {sheet=} {bounds=} {region=}')
    sections.stats['read', 'permanent ray' if region else 'other'] += 1
    sections.stats['read', 'nonhomogeneous' if any(intercept) else 'homogeneous'] += 1
    if kind == 'terminal':
        sections.stats['terminal', 'short'] += 1
        return
    short, _, period, failures, rank = details
    sections.stats['max shorts', rank] = max(sections.stats['max shorts', rank], short)
    sections.stats['stage', rank, period, bool(failures)] += 1
    distances = [gap_row(exit[0]) for exit in exits]
    for (time, _, out, children), (slope, offset) in zip(exits, distances):
        winning = restrict(inside, [(a - slope, b - offset) for a, b in distances if (a, b) != (slope, offset)])
        if not winning:
            continue
        if region:
            assert not failures and slope * (-lo / (hi - lo)) + offset == 0
            g = sum(winning) / 2
            t = lo + g * (hi - lo)
            p = (0, *(t * a for a in direction), 0)
            r = max(b - a for a, b in zip(p, p[1:]))
            ratio = period * (slope * g + offset) / r
            assert ratio > 0
            if sections.contraction is None or ratio < sections.contraction:
                sections.contraction = ratio
        if not children or all(out[i] == 0 for i in range(35) if i % 7 != 6):
            sections.stats['terminal', 'exit'] += 1
            continue
        for child in {x[0] for x in children}:
            new_sheet, a, v, axis = canonical(child)
            interval = output_interval(out, axis, winning)
            if axis is not None:
                check_ingress_chart(out, a, v, axis)
            else:
                assert all(tuple(Fraction(x, out[i + 6]) for x in out[i:i + 6]) == (0,) * 5 + (entry,)
                           for i, entry in zip(range(0, 35, 7), a))
            for target in sections.child(key, new_sheet, interval):
                record_boundary(sections, key, target, bounds, winning, (slope, offset), out, axis)


def cycle_remainder(nodes, edges):
    incoming = dict.fromkeys(nodes, 0)
    for source in nodes:
        for target in edges.get(source, ()):
            if target in nodes:
                incoming[target] += 1
    roots = deque(key for key, count in incoming.items() if count == 0)
    visited = 0
    while roots:
        source = roots.popleft()
        visited += 1
        for target in edges.get(source, ()):
            if target in nodes:
                incoming[target] -= 1
                if incoming[target] == 0:
                    roots.append(target)
    return len(nodes) - visited


def walk(reader, sections, args):
    seen = 0
    print('rank-one ingress', len(sections.coverage), 'sections', len(sections.pending), 'intervals',
          sorted(sections.stats.items()), flush=True)
    while sections.pending:
        key, bounds = sections.pending.pop()
        unseen = []
        certify_read(reader, sections, key, bounds, args, unseen)
        sections.pending.extend((key, rest) for rest in unseen)
        seen += 1
        if seen % 1000 == 0:
            print('rank-one progress', seen, 'sections', len(sections.coverage),
                  'pending', len(sections.pending), sorted(sections.stats.items()), flush=True)
        if seen >= args.line_budget and sections.pending:
            print('rank-one budget', seen, 'sections', len(sections.coverage),
                  'pending', len(sections.pending), sorted(sections.stats.items()), flush=True)
            return False
    outside = {key for key in sections.coverage if key[1] == UNSAFE}
    cyclic = cycle_remainder(outside, sections.edges)
    boundary_cyclic = cycle_remainder(set(sections.boundary_edges), sections.boundary_edges)
    print('rank-one closed', seen, 'sections', len(sections.coverage), 'edges',
          sum(map(len, sections.edges.values())), 'outside topological remainder', cyclic,
          'zero-duration boundary topological remainder', boundary_cyclic,
          'min elapsed/radius on rays', sections.contraction, sorted(sections.stats.items()), flush=True)
    return boundary_cyclic == 0


def reader_process(executable):
    return subprocess.Popen([str(executable)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)


def run_root(executable, continuation, root, args):
    checkpoint = args.state_dir / f'root-{root}.ingress' if args.state_dir else None
    if args.resume and checkpoint and checkpoint.exists():
        with checkpoint.open('rb') as f:
            sections = pickle.load(f)
        print('resumed validated rank-one ingress for root', root, flush=True)
    else:
        sections = capture_ingress(executable, root, args)
        if checkpoint:
            with checkpoint.with_suffix('.tmp').open('wb') as f:
                pickle.dump(sections, f)
            checkpoint.with_suffix('.tmp').replace(checkpoint)
    if args.ingress_only:
        return 0
    with reader_process(continuation) as reader:
        try:
            ok = walk(reader, sections, args)
        finally:
            reader.stdin.close()
        assert reader.wait() == 0
    return 0 if ok else 3


def capture_ingress(executable, root, args):
    sections = Sections()
    IngressQueue.sections = sections
    original_queue = cover.deque
    cover.deque = IngressQueue
    with reader_process(executable) as reader:
        config = argparse.Namespace(offsets=[root], rank=1, compress_planes=True, saturate_planes=False,
                                    compress_lines=False, shorts=args.shorts, limit=args.ingress_budget,
                                    keep_limit=2)
        try:
            report = StringIO()
            with redirect_stdout(report):
                result = cover.search(reader, config)
        finally:
            cover.deque = original_queue
            reader.stdin.close()
        assert reader.wait() == 0
    output = report.getvalue().splitlines()
    print(output[-1], flush=True)
    if result:
        raise RuntimeError(f'rank-one premise failed to close for root {root}: {result}')
    expected = next(line for line in (ACCESS / 'z1_stages/rank-two-inherited.out').read_text().splitlines()
                    if line.startswith(f'closed [{root}] '))
    assert output[-1] == expected, (output[-1], expected)
    assert {history for (label, history) in sections.stats if label == 'arrivals'} == {
        (5, 4, 3, 2, 1), (5, 4, 3, 2, 2, 1), (5, 4, 3, 2, 2, 2, 1)}
    print('rank-one ingress replay matches root', root, flush=True)
    return sections


def build(snapshot, paired_fast):
    snapshot.mkdir()
    for path in (ACCESS / 'z1_stages/chart.cpp', ACCESS / 'stage_graph/graph.cpp',
                 ACCESS / 'stage_graph/lottery.cpp'):
        shutil.copyfile(path, snapshot / path.name)
    cover.record_stage_probes(snapshot / 'graph.cpp', paired_fast)
    subprocess.run(['c++', '-O3', '-std=c++17', '-Wno-return-type', 'chart.cpp', '-o', 'chart'],
                   cwd=snapshot, check=True)
    return snapshot / 'chart'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roots', type=int, nargs='+', choices=range(6), default=list(range(6)))
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--shorts', type=int, default=500)
    parser.add_argument('--ingress-budget', type=int, default=50000)
    parser.add_argument('--line-budget', type=int, default=100000)
    parser.add_argument('--ingress-only', action='store_true')
    parser.add_argument('--state-dir', type=Path, help='private directory for validated ingress checkpoints')
    parser.add_argument('--resume', action='store_true', help='reuse checkpoints from this task in --state-dir')
    args = parser.parse_args()
    if args.resume and not args.state_dir:
        parser.error('a private --state-dir is required to resume')
    if args.state_dir:
        args.state_dir.mkdir(exist_ok=True)
    with (ROOT / '.task-monitor/compute.lock').open() as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix='task142-') as temp:
            root = Path(temp)
            executable = build(root / 'ingress', False)
            continuation = executable if args.ingress_only else build(root / 'continuation', True)
            with ProcessPoolExecutor(max_workers=min(args.jobs, len(args.roots))) as pool:
                results = [pool.submit(run_root, executable, continuation, root, args) for root in args.roots]
                return max(result.result() for result in results)


if __name__ == '__main__':
    raise SystemExit(main())
