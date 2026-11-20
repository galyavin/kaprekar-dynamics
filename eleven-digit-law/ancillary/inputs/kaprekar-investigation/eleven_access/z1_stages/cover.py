"""Partition the strict Pi source by exact stage guards and inherited lanes."""
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import argparse
from ast import literal_eval
import fcntl
from functools import lru_cache
from itertools import chain, combinations
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
sys.path.insert(0, str(ACCESS / "stage_graph"))
from pi_geometry import Poly, dirichlet_mass, within_interlevel  # noqa: E402
from geometry import integer_row, xchart  # noqa: E402


def flatten(rows):
    return tuple(chain.from_iterable(rows))


def centered_row(form):
    weights, offset = form[:5], form[5]
    return integer_row(offset + sum(w * (i - 6 * (j < i))
                                    for i, w in enumerate(weights, 1)) for j in range(6))


START = flatten(tuple(int(i == j) for i in range(5)) + (0, 1) for j in range(5))
ID = flatten(tuple(6 * int(i == j) for i in range(5)) for j in range(5))


def exact_point(poly, chart, rng, retry):
    weights = [rng.randrange(1, 51) if retry else 1 for _ in poly.vertices]
    gaps = tuple(sum(w * v[i] for w, v in zip(weights, poly.vertices)) / sum(weights)
                 for i in range(6))
    initial = xchart(gaps)
    point = [(sum(chart[i + k] * initial[k] for k in range(5)) + chart[i + 5]) / chart[i + 6]
             for i in range(0, 35, 7)]
    scale = lcm(*(p.denominator for p in point))
    return (scale, *(int(scale * p) for p in point))


def padded_form(row):
    denominator = lcm(*(a.denominator for a in row))
    return (*(int(denominator * a) for a in row), denominator)


def project_line(poly, chart, sheet):
    direction = sheet[13:18]
    axis = next(i for i, a in enumerate(direction) if a)
    basis = [Fraction(a, direction[axis]) for a in direction]
    forms = [tuple(Fraction(a, chart[i + 6]) for a in chart[i:i + 6])
             for i in range(0, 35, 7)]
    arg = forms[axis]
    intercept = [f[-1] - v * arg[-1] for f, v in zip(forms, basis)]
    assert all(tuple(a - v * b for a, b in zip(f[:-1], arg[:-1])) == (0,) * 5
               for f, v in zip(forms, basis))
    values = [sum(a * x for a, x in zip(arg[:-1], xchart(g))) + arg[-1] for g in poly.vertices]
    return tuple(intercept), tuple(basis), (min(values), max(values))


def saturate_line(sheet):
    slope, point, column = sheet[3:8], sheet[8:13], sheet[13:18]
    divisor = gcd(*column)
    assert divisor % 6 == 0
    primitive = tuple(6 * a // divisor for a in column)
    axis = next(i for i, a in enumerate(primitive) if a)
    assert primitive[axis] > 0

    def reduce_on_line(vector):
        shift = vector[axis] // primitive[axis]
        return tuple(a - shift * b for a, b in zip(vector, primitive))

    return (*sheet[:3], *reduce_on_line(slope), *reduce_on_line(point), *primitive)


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def convex_hull(points):
    ordered = sorted(set(points))

    def half(rows):
        boundary = []
        for p in rows:
            while len(boundary) >= 2 and cross(*boundary[-2:], p) <= 0:
                boundary.pop()
            boundary.append(p)
        return boundary[:-1]

    return half(ordered) + half(reversed(ordered))


def extended_gcd(a, b):
    if a == b == 0:
        return 0, 0, 0
    aa, bb = abs(a), abs(b)
    x, y, u, v = 1, 0, 0, 1
    while bb:
        q, remainder = divmod(aa, bb)
        aa, bb = bb, remainder
        x, y, u, v = u, v, x - q * u, y - q * v
    return aa, x if a >= 0 else -x, y if b >= 0 else -y


def combine(a, left, b, right):
    return tuple(a * x + b * y for x, y in zip(left, right))


def saturate_plane(sheet, i, j, ratio_i, ratio_j):
    one, two = (1, 0), (0, 1)
    for v, w in zip(ratio_i, ratio_j):
        modulus = lcm(v.denominator, w.denominator)
        if modulus == 1:
            continue
        a, b = int(modulus * v), int(modulus * w)
        first, second = a * one[0] + b * one[1], a * two[0] + b * two[1]
        common, x, y = extended_gcd(first, second)
        if not common:
            continue
        index = modulus // gcd(modulus, common)
        one, two = combine(index * x, one, index * y, two), combine(-second // common, one, first // common, two)

    def lift(pair):
        values = [a * pair[0] + b * pair[1] for a, b in zip(ratio_i, ratio_j)]
        assert all(a.denominator == 1 for a in values)
        return tuple(6 * int(a) for a in values)

    one, two = lift(one), lift(two)
    common, x, y = extended_gcd(one[i], two[i])
    one, two = combine(x, one, y, two), combine(-two[i] // common, one, one[i] // common, two)
    if two[j] < 0:
        two = tuple(-a for a in two)
    one = combine(1, one, - (one[j] // two[j]), two)

    def reduce(a):
        a = combine(1, a, - (a[i] // one[i]), one)
        return combine(1, a, - (a[j] // two[j]), two)

    assert not any(reduce(sheet[13:18])) and not any(reduce(sheet[18:23]))
    return (*sheet[:3], *reduce(sheet[3:8]), *reduce(sheet[8:13]), *one, *two)


def project_plane(poly, chart, sheet, simplex, saturate):
    a, b = sheet[13:18], sheet[18:23]
    pairs = ((i, j, a[i] * b[j] - a[j] * b[i]) for i, j in combinations(range(5), 2))
    i, j, determinant = next(p for p in pairs if p[2])
    forms = [tuple(Fraction(a, chart[k + 6]) for a in chart[k:k + 6]) for k in range(0, 35, 7)]
    ratio_i = [Fraction(a[k] * b[j] - b[k] * a[j], determinant) for k in range(5)]
    ratio_j = [Fraction(b[k] * a[i] - a[k] * b[i], determinant) for k in range(5)]
    projected = []
    for f, v, w in zip(forms, ratio_i, ratio_j):
        aff = tuple(x - v * y - w * z for x, y, z in zip(f, forms[i], forms[j]))
        assert aff[:-1] == (0,) * 5
        # t = -6+24g0 and u = -6+24g1, a triangle containing the physical square.
        projected.append(padded_form((-4 * v + 4 * w, -4 * w, 0, 0, 0, aff[-1] - 2 * (v + w))))

    points = []
    for g in poly.vertices:
        centered = xchart(g)
        points.append(tuple(f[-1] + sum(a * x for a, x in zip(f, centered)) for f in (forms[i], forms[j])))
    hull = convex_hull(points)
    assert len(hull) >= 3 and all(-6 < x < 6 for p in hull for x in p)
    triangle = ((18, -6), (-6, 18), *[(-6, -6)] * 4)
    cuts = sorted(integer_row(cross(v, w, p) for p in triangle) for v, w in zip(hull, hull[1:] + hull[:1]))
    domain = simplex.restrict(cuts)
    assert domain is not None
    if saturate:
        sheet = saturate_plane(sheet, i, j, ratio_i, ratio_j)
    return domain, flatten(projected), sheet, tuple(cuts)


def new_line_chart(intercept, direction, interval):
    low, high = interval
    width = high - low
    assert width > 0
    # The new free parameter ranges over the first gap of a standard simplex.
    # In its centered coordinates, that gap is (1-x0)/6.
    return flatten(padded_form((-v * width / 6, 0, 0, 0, 0, a + v * (low + width / 6)))
                   for a, v in zip(intercept, direction))


def extend_central_ray(intercept, direction, interval):
    lo, hi = interval
    if any(intercept) or lo < 0 < hi:
        return interval
    sign = 1 if lo >= 0 else -1
    coordinates = (0, *direction, 0)
    radius = max(sign * (b - a) for a, b in zip(coordinates, coordinates[1:]))
    end = 1 / (7 * radius)
    if abs(lo) <= end and abs(hi) <= end:
        return (0, end) if sign > 0 else (-end, 0)
    return interval


def central_line_domain(intercept, direction, interval):
    if any(intercept):
        return False
    path = (0, *direction, 0)
    return all(7 * t * (b - a) <= 1 for t in interval for a, b in zip(path, path[1:]))


def remaining_intervals(interval, previous):
    lo, hi = interval
    for low, high in previous:
        if high <= lo or hi <= low:
            continue
        if lo < low:
            yield lo, low
        lo = max(lo, high)
        if lo >= hi:
            return
    if lo < hi:
        yield lo, hi


def add_interval(interval, previous):
    union = []
    for low, high in sorted([interval, *previous]):
        if union and low <= union[-1][1]:
            union[-1] = (union[-1][0], max(high, union[-1][1]))
        else:
            union.append((low, high))
    return union


def reply(reader, point, chart, sheet, max_short):
    base, period, rank, *lat = sheet
    reader.stdin.write(" ".join(map(str, (*point, *chart, base, period, rank, max_short, *lat))) + "\n")
    reader.stdin.flush()
    row = reader.stdout.readline().split()
    if not row:
        raise RuntimeError(f"chart reader exited for {point=} {chart=} {sheet=}")
    kind, *details = row
    exits = []
    if kind == "stage":
        for _ in range(int(reader.stdout.readline())):
            *distance, denominator = map(int, reader.stdout.readline().split())
            time_row = tuple(Fraction(a, denominator) for a in distance)
            descriptor = tuple(map(int, reader.stdout.readline().split()))
            next_chart = tuple(map(int, reader.stdout.readline().split()))
            children = []
            for _ in range(int(reader.stdout.readline())):
                p, period, rem, index, r, *lat = map(int, reader.stdout.readline().split())
                children.append(((p, period, r, *lat), rem, index))
            exits.append((time_row, descriptor, next_chart, children))
    guards = {centered_row(tuple(map(int, reader.stdout.readline().split())))
              for _ in range(int(reader.stdout.readline()))}
    return kind, tuple(map(int, details)), exits, guards


def partition(poly, rows):
    # Work with the convex closure. A newly imposed proper equality has zero
    # original measure. The exact cuts partition the full-dimensional pieces.
    for row in sorted(rows):
        values = [sum(a * b for a, b in zip(row, vertex)) for vertex in integer_vertices(poly.vertices)]
        if min(values) >= 0:
            continue
        assert max(values) > 0
        other = poly.cut(tuple(-a for a in row))
        exact = poly.cut(row)
        assert exact is not None
        if other is not None and other is not poly:
            yield False, other
        poly = exact
    yield True, poly


@lru_cache(maxsize=2048)
def integer_vertices(vertices):
    points = []
    for vertex in vertices:
        scale = lcm(*(coordinate.denominator for coordinate in vertex))
        points.append(tuple(int(scale * coordinate) for coordinate in vertex))
    return points


@lru_cache(maxsize=2048)
def cached_partition(poly, rows):
    return tuple(partition(poly, rows))


@lru_cache(maxsize=2048)
def race_cell(poly, rows, winner):
    cuts = [centered_row(tuple(a - b for a, b in zip(other, winner)))
            for other in rows if other != winner]
    return poly.restrict(cuts)


def search(reader, args):
    desc, inter = within_interlevel()
    pi = Poly.root().restrict([*sorted(inter), *sorted(desc)])
    energy = [sum((g - Fraction(1, 6)) ** 2 for g in x) for x in pi.vertices]
    assert len(energy) == 77 and max(energy) == Fraction(1, 30) and energy.count(max(energy)) == 1
    pending = deque()
    for offset in args.offsets:
        shift = tuple(((i - 5) * offset) % 6 for i in range(5))
        pending.append((pi, START, (offset, 6, 5, *([0] * 5), *shift, *ID), (5,), (), False))
    rng = random.Random(137049)
    stats = Counter()
    seen = set()
    line_coverage = {}
    plane_coverage = {}
    simplex = Poly.root()
    wanted_prefix = getattr(args, "prefix", ())
    chart_reads = 0
    while pending:
        poly, chart, sheet, ranks, history, compressed = pending.pop()
        if args.compress_planes and sheet[2] == 2 and not compressed:
            poly, chart, sheet, cuts = project_plane(poly, chart, sheet, simplex, args.saturate_planes)
            plane_key = sheet, ranks, chart
            previous = plane_coverage.setdefault(plane_key, [])
            vertices = integer_vertices(poly.vertices)
            if any(all(sum(a * b for a, b in zip(row, vertex)) >= 0
                       for row in region for vertex in vertices) for region in previous):
                continue
            previous.append(cuts)
            compressed = True
        if args.compress_lines and sheet[2] == 1 and not compressed:
            intercept, direction, interval = project_line(poly, chart, sheet)
            if args.central_rays:
                interval = extend_central_ray(intercept, direction, interval)
            if args.saturate_lines:
                sheet = saturate_line(sheet)
                if args.central_rays and central_line_domain(intercept, direction, interval):
                    assert not any(sheet[3:8])
                    sheet = (sheet[0] % 6, 6, *sheet[2:])
            line_key = sheet, ranks, intercept, direction
            previous = line_coverage.get(line_key, ())
            for uncovered in remaining_intervals(interval, previous):
                pending.append((simplex, new_line_chart(intercept, direction, uncovered), sheet,
                                ranks, history, True))
            line_coverage[line_key] = add_interval(interval, previous)
            continue
        key = poly, chart, sheet, ranks
        if key in seen:
            continue
        seen.add(key)
        if sheet[2] == 0:
            raise RuntimeError(f"noncentral rank-zero section at {sheet=} {ranks=} {chart=} {history=}")
        for retry in range(20):
            point = exact_point(poly, chart, rng, retry)
            if max(map(abs, point)) >= 2**63:
                continue
            kind, details, exits, guards = reply(reader, point, chart, sheet, args.shorts)
            if kind != "ambiguous":
                break
        else:
            raise RuntimeError(f"no generic interior representative at {ranks=} {sheet=} {point=}")
        chart_reads += 1
        stats[kind, ranks] += 1
        inside = None
        for matches, part in cached_partition(poly, tuple(sorted(guards))):
            if matches:
                inside = part
            else:
                pending.append((part, chart, sheet, ranks, history, compressed))
        if kind != "stage" and kind != "terminal":
            print("unresolved", kind, details, "chart", chart, "sheet", sheet,
                  "ranks", ranks, "point", point, "history", history,
                  "compressed", compressed, flush=True)
            return 1
        if kind == "stage":
            n_short, search_length, period, failures, _ = details
            assert details[-1] == ranks[-1]
            stats["shortest-max", ranks] = max(stats["shortest-max", ranks], n_short)
            stats["physical period", ranks[-1], period, bool(failures)] += 1
            stats["probe max", ranks[-1]] = max(stats["probe max", ranks[-1]], search_length)
            stats["total short max", ranks] = max(stats["total short max", ranks],
                                                  n_short + sum(edge[4][0] for edge in history))
            rows = tuple(option[0] for option in exits)
            for winner, descriptor, out, children in exits:
                winning_cell = race_cell(inside, rows, winner)
                if winning_cell is None:
                    continue
                stats["exit", ranks] += 1
                by_state = {}
                if not children:
                    assert ranks[-1] == 1 and all(out[i] == 0 for i in range(35) if i % 7 != 6)
                    stats["terminal exit", (*ranks, 0)] += 1
                    continue
                for child, remainder, index in children:
                    by_state.setdefault(child, (remainder, index))
                for child, (remainder, index) in by_state.items():
                    rank = child[2]
                    assert ranks[-1] - rank in (0, 1)
                    next_ranks = (*ranks, rank)
                    edge = (sheet[:2], child[:2], remainder, index, (*details, *descriptor))
                    if len(history) < len(wanted_prefix):
                        wanted = wanted_prefix[len(history)]
                        # The probe horizon itself can differ after a chart reparametrization.
                        if edge[:4] != wanted[:4] or edge[4][:1] + edge[4][2:] != wanted[4][:1] + wanted[4][2:]:
                            continue
                    next_history = (*history, edge)
                    keeps = sum(a == b for a, b in zip(next_ranks, next_ranks[1:]))
                    desired = getattr(args, "stop_signature", "")
                    cap = getattr(args, "keep_limit", 2)
                    if rank == ranks[-1] >= 3 or keeps > cap or desired == "".join(map(str, next_ranks)):
                        synthetic = compressed or (args.compress_planes and 2 in ranks)
                        mass = None if synthetic else dirichlet_mass(winning_cell)
                        lane = Fraction(1)
                        for record in next_history:
                            lane /= record[3]
                        print("keep witness", next_ranks, "offset/period", child[:2],
                              "Pi mass", mass, "lane lower bound", lane,
                              "original direct mass lower bound", mass * lane if mass else None,
                              "vertices", winning_cell.vertices, "chart", out,
                              "sheet", child, "history", next_history, flush=True)
                        return 2
                    if all(out[i] == 0 for i in range(35) if i % 7 != 6):
                        stats["terminal exit", next_ranks] += 1
                    elif rank >= args.rank and (rank > 1 or ranks.count(2) - 1 >= getattr(args, "rank_one_min_rank_two_keeps", 0)):
                        # A new rank-one exit is reparametrized to its common line.
                        pending.append((winning_cell, out, child, next_ranks, next_history, False))
        if chart_reads % 1000 == 0:
            print("progress", args.offsets, chart_reads, "pending", len(pending),
                  "rank profile", sorted(stats.items()), flush=True)
        if chart_reads == args.limit:
            print("budget", args.offsets, chart_reads, "pending", len(pending),
                  "pending ranks", sorted(Counter(item[3] for item in pending).items()),
                  "rank profile", sorted(stats.items()), flush=True)
            return 3
    print("closed", args.offsets, chart_reads, "rank profile", sorted(stats.items()), flush=True)
    return 0


def run_reader(executable, args):
    with subprocess.Popen([str(executable)], stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE, text=True) as reader:
        result = search(reader, args)
        reader.stdin.close()
        if reader.wait():
            raise RuntimeError("chart reader failed")
    return result


def record_stage_probes(source, paired_fast=False):
    text = source.read_text()
    declaration = "Stage stage(const Context &ctx) {"
    phase = "Phase phase = primitive(f, ctx);"
    assert text.count(declaration) == 1 and text.count(phase) == 1
    text = text.replace(declaration, "#define STAGE_GUARD_PROBE\n"
                        "std::vector<Phase> *recorded_stage_phases = nullptr;\n" + declaration)
    text = text.replace(phase, phase + "\n        if (recorded_stage_phases) recorded_stage_phases->push_back(phase);")
    if paired_fast:
        sort = """    std::stable_sort(order.begin(), order.end(), [&](int i, int k) {
        return ctx.nonnegative(digits[i] - digits[k] - Form{{}, 0, i < k ? 0 : 6});
    });"""
        paired_sort = """    const auto before = [&](int i, int k) {
        return ctx.nonnegative(digits[i] - digits[k] - Form{{}, 0, i < k ? 0 : 6});
    };
    W old = 0, radius = 0;
    for (const W z : ctx.x.n) {
        radius = std::max(radius, z - old);
        old = z;
    }
    radius = std::max(radius, -old);
    if (r == 5 && 7 * radius < ctx.x.den) {
        order[0] = 5;
        for (int i = 0; i < 5; ++i) {
            const bool positive = before(i, 6 + i);
            order[2 * i + 1] = positive ? i : 6 + i;
            order[2 * i + 2] = positive ? 6 + i : i;
        }
    } else {
        std::stable_sort(order.begin(), order.end(), before);
    }"""
        assert text.count(sort) == 1
        text = text.replace(sort, paired_sort)
    source.write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offsets", nargs="+", type=int, default=list(range(6)))
    parser.add_argument("--rank", type=int, default=3)
    parser.add_argument("--shorts", type=int, default=500)
    parser.add_argument("--limit", type=int, default=10_000)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--compress-lines", action="store_true")
    parser.add_argument("--compress-planes", action="store_true")
    parser.add_argument("--saturate-lines", action="store_true")
    parser.add_argument("--saturate-planes", action="store_true")
    parser.add_argument("--central-rays", action="store_true")
    parser.add_argument("--stop-signature", default="")
    parser.add_argument("--paired-fast", action="store_true")
    parser.add_argument("--keep-limit", type=int, default=2)
    parser.add_argument("--rank-one-min-rank-two-keeps", type=int, default=0)
    parser.add_argument("--prefix", type=Path, help="restrict to this literal list of recorded edges")
    args = parser.parse_args()
    args.prefix = literal_eval(args.prefix.read_text()) if args.prefix else ()
    assert all(0 <= x < 6 for x in args.offsets)
    assert args.compress_lines or not (args.saturate_lines or args.central_rays)
    assert args.compress_planes or not args.saturate_planes
    with (ROOT / ".task-monitor" / "compute.lock").open() as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix="task137-cover-") as name:
            snapshot = Path(name)
            for source in [HERE / "chart.cpp", ACCESS / "stage_graph" / "graph.cpp",
                           ACCESS / "stage_graph" / "lottery.cpp"]:
                shutil.copyfile(source, snapshot / source.name)
            record_stage_probes(snapshot / "graph.cpp", args.paired_fast)
            subprocess.run(["c++", "-O3", "-std=c++17", "-Wno-return-type", "chart.cpp", "-o", "chart"],
                           cwd=snapshot, check=True)
            if args.jobs == 1:
                return run_reader(snapshot / "chart", args)
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                calls = [pool.submit(run_reader, snapshot / "chart",
                                     argparse.Namespace(**{**vars(args), "offsets": [offset]}))
                         for offset in args.offsets]
                return max(call.result() for call in calls)


if __name__ == "__main__":
    raise SystemExit(main())
