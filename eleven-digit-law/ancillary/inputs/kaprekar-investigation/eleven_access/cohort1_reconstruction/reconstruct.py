"""Reconstruct the cohort-1 orbit quotient from the Task 029 endpoint record."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cache
from itertools import combinations
import json
from pathlib import Path
from random import Random
import sys

SCALE = 90_000
BASE_MIN = 1_000_001
V = (5, 4, 3, 2, 1)
BOX = tuple((Q(a, 100_000), Q(b, 100_000)) for a, b in (
    (94040, 94130), (73770, 73850), (57940, 58015), (20, 80), (200, 300)
))
SECOND_DRIFT = (-20, -4, 0, 4, -4)
DETECTION_DRIFT = (-10, -11, -6, -10, -11)
ORDER_NAMES = {
    "5061728394a": "S",
    "5061782394a": "U",
    "5067128934a": "C",
    "5067182934a": "W",
}


@dataclass(frozen=True)
class Row:
    """Coefficients of (B, p, r, 1), multiplied by SCALE."""

    terms: tuple

    def __add__(self, other):
        if isinstance(other, int):
            a, b, c, d = self.terms
            return Row((a, b, c, d + other * SCALE))
        return Row(tuple(a + b for a, b in zip(self.terms, other.terms)))

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, integer):
        return Row(tuple(integer * a for a in self.terms))

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        assert all(a % divisor == 0 for a in self.terms)
        return Row(tuple(a // divisor for a in self.terms))

    def eval(self, base, p, r):
        return Q(sum(a * b for a, b in zip(self.terms, (base, p, r, 1))), SCALE)


B = Row((SCALE, 0, 0, 0))
P = Row((0, SCALE, 0, 0))
R = Row((0, 0, SCALE, 0))


def digit_rows(base, q):
    a, b, c, d, e = q
    return a, b, c, d, e - 1, base - 1, base - e - 1, base - d - 1, base - c - 1, base - b - 1, base - a


@cache
def indices(order):
    return tuple(int(index, 16) for index in order)


def subtract_in_order(base, q, order):
    digits = digit_rows(base, q)
    ix = indices(order)
    return tuple(digits[ix[i]] - digits[ix[-i - 1]] for i in range(5))


def step_many(base, q, orders):
    for order in orders:
        q = subtract_in_order(base, q, order)
    return q


def translate(q, count, delta):
    return tuple(row + count * d for row, d in zip(q, delta))


def difference(left, right):
    return tuple(a - b for a, b in zip(left, right))


def prefix_coordinates(base, a, b, c, h, g, z, r1):
    n = (h + 3) / 6
    m = (g - 1 + 5 * n + z) / 2
    N = (a + c - 8 * n - 10 * m - 3 * (base - 1) / 2 - r1) / 22
    return c - 3 * n - 2 * N, b - 4 * n - 2 * m - 4 * N, N


@cache
def chart_constants(z, r1):
    p, r, _ = prefix_coordinates(Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), z, r1)
    return p, r


@cache
def normalized_bounds(base, p, r):
    """Extrema of (base*B + p*p_coordinate + r*r_coordinate)/B at leading order."""
    p0, r0 = chart_constants(1, 0)
    p1, r1, _ = prefix_coordinates(Q(1), Q(0), Q(0), Q(0), Q(0), Q(0), 1, 0)
    low = high = base + p * (p1 - p0) + r * (r1 - r0)
    for axis, (lo, hi) in enumerate(BOX):
        chart = [Q(int(i == axis)) for i in range(5)]
        p1, r1, _ = prefix_coordinates(Q(0), *chart, 1, 0)
        coefficient = p * (p1 - p0) + r * (r1 - r0)
        lo, hi = sorted((coefficient * lo, coefficient * hi))
        low += lo
        high += hi
    return low, high


class BoxGuards:
    def __init__(self, z, r1):
        self.p0, self.r0 = chart_constants(z, r1)

    def bounds(self, row):
        a, b, c, d = row.terms
        low, high = normalized_bounds(a, b, c)
        intercept = b * self.p0 + c * self.r0 + d
        return low, high, intercept

    def nonnegative(self, row):
        low, _, intercept = self.bounds(row)
        return low >= 0 and low * BASE_MIN + intercept >= 0

    def paired_order(self, q):
        digits = digit_rows(B, q)
        order = [5]
        for a, b in zip(range(5), range(6, 11)):
            low, high, intercept = self.bounds(digits[a] - digits[b])
            if low >= 0 and low * BASE_MIN + intercept >= 0:
                order.extend((a, b))
            else:
                assert high <= 0 and high * BASE_MIN + intercept < 0, (q, a, b)
                order.extend((b, a))
        return "".join(format(d, "x") for d in order)

    def walk(self, q, orders):
        for order in orders:
            digits = digit_rows(B, q)
            ix = indices(order)
            assert all(self.nonnegative(digits[a] - digits[b]) for a, b in zip(ix, ix[1:])), (order, q)
            assert all(self.nonnegative(gap - 1) for gap in (
                B - q[0], *difference(q[:-1], q[1:]), q[-1]
            ))
            q = tuple(digits[ix[i]] - digits[ix[-i - 1]] for i in range(5))
        return q


@dataclass(frozen=True)
class Section:
    label: tuple
    word: tuple
    checkpoint: tuple
    displacement: tuple
    columns: tuple
    lift: tuple
    clock_constant: int
    from_q3: int
    pin: int

    def origin(self, K):
        base = 6 * K + 5
        h0 = ((self.label[2] - self.clock_constant) // 2 - base) % 36
        return tuple(integer(row.eval(base, h0, 0) - K * v)
                     for row, v in zip(self.checkpoint, V))

    @property
    def group(self):
        return self.word, self.columns, self.lift


def integer(x):
    assert x == int(x), x
    return int(x)


def checkpoint_from(q3, block, guard):
    for k in (0, 7):
        start = translate(q3, k, SECOND_DRIFT)
        assert guard.walk(start, block) == translate(start, 1, SECOND_DRIFT)
    states = [translate(q3, 8, SECOND_DRIFT)]
    orders = []
    for _ in range(160):
        order = guard.paired_order(states[-1])
        orders.append(order)
        states.append(guard.walk(states[-1], (order,)))
        for period in (10, 20):
            if len(orders) < 2 * period or orders[-period:] != orders[-2 * period:-period]:
                continue
            detection = tuple((period // 10) * d for d in DETECTION_DRIFT)
            before = states[-period - 1]
            if before != translate(states[-2 * period - 1], 1, detection):
                continue
            current = states[-1]
            if current != translate(before, 1, detection):
                continue
            tail = tuple(orders[-period:])
            phase = min(range(period), key=lambda k: tail[k:] + tail[:k])
            word = tail[phase:] + tail[:phase]
            checkpoint = guard.walk(current, tail[:phase])
            after = guard.walk(checkpoint, word)
            delta = tuple(integer(row.terms[-1] / Q(SCALE)) for row in difference(after, checkpoint))
            assert after == translate(checkpoint, 1, delta)
            far = translate(checkpoint, B / (500 * period), delta)
            assert guard.walk(far, word) == translate(far, 1, delta)
            assert guard.nonnegative(3 * B / 20 - (B - checkpoint[0]))
            for endpoint in (checkpoint, translate(checkpoint, 10, delta)):
                assert guard.nonnegative(endpoint[0] - 5 * B / 6 - 85)
            return checkpoint, word, delta, 160 + len(orders) + phase
    raise AssertionError(("return was not found", q3))


def build_sections(records):
    sections = {}
    for z in (1, 2):
        for r1 in range(-4, 18):
            guard = BoxGuards(z, r1)
            record = records[f"{z},{r1}"]
            qR = (3 * (B - 1) / 2 + r1 - P, R, P, B - 1 + z - P, R - (B - 3) / 2)
            q2 = step_many(B, qR, record["orders_to_first_paired_block"])
            a, b, _, d, e = q2
            J = 2 * a + 2 * b - 2 * d + 4 * e - 3 * B
            kappa = record["second_clock_constant"]
            assert J == 2 * B + 2 * P - 4 * R + kappa and kappa % 2 == 1
            pin = r1 - 11 * (max(record["failed_times_in_this_prefix"]) // 10)
            for r2 in range(601, 672, 2):
                label = z, r1, r2
                q3 = translate(q2, (J - r2) / 72, SECOND_DRIFT)
                args = q3, record["paired_block"], guard
                checkpoint, word, delta, elapsed = checkpoint_from(*args)
                b, p, r, _ = zip(*(row.terms for row in checkpoint))
                columns = (
                    tuple(integer(Q(36 * x, SCALE)) for x in p),
                    tuple(integer(Q(2 * x + y, SCALE)) for x, y in zip(p, r)),
                )
                lift = tuple(integer(Q(36 * x, SCALE) - 6 * v) for x, v in zip(b, V))
                sections[label] = Section(label, word, checkpoint, delta, columns,
                                          lift, kappa, elapsed, pin)
    return sections


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def component(scale, vector):
    return tuple(scale * a for a in vector)


@cache
def first_minor(columns):
    a, b = columns
    for i, j in combinations(range(5), 2):
        determinant = a[i] * b[j] - a[j] * b[i]
        if determinant:
            return i, j, determinant
    raise AssertionError(("section has rank below two", columns))


def decompose(point, columns):
    a, b = columns
    i, j, determinant = first_minor(columns)
    u = Q(point[i] * b[j] - point[j] * b[i], determinant)
    v = Q(a[i] * point[j] - a[j] * point[i], determinant)
    normal = difference(point, add(component(u, a), component(v, b)))
    return (u, v), normal


def lattice_fraction(coords):
    return tuple(x % 1 for x in coords)


def orbit_coset(coords, delta):
    return min(lattice_fraction(add(coords, component(j, delta))) for j in range(11))


def witness(before, after, K, shift):
    offsets = difference(after.origin(K + shift), before.origin(K))
    coords, normal = decompose(offsets, before.columns)
    assert not any(normal)
    delta, _ = decompose(before.displacement, before.columns)
    for j in range(11):
        u, v = difference(coords, component(j, delta))
        if u.denominator == v.denominator == 1:
            return int(u), int(v), j
    raise AssertionError(("matching key without integer witness", before.label, after.label, K, shift))


def keys_for(sections, K, delta_by_group):
    keys = defaultdict(list)
    for label, section in sections.items():
        coords, normal = decompose(section.origin(K), section.columns)
        keys[section.group, normal, orbit_coset(coords, delta_by_group[section.group])].append(label)
    return keys


def compare_keys(sections, before, after, K, shift):
    pairs = []
    missing_before, missing_after = [], []
    for key in sorted(before.keys() | after.keys()):
        left, right = before.get(key, []), after.get(key, [])
        pairs.extend((a, b, *witness(sections[a], sections[b], K, shift)) for a, b in zip(left, right))
        missing_before.extend(left[len(right):])
        missing_after.extend(right[len(left):])
    return pairs, sorted(missing_before), sorted(missing_after)


def cycle_profiles(keys, delta_by_group, sigma):
    cycles = set()
    for group, normal, coords in keys:
        delta = delta_by_group[group]
        rotations = [orbit_coset(add(coords, component(j, sigma)), delta) for j in range(4)]
        phase = min(range(4), key=rotations.__getitem__)
        rotated = tuple(rotations[phase:] + rotations[:phase])
        assert len(set(rotated)) == 4
        cycles.add((group, normal, rotated))
    return Counter(tuple(len(keys.get((group, normal, p), ())) for p in phases)
                   for group, normal, phases in cycles)


def display(word):
    return "".join(ORDER_NAMES[order] for order in word)


def enumerate_quotient(sections):
    delta_by_group = {}
    word_counts = Counter()
    word_details = set()
    shift_residues = set()
    sigma = (-Q(1, 44), -Q(48, 11))
    for section in sections.values():
        delta, normal = decompose(section.displacement, section.columns)
        assert not any(normal), (section.label, normal)
        assert lattice_fraction(component(11, delta)) == (0, 0)
        assert len({lattice_fraction(component(j, delta)) for j in range(11)}) == 11
        assert delta_by_group.setdefault(section.group, delta) == delta
        word_counts[section.word] += 1
        word_details.add((section.word, section.columns, section.displacement, delta))
        for K in range(6):
            assert difference(section.origin(K + 6), section.origin(K)) == section.lift
            coords, normal = decompose(difference(section.origin(K + 1), section.origin(K)),
                                       section.columns)
            assert not any(normal), (section.label, K, normal)
            shift_residues.add(lattice_fraction(coords))
        zero = orbit_coset((Q(0), Q(0)), delta)
        assert [j for j in range(1, 5) if orbit_coset(component(j, sigma), delta) == zero] == [4]
    assert shift_residues == {lattice_fraction(sigma)}, shift_residues

    key_table = [keys_for(sections, k, delta_by_group) for k in range(10)]
    match_tables, match_summaries = {}, {}
    for shift in (1, 2, 3, 4):
        match_tables[shift] = {}
        match_summaries[str(6 * shift)] = rows = []
        for K in range(6):
            pairs, a, b = compare_keys(sections, key_table[K], key_table[K + shift], K, shift)
            match_tables[shift][K] = pairs
            hist = sorted(Counter(pair[-1] for pair in pairs).items())
            rows.append({"K_mod_6": K, "matched": len(pairs), "j_counts": hist,
                         "left_residual": a, "right_residual": b,
                         "left_orphan_count": sum(sections[label].pin not in range(-4, 6) for label in a),
                         "right_orphan_count": sum(sections[label].pin not in range(-4, 6) for label in b)})
    cycle_summaries = []
    for K in range(6):
        profiles = cycle_profiles(key_table[K], delta_by_group, sigma)
        assert sum(profiles.values()) > 0
        assert sum(sum(profile) * count for profile, count in profiles.items()) == len(sections)
        assert all(a == c and b == d for a, b, c, d in profiles), profiles
        cycle_summaries.append({
            "K_mod_6": K,
            "occupied": sum(profiles.values()),
            "imbalance_histogram": sorted(sum_counts(profiles).items()),
            "profiles": sorted(profiles.items()),
        })

    summary = {
        "sections": len(sections),
        "phase_displacements": [
            {"word": display(word), "classes": word_counts[word], "columns": columns,
             "delta": delta, "plane_coordinates": coords}
            for word, columns, delta, coords in sorted(word_details)
        ],
        "base_shift_residue": sorted(shift_residues),
        "base_shift_order": 4,
        "matches": match_summaries,
        "cycles": cycle_summaries,
    }
    return match_tables, summary


def sum_counts(profiles):
    histogram = Counter()
    for (a, b, _, _), count in profiles.items():
        histogram[abs(a - b)] += count
    return histogram


def in_box(base, q):
    a, b, c, d, e = q
    chart = (a, b, c, 2 * (b - e) - base, base - c - d)
    return all(lo * base <= x <= hi * base for x, (lo, hi) in zip(chart, BOX))


def source_from(base, section, K, h, r, N, n, m):
    z, r1, r2 = section.label
    p0 = ((r2 - section.clock_constant) // 2 - base) % 36
    p = p0 + 36 * h + 2 * r
    c = p + 3 * n + 2 * N
    b = r + 4 * n + 2 * m + 4 * N
    a = 3 * (base - 1) // 2 + r1 - c + 8 * n + 10 * m + 22 * N
    d = base - c - 1 + 5 * n + z - 2 * m
    e = b - 3 * n - (base - 3) // 2
    q = a, b, c, d, e
    checkpoint = add(section.origin(K), add(component(h, section.columns[0]),
                                           component(r, section.columns[1])))
    return q, checkpoint


def random_interior(base, rng):
    for lo, hi in BOX:
        margin = (hi - lo) / 10
        lower = -(-(lo + margin) * base // 1)
        upper = (hi - margin) * base // 1
        yield rng.randint(lower, upper)


def draw_pair(K, shift, a, b, s, t, rng):
    base = 6 * K + 5
    other = base + 6 * shift
    for _ in range(1_000):
        A, B2, C, H, G = random_interior(base, rng)
        p, r, N = prefix_coordinates(Q(base), Q(A), Q(B2), Q(C), Q(H), Q(G), *a.label[:2])
        h = round((p - 2 * r - ((a.label[2] - a.clock_constant) // 2 - base) % 36) / 36)
        r, N = round(r), round(N)
        n = round(Q(H + 3, 6))
        m = round(Q(G - 1 + 5 * n + a.label[0], 2))
        q1, w1 = source_from(base, a, K, h, r, N, n, m)
        q2, w2 = source_from(other, b, K + shift, h - s, r - t, N, n, m)
        if in_box(base, q1) and in_box(other, q2):
            return q1, w1, q2, w2
    raise AssertionError(("no interior pair found", K, shift, a.label, b.label))


def replay_endpoint(base, source, section, checkpoint, predict, centered_step):
    output = predict(base, source)
    label = output["first_overshoot"], output["first_remainder"], output["second_remainder"]
    assert label == section.label, (label, section.label)
    x = tuple(6 * q - base * v for q, v in zip(output["bridge_state"], V))
    for _ in range(section.from_q3):
        x = centered_step(x)
    assert x == tuple(6 * w - 5 * v for w, v in zip(checkpoint, V))
    at_checkpoint = output["bridge_start"] + section.from_q3
    assert output["difference_time"] > at_checkpoint
    return {
        "source": source,
        "checkpoint_time": at_checkpoint,
        "remaining_difference_steps": output["difference_time"] - at_checkpoint,
        "period": output["period"],
        "cycle": output["first_cyclic_centered_state"],
        "word_time": output["word_time"],
    }


def exact_replays(sections, match_tables):
    from certify import centered_step, predict

    seed = 0x054B3C01
    rng = Random(seed)
    results = []
    for scale in (204_000, 400_002):
        for residue in range(6):
            K = scale + (residue - scale) % 6
            for shift in (1, 2):
                pairs = match_tables[shift][residue]
                zero = [pair for pair in pairs if pair[-1] == 0]
                nonzero = [pair for pair in pairs if pair[-1] != 0]
                for a_label, b_label, s, t, j in [rng.choice(zero), *rng.sample(nonzero, 2)]:
                    a, b = sections[a_label], sections[b_label]
                    q1, w1, q2, w2 = draw_pair(K, shift, a, b, s, t, rng)
                    assert w2 == add(w1, component(j, a.displacement))
                    base = 6 * K + 5
                    left = replay_endpoint(base, q1, a, w1, predict, centered_step)
                    right = replay_endpoint(base + 6 * shift, q2, b, w2, predict, centered_step)
                    assert left["cycle"] == right["cycle"] and left["period"] == right["period"]
                    assert left["remaining_difference_steps"] - right["remaining_difference_steps"] == j * len(a.word)
                    results.append({"base": base, "shift": 6 * shift, "labels": (a_label, b_label),
                                    "witness": (s, t, j), "left": left, "right": right})
        print(f"Exact predictor completed scale K~{scale}", file=sys.stderr, flush=True)
    return {"seed": seed, "pairs": results}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, default=str, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    evidence, output_dir = map(Path, sys.argv[1:])
    sections = build_sections(json.loads(evidence.read_text())["residues"])
    print(f"Built {len(sections)} sections from the certified short passage", file=sys.stderr, flush=True)
    matches, statistics = enumerate_quotient(sections)
    write_json(output_dir / "quotient.json", statistics)
    print("Orbit quotient written", file=sys.stderr, flush=True)
    write_json(output_dir / "replays.json", exact_replays(sections, matches))
