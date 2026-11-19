"""Exact guarded affine verification on the residual box chart."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cache
from itertools import product
import json
import sys
from math import ceil, floor, gcd, lcm
from pathlib import Path

import certificate as prefix
import channel as kap
import cohort1 as orbit
import quotient as box


# Affine operations stay on exact integer numerators. Every division checks
# that its result fits this scale; the verifier rejects unsupported divisors.
SCALE = 90_000 * lcm(*range(1,97))**4


@dataclass(frozen=True)
class Form:
    """Affine form in (base, third q_R coordinate, second q_R coordinate)."""
    v: tuple

    @classmethod
    def create(cls, a=(0, 0, 0), c=0):
        return cls(tuple(orbit.integer(x * SCALE) for x in (*a, c)))

    @property
    def linear(self):
        return self.v[:3]

    def __add__(self, other):
        a,b,c,d = self.v
        if isinstance(other, Form):
            e,f,g,h = other.v
            return Form((a+e,b+f,c+g,d+h))
        return Form((a,b,c,d+orbit.integer(SCALE*other)))

    __radd__ = __add__

    def __mul__(self, value):
        if isinstance(value, int):
            return Form(tuple(n*value for n in self.v))
        return self.__mul__(value.numerator).__truediv__(value.denominator)

    __rmul__ = __mul__

    def __neg__(self):
        a,b,c,d = self.v
        return Form((-a,-b,-c,-d))

    def __sub__(self, other):
        if isinstance(other, Form):
            a,b,c,d = self.v
            e,f,g,h = other.v
            return Form((a-e,b-f,c-g,d-h))
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __truediv__(self, value):
        assert all(n%value == 0 for n in self.v), value
        return Form(tuple(n//value for n in self.v))

    def scaled_at(self, point):
        a,b,c,d = self.v
        B,p,r = point
        return a*B+b*p+c*r+d

    def at(self, point):
        return Q(self.scaled_at(point), SCALE)

    def constant(self):
        assert not any(self.linear), self
        return Q(self.v[3], SCALE)

    def leading(self):
        return Form((*self.linear, 0))

    def record(self):
        return [str(Q(x, SCALE)) for x in self.v]


BASE = Form.create((1, 0, 0))
NEG_SIX = Form.create(c=-6)
V = kap.V
CORNERS = tuple(product(*kap.BOX))
MID = tuple((lo + hi) / 2 for lo, hi in kap.BOX)
START = 1_000_001
PROBE_BASE = 36_000_000_000
PHASE_LIMIT = 400
LANES = 'S12L'
FROZEN = ((542, 52, 82, 380), (598, 52, 94, 312))
# Prior classification separators at the chart midpoint.
MID_SEPARATORS = (Q(1160,10000), Q(1177,10000), Q(1205,10000))
CENTER = (940800, 738087, 579763, 484, 2491)
SCALES = (446, 365, 342, 280, 449)
COEFFICIENTS = (590, -307, -447, 127, -296)


def centered(row, v):
    form = Form.create(tuple(Q(x, orbit.SCALE) for x in row.terms[:3]), Q(row.terms[3], orbit.SCALE))
    return 6 * form - v * BASE


def shift(state, count, drift):
    return tuple(x + count * d for x, d in zip(state, drift))


def diff(left, right):
    return tuple(a - b for a, b in zip(left, right))


def offsets(state):
    return tuple((x - 6 * (i == 4), -state[4-i] - 6 * (i < 4))
                 for i, x in enumerate(state))


def step(state, word):
    rows = (NEG_SIX, *(v for pair, swap in zip(offsets(state), word)
                          for v in (pair[::-1] if swap else pair)))
    return tuple(rows[i] - rows[10-i] for i in range(5))


def comparisons(state, word):
    """For a recorded order, weak comparisons suffice when digit rows tie."""
    return tuple((b - a if swap else a - b) for (a, b), swap in zip(offsets(state), word))


def paired_order(order):
    assert order[0] == '5'
    for i in range(5):
        assert set(order[2*i+1:2*i+3]) == {format(i, 'x'), format(i+6, 'x')}
    return tuple(int(order[2*i+1], 16) != i for i in range(5))


def word_record(word):
    return [''.join(str(int(swap)) for swap in order) for order in word]


def rotate(word):
    return min(word[i:] + word[:i] for i in range(len(word)))


@dataclass(frozen=True)
class Lattice:
    """Parameterization of the residue class in the three route coordinates."""
    point: tuple
    columns: tuple

    def residue(self, row):
        offset = row.scaled_at(self.point)
        increments = [sum(x*y for x,y in zip(row.linear,col)) for col in self.columns]
        if all(x%SCALE == 0 for x in increments):
            return Q(offset%SCALE, SCALE)
        return None

    def split(self, row):
        """Partition once by the image of the quotient in the integers mod 1."""
        increments = [Q(sum(x*y for x, y in zip(row.linear, col)), SCALE) for col in self.columns]
        denominator = lcm(*(x.denominator for x in increments))
        numerators = [int(x*denominator) for x in increments]
        columns = list(self.columns)
        for j in (1, 2):
            a, b = numerators[0], numerators[j]
            if not b:
                continue
            g, u, v = extended_gcd(a, b)
            c, d = columns[0], columns[j]
            columns[0] = tuple(u*x + v*y for x, y in zip(c, d))
            columns[j] = tuple((-b//g)*x + (a//g)*y for x, y in zip(c, d))
            numerators[0] = g
        modulus = denominator // gcd(denominator, numerators[0])
        first = columns[0]
        columns[0] = tuple(modulus*x for x in first)
        for k in range(modulus):
            yield Lattice(tuple(x + k*y for x, y in zip(self.point, first)), tuple(columns))

    def record(self):
        return {'point': self.point, 'columns': self.columns}

    def base_phases(self, base):
        modulus = gcd(72, *(c[0] for c in self.columns))
        return tuple(phase for phase in (0,1) if (base+36*phase-self.point[0]) % modulus == 0)

    def probe(self, bounds):
        """Round the real-chart midpoint in the current lattice coordinates."""
        p_mid, r_mid = bounds.coordinates(MID, PROBE_BASE)
        target = (PROBE_BASE, p_mid, r_mid)
        offset = tuple(x-y for x, y in zip(target, self.point))
        columns = self.columns
        divisor = det(*columns)
        coordinate = [round(det(*(columns[:i] + (offset,) + columns[i+1:])) / Q(divisor))
                      for i in range(3)]
        return tuple(x + sum(k*c[i] for k, c in zip(coordinate, columns))
                     for i, x in enumerate(self.point))


def det(c, d, e):
    return sum(c[i] * (d[(i+1)%3]*e[(i+2)%3] - d[(i+2)%3]*e[(i+1)%3]) for i in range(3))


def extended_gcd(a, b):
    if b == 0:
        return abs(a), (1 if a >= 0 else -1), 0
    g, x, y = extended_gcd(b, a % b)
    return g, y, x-(a//b)*y


@cache
def coefficient_bounds(a):
    return orbit.normalized_bounds(*(Q(x,SCALE) for x in a))


@cache
def scaled_minimum(a):
    low = coefficient_bounds(a)[0]*SCALE
    return int(low) if low.denominator == 1 else low


@cache
def input_bounds(cohort, z, r1):
    if cohort == 1:
        return orbit.BoxGuards(z, r1)
    return box.Cohort2Guards(prefix.first_failure(cohort, z)[0], z, r1)


class Guards:
    def __init__(self, cohort, z, r1):
        origin = input_bounds(cohort, z, r1)
        self.p0, self.r0 = Q(origin.p0), Q(origin.r0)
        self.intercept_den = lcm(self.p0.denominator, self.r0.denominator)
        self.p_int = orbit.integer(self.p0*self.intercept_den)
        self.r_int = orbit.integer(self.r0*self.intercept_den)
        self.cohort, self.z, self.r1 = cohort, z, r1
        self.minimum_base = START
        self.current = None
        self.lattice = None

    def coordinates(self, normalized, base):
        p, r, _ = orbit.prefix_coordinates(base, *(x*base for x in normalized), self.z, self.r1)
        p0, r0 = orbit.chart_constants(self.z, self.r1)
        return p + self.p0 - p0, r + self.r0 - r0

    def parts(self, guard):
        lo, hi = coefficient_bounds(guard.linear)
        _, b, r, c = guard.v
        constant = Q(c + b*self.p0 + r*self.r0, SCALE)
        return lo, hi, constant

    def nonnegative(self, guard, location):
        lo = scaled_minimum(guard.linear)
        _, b, r, c = guard.v
        intercept = c*self.intercept_den + b*self.p_int + r*self.r_int
        if lo > 0:
            slope = lo*self.intercept_den
            if slope*self.minimum_base + intercept < 0:
                self.minimum_base = ceil(Q(-intercept, slope))
            return
        if lo == 0 and intercept >= 0:
            return
        lo, hi, constant = self.parts(guard)
        corner_values = []
        for x in CORNERS:
            p, r = self.coordinates(x, 1)
            corner_values.append((guard.at((1, p, r)) - constant, x))
        assert min(corner_values)[0] == lo and max(corner_values)[0] == hi
        raise GuardFailure({
            'class': self.current,
            'location': location,
            'lattice': self.lattice.record(),
            'affine_B_p_r_1': guard.record(),
            'chart_intercept': str(constant),
            'normalized_slope_range': [str(lo), str(hi)],
            'min_corner': [str(x) for x in min(corner_values)[1]],
            'max_corner': [str(x) for x in max(corner_values)[1]],
        })


class GuardFailure(Exception):
    pass


def propose(state, probe):
    return tuple(a.scaled_at(probe) < b.scaled_at(probe) for a, b in offsets(state))


@cache
def displacement_slopes(word, drift):
    """Propagate a displacement through the linear part of a return word."""
    current = drift
    slopes = []
    for order in word:
        pairs = tuple((x, -current[4-i]) for i, x in enumerate(current))
        slopes.append(tuple((b-a if swap else a-b) for (a,b),swap in zip(pairs,order)))
        rows = (0, *(v for pair, swap in zip(pairs, order)
                     for v in (pair[::-1] if swap else pair)))
        current = tuple(rows[i]-rows[10-i] for i in range(5))
    assert current == drift
    return tuple(slopes)


def repeat_guards(state, word, drift):
    original = state
    checks = []
    for phase, (order, slopes) in enumerate(zip(word, displacement_slopes(word, drift))):
        checks.extend((guard,slope,phase,pair) for pair,(guard,slope) in
                      enumerate(zip(comparisons(state,order),slopes)))
        state = step(state, order)
    assert state == shift(original, 1, drift)
    return checks


def repeat_candidates(checks):
    return [(guard/-slope,phase,pair) for guard,slope,phase,pair in checks if slope < 0]


@cache
def capture(point):
    return kap.return_to_cycle(point)


def terminal_time_bounds(state, word, drift, checks, limits, control, clock, guards):
    """Stop a fixed number of repetitions before a uniformly bounded tail."""
    at_root = shift(state, control, drift)
    if any(any(row.linear) for row in at_root):
        return None

    # Truncating at floor(control)-offset keeps a fixed distance from the root,
    # and keeps every guard of the last repetition nonnegative. Positive
    # leading slacks absorb finite intercepts after increasing the base cutoff.
    offset = 0
    for other, phase, pair in limits:
        lo, _, intercept = guards.parts(other-control)
        if lo == 0:
            offset = max(offset, ceil(-intercept))
        elif lo < 0:
            guards.nonnegative(other-control, ('terminal controller', phase, pair))
    before = control-offset
    guards.nonnegative(before-1, ('terminal loop positive',))
    for guard,slope,phase,pair in checks:
        guards.nonnegative(guard, ('terminal loop', 'first', phase, pair))
        guards.nonnegative(guard+before*slope, ('terminal loop', 'truncated', phase, pair))

    # At floor(control)-offset, the state is on a bounded line segment.
    # All physical points on it can be enumerated without further lattice splits.
    center = [row.constant() for row in shift(at_root, -offset, drift)]
    pivot = next(i for i, d in enumerate(drift) if d)
    a, b = sorted((center[pivot]-drift[pivot], center[pivot]))
    tail = []
    for x in range(ceil(a), floor(b)+1):
        fractional = Q(center[pivot]-x, drift[pivot])
        if not 0 <= fractional < 1:
            continue
        point = tuple(v - fractional*d for v, d in zip(center, drift))
        if all(q == int(q) and (int(q)+5*v) % 6 == 0 for q, v in zip(point, V)):
            tail.append(capture(tuple(int(q) for q in point)))
    assert tail

    # A state 5041 returns before the bounded endpoint and the endpoint
    # cannot both lie on a centered cycle: all its coordinates have absolute
    # value at most 2520. This supplies a uniform lower clock error bound.
    guards.nonnegative(before-5042, ('terminal lower clock',))
    return clock+len(word)*before, offset, tail


def walk_macro(state, clock, word, drift, lattice, guards, jumps, short=(), checks=None):
    """Traverse bounded switches and exact repeated words, including the capture tail."""
    guards.lattice = lattice
    depth = len(jumps)
    if depth >= 8:
        raise RuntimeError('too many macroscopic legs')
    probe = lattice.probe(guards)
    if checks is None:
        checks = repeat_guards(state, word, drift)
    limits = repeat_candidates(checks)
    assert limits
    control, phase, pair = min(limits, key=lambda item: item[0].scaled_at(probe))
    terminal = terminal_time_bounds(state, word, drift, checks, limits, control, clock, guards)
    if terminal is not None:
        end_time, offset, tail = terminal
        jumps = (*jumps, (word, control, tuple(drift), tuple(short), (phase, pair), offset))
        yield end_time, jumps, tail, lattice
        return
    remainder = lattice.residue(control)
    if remainder is None:
        # A finer congruence does not change the real chart; certify the same
        # path on every child, with that child's exact remainder.
        for child in lattice.split(control):
            yield from walk_macro(state, clock, word, drift, child, guards, jumps, short, checks)
        return
    count = control - remainder + 1
    guards.nonnegative(count-1, ('positive loop count', depth))
    end = shift(state, count, drift)
    for guard,slope,i,j in checks:
        guards.nonnegative(guard, ('loop', depth, 'first', i, j))
        guards.nonnegative(guard+(count-1)*slope, ('loop', depth, 'last', i, j))
    assert Q(0) <= remainder < 1
    # The controlling guard's next value is -speed * (1 - remainder).
    jumps = (*jumps, (word, count, tuple(drift), tuple(short), (phase, pair), remainder))
    clock += len(word) * count
    state = end
    short = []
    history = [state]
    for _ in range(PHASE_LIMIT):
        if not any(any(x.linear) for x in state):
            raise RuntimeError(('bounded state in a short transition', depth, state))
        order = propose(state, probe)
        for j, guard in enumerate(comparisons(state, order)):
            guards.nonnegative(guard, ('switch', depth, len(short), j))
        state = step(state, order)
        clock += 1
        short.append(order)
        history.append(state)
        for period in range(1, min(100, len(short)//2) + 1):
            if short[-period:] != short[-2*period:-period]:
                continue
            delta = diff(state, history[-period-1])
            if not any(x.v[3] for x in delta) or any(any(x.linear) for x in delta):
                continue
            if delta != diff(history[-period-1], history[-2*period-1]):
                continue
            proposed = tuple(short[-period:])
            d = tuple(orbit.integer(row.constant()) for row in delta)
            checks = repeat_guards(state, proposed, d)
            limits = repeat_candidates(checks)
            potential = min(candidate.scaled_at(probe) for candidate, _, _ in limits)
            if potential > 100*SCALE:
                yield from walk_macro(state, clock, proposed, d, lattice, guards, jumps, short, checks)
                return
    raise RuntimeError(('no macroscopic continuation or bounded terminal state', guards.current, depth, [x.record() for x in state], word_record(tuple(short[-60:]))))


def branches(entry, clock, section, lattice, guards):
    word = tuple(paired_order(w) for w in section.word)
    drift = tuple(6 * x for x in section.displacement)
    return walk_macro(entry, clock, word, drift, lattice, guards, ())


def initial_sections():
    quotient = json.loads(Path('quotient.json').read_text())['by_K_mod_6']
    labels = sorted(tuple(label) for label in quotient['0']['whole_box_left_residual'])
    assert len(labels) == len(set(labels)) == 176
    for k in ('0', '2', '4'):
        for side in ('left', 'right'):
            assert sorted(tuple(label) for label in quotient[k][f'whole_box_{side}_residual']) == labels
    zero_two = {(x['cohort'], x['z'], x['r1']): x for x in json.loads(Path('routes.json').read_text())}
    one = json.loads(Path('endpoint_verification.json').read_text())['residues']
    for cohort, z, r1, r2 in labels:
        guards = Guards(cohort, z, r1)
        if cohort == 2:
            record = zero_two[2,z,r1]
            switch, block, failures = prefix.word(record['switch']), prefix.word(record['paired']), record['failures']
        else:
            record = one[f'{z},{r1}']
            switch, block, failures = (record['orders_to_first_paired_block'],
                                       record['paired_block'], record['failed_times_in_this_prefix'])
        use_cohort = 2 if cohort == 2 else 1
        source = box.rows_at_switch(use_cohort, z, r1)
        bounds = input_bounds(cohort,z,r1)
        q2 = bounds.walk(source, switch)
        j = kap.second_clock(orbit.B, q2)
        count = (j-r2)/72
        assert bounds.nonnegative(count-1)
        q3 = orbit.translate(q2, count, orbit.SECOND_DRIFT)
        last = orbit.translate(q3, -1, orbit.SECOND_DRIFT)
        assert bounds.walk(q2, block) == orbit.translate(q2, 1, orbit.SECOND_DRIFT)
        assert bounds.walk(last, block) == q3
        section = box.build_section(use_cohort, (z,r1,r2), switch, block, failures, bounds)
        entry = tuple(centered(row, v) for row, v in zip(section.checkpoint, V))
        # The prefix certificate verifies that the original first coordinate decreases by
        # one at each step through q3. Add all explicit updates from q3 to q4.
        q3_a = (centered(q3[0], 5) + 5 * BASE)/6
        clock = -q3_a + section.from_q3
        yield (cohort, z, r1, r2), guards, section, entry, clock


def leading_coordinates(chart):
    p, r, _ = orbit.prefix_coordinates(Q(1), *chart, 1, 0)
    p0, r0 = orbit.chart_constants(1, 0)
    return 1, p-p0, r-r0


MID_COORDS = leading_coordinates(MID)
CORNER_COORDS = tuple(leading_coordinates(x) for x in CORNERS)


def theta(x):
    return Q(121049,10**6) + sum(Q(v, s)*(y-Q(c,10**6))
                                 for y, v, s, c in zip(x, COEFFICIENTS, SCALES, CENTER))


@cache
def profile(clock):
    at_mid = MID[0] + clock.at(MID_COORDS)
    lane = LANES[sum(at_mid >= t for t in MID_SEPARATORS)]
    corners = tuple(x[0] + clock.at(coords) for x, coords in zip(CORNERS, CORNER_COORDS))
    plateau = min((t-theta(x) if lane == 'L' else theta(x)-t) - Q(1,1000)
                  for x, t in zip(CORNERS, corners))
    return lane, corners, plateau, at_mid


def main():
    records = initial_sections()
    count = 0
    arithmetic_leaves = 0
    clock_rows = defaultdict(set)
    counts = defaultdict(Counter)
    lane_words = defaultdict(set)
    stage_words = defaultdict(set)
    classes = []
    greatest_delay = greatest_early = base_cutoff = 0
    plateau_min = Q(1)
    for label, guards, section, entry, clock in records:
        _, _, _, r2 = label
        for k in (0, 2, 4):
            for arm in (0, 1):
                base = 6 * (k + arm) + 5
                p0 = ((r2 - section.clock_constant)//2 - base) % 36
                for residue in range(6):
                    guards.current = (label, k, arm, residue)
                    start = Lattice((base, p0 + 2 * residue, residue), ((36,0,0),(0,36,0),(0,12,6)))
                    phase_lanes = defaultdict(set)
                    phase_profiles = defaultdict(set)
                    try:
                        for end, jumps, tail, lattice in branches(entry, clock, section, start, guards):
                            arithmetic_leaves += 1
                            assert len(jumps) in (2,3)
                            leading = end.leading()
                            lane, _, slack, _ = profile(leading)
                            assert slack > 0, (guards.current, end.record(), str(slack))
                            plateau_min = min(plateau_min, slack)
                            clock_rows[lane].add(leading)
                            lane_words[lane].add(rotate(jumps[1][0]))
                            for i, jump in enumerate(jumps):
                                stage_words[i].add(rotate(jump[0]))
                            for phase in lattice.base_phases(base):
                                phase_lanes[phase].add(lane)
                                phase_profiles[phase].add(tuple(leading.record()[:3]))
                            _, _, intercept = guards.parts(end)
                            greatest_delay = max(greatest_delay, ceil(intercept + 1 + max(t[0] for t in tail)))
                            greatest_early = max(greatest_early, ceil(len(jumps[-1][0])*5042 - intercept))
                    except GuardFailure as failure:
                        print(json.dumps({'completed': count, 'leaves': arithmetic_leaves,
                                          'counts': render_counts(counts), 'failure': failure.args[0]},
                                         indent=2))
                        return
                    assert all(len(phase_lanes[phase])==1 for phase in (0,1)), (guards.current,phase_lanes)
                    row = tuple(next(iter(phase_lanes[phase])) for phase in (0,1))
                    for phase, lane in enumerate(row):
                        counts[k,arm,phase][lane] += 1
                    classes.append((*guards.current, row, tuple(sorted(phase_profiles[0])),
                                    tuple(sorted(phase_profiles[1]))))
                    base_cutoff = max(base_cutoff, guards.minimum_base)
                    count += 1
                    if count % 360 == 0:
                        print(f'{count} primary class/buckets, {arithmetic_leaves} affine leaves',
                              file=sys.stderr, flush=True)
    gaps = []
    for low, high in zip(LANES, LANES[1:]):
        gap = min(y-x for a in clock_rows[low] for b in clock_rows[high]
                  for x,y in zip(profile(a)[1],profile(b)[1]))
        assert gap > 0, (low,high,gap)
        gaps.append(str(gap))
    summary = render_counts(counts)
    assert all(tuple(cnt) == FROZEN[arm] for (k, arm, phase), cnt in summary)
    assert all(len(v)>0 for v in lane_words.values())
    print(json.dumps({'completed': count, 'leaves': arithmetic_leaves, 'counts': summary,
                      'minimum_valid_base': base_cutoff, 'all_plateau_slack': str(plateau_min),
                      'max_clock_delay': greatest_delay, 'max_clock_early': greatest_early,
                      'adjacent_lane_gaps': gaps,
                      'lane_words': {lane: [word_record(w) for w in sorted(words)]
                                     for lane,words in lane_words.items()},
                      'stage_words': {str(i): len(words) for i,words in sorted(stage_words.items())},
                      'profiles': {lane: [r.record()[:3] for r in sorted(rows,key=lambda r: r.v)]
                                   for lane,rows in clock_rows.items()},
                      'classes': classes}, separators=(',', ':')))


def render_counts(counts):
    return [(key, [cnt[lane] for lane in LANES]) for key,cnt in sorted(counts.items())]


if __name__ == '__main__':
    main()
