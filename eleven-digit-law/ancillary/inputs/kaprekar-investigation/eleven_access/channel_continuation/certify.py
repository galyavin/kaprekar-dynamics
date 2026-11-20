"""Rational endpoint certification and exact capture for the Task 029 channel."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import json


BOX = tuple((F(lo, 100000), F(hi, 100000)) for lo, hi in
            ((94040, 94130), (73770, 73850), (57940, 58015),
             (20, 80), (200, 300)))
BASE_MIN = 1000001
V = (5, 4, 3, 2, 1)

S = "5061728394a"
U = "5061782394a"
FAILURE = "5061273894a"
W = "5067182934a"
FIRST_WORD = (FAILURE, U, S, S, S, S, U, S, S, S) * 2
FIRST_DELTA = (-20, -4, -2, 2, -4)
SECOND_DELTA = (-20, -4, 0, 4, -4)
R1_RANGE = range(-4, 18)
R2_RANGE = range(601, 672, 2)


@dataclass(frozen=True)
class Affine:
    coefficients: tuple

    @staticmethod
    def constant(value):
        return Affine((0, 0, 0, 0, 0, 0, value))

    def __add__(self, rhs):
        if not isinstance(rhs, Affine):
            rhs = self.constant(rhs)
        return Affine(tuple(x + y for x, y in zip(self.coefficients, rhs.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, rhs):
        return self + -rhs

    def __rsub__(self, lhs):
        return -self + lhs

    def __mul__(self, factor):
        return Affine(tuple(factor * x for x in self.coefficients))

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        return self * F(1, divisor)

    def normalized_minimum(self):
        return self.coefficients[5] + sum(
            weight * (lo if weight >= 0 else hi)
            for weight, (lo, hi) in zip(self.coefficients, BOX)
        )

    def at(self, values):
        return sum(weight * value for weight, value in zip(self.coefficients, (*values, 1)))


A, B2, C, H0, G0, BASE = tuple(
    Affine(tuple(int(i == j) for i in range(7))) for j in range(6)
)
SOURCE = (A, B2, C, BASE - C - G0, B2 - (BASE + H0)/2)
MIDPOINT = tuple((lo + hi)*BASE_MIN/2 for lo, hi in BOX) + (BASE_MIN,)


def full_gaps(base, q):
    return (base - q[0], *(x-y for x, y in zip(q, q[1:])), q[-1])


def digits(base, q):
    a, b, c, d, e = q
    return (a, b, c, d, e - 1, base - 1, base - e - 1,
            base - d - 1, base - c - 1, base - b - 1, base - a)


def ordered_digits(base, q, word):
    values = digits(base, q)
    return tuple(values[int(index, 16)] for index in word)


def apply_order(base, q, word):
    ordered = ordered_digits(base, q, word)
    return tuple(ordered[i] - ordered[-i-1] for i in range(5))


def known_sorted(values):
    if isinstance(values[0], Affine):
        scores = [value.at(MIDPOINT) for value in values]
    else:
        scores = values
    return ''.join(format(i, 'x') for i in sorted(range(11), key=lambda i: (-scores[i], i)))


def nonnegative(row, label):
    if not isinstance(row, Affine):
        assert row >= 0, (label, row)
        return
    slope = row.normalized_minimum()
    constant = row.coefficients[-1]
    assert slope >= 0 and slope*BASE_MIN + constant >= 0, (label, slope, constant, row)


def shift(q, count, delta):
    return tuple(x + count*y for x, y in zip(q, delta))


def walk(base, q, words, label):
    for phase, word in enumerate(words):
        for gap in full_gaps(base, q):
            nonnegative(gap - 1, (label, phase, 'physical'))
        ordered = ordered_digits(base, q, word)
        for upper, lower in zip(ordered, ordered[1:]):
            nonnegative(upper-lower, (label, phase, word))
        q = apply_order(base, q, word)
    return q


def incoming(base, q, z=None):
    a, b, c, d, e = q
    h = 2*(b-e) - base + 1
    g = base - c - d - 1
    n = (h+2)/6 if z is not None else (h+2)//6
    prep = shift(q, n, (-5, -4, -3, -2, -1))
    g1 = g + 5*n
    if z is None:
        z = 2 - g1 % 2
        m = (g1+z)//2
    else:
        m = (g1+z)/2
    return shift(prep, m, (-10, -2, 0, 2, -2)), 5*n + 10*m, z


def first_clock(base, q):
    return q[0] + q[2] - (3*base-3)*F(1, 2)


def second_clock(base, q):
    a, b, _, d, e = q
    return 2*a + 2*b - 2*d + 4*e - 3*base


def first_switch(base, q):
    """Run the bounded remainder until the first W order, without taking it."""
    orders = []
    states = []
    for _ in range(26):
        word = known_sorted(digits(base, q))
        if word == W:
            return q, tuple(orders), states
        states.append(q)
        orders.append(word)
        q = walk(base, q, (word,), "first switch")
    raise AssertionError("First switch takes more than 25 steps")


def second_word(base, q):
    orders = []
    for _ in range(20):
        word = known_sorted(digits(base, q))
        orders.append(word)
        q = walk(base, q, (word,), "second word discovery")
    return tuple(orders), q


def failure_slack(base, q, word):
    if word != FAILURE:
        assert word in (S, U, W)
        return 0
    row = base - q[2] - q[3] - 1
    if isinstance(row, Affine):
        assert row == Affine.constant(row.coefficients[-1])
        row = row.coefficients[-1]
    assert row in (0, -1, -2)
    return row


def check_bridge(base, q, label):
    gaps = full_gaps(base, q)
    ceiling_upper = base/7 + 65
    bridge_upper = ceiling_upper - gaps[0]
    nonnegative(base/7 + 64 - gaps[0], (label, 'bridge nonnegative'))
    nonnegative(base/250 - bridge_upper, (label, 'bridge upper bound'))
    for i, gap in enumerate(gaps):
        nonnegative(gap - base*F(138, 1000), (label, 'bridge gap low', i))
        nonnegative(base*F(190, 1000) - gap, (label, 'bridge gap high', i))
    for i, gap in enumerate(gaps[1:], 1):
        nonnegative(gap - ceiling_upper - bridge_upper/3 - 3,
                    (label, 'bridge non-anchor margin', i))


def correction_control():
    for overshoot in (1, 2):
        q = (A, B2, C, BASE-C-1+overshoot, B2-(BASE-3)/2)
        reference = apply_order(BASE, q, S)
        for j in (0, 1, 2, 3):
            translated = shift(q, j, V)
            actual_order = FAILURE if j < overshoot else S
            actual = apply_order(BASE+6*j, translated, actual_order)
            expected = shift(reference, j, V)
            assert actual == shift(expected, max(0, overshoot-j), (0, 0, 0, 1, 0))


def certify():
    correction_control()
    section_words = Counter()
    switch_lengths = Counter()
    last_offsets = []
    fallbacks = {}
    for z, remainder in product((1, 2), R1_RANGE):
        q, _, _ = incoming(BASE, SOURCE, z)
        u = first_clock(BASE, q)
        count = (u - remainder)/22
        nonnegative(count - 1, ('first clock count', z, remainder))
        for end, j in (('first', 0), ('last', count - 1)):
            at = shift(q, j, FIRST_DELTA)
            assert walk(BASE, at, FIRST_WORD, (end, z, remainder)) == shift(at, 1, FIRST_DELTA)
            halfway = walk(BASE, at, FIRST_WORD[:10], 'failure slack at half block')
            assert failure_slack(BASE, at, FAILURE) == -z
            assert failure_slack(BASE, halfway, FAILURE) == z - 2
        at = shift(q, count, FIRST_DELTA)
        assert first_clock(BASE, at) == Affine.constant(remainder)
        q2, transition, states = first_switch(BASE, at)
        assert all(word[0] == '5' and word[-1] == 'a' for word in transition)
        failure_times = [i for i, (state, word) in enumerate(zip(states, transition))
                         if failure_slack(BASE, state, word) < 0]
        assert failure_times[0] == 0
        last_offsets.append(F(10 + 5*z - 10*remainder, 11) + failure_times[-1])
        switch_lengths[len(transition)] += 1

        word, end_q2 = second_word(BASE, q2)
        assert set(word) <= {S, U, W}
        assert end_q2 == shift(q2, 1, SECOND_DELTA)
        section_words[''.join({S: 'S', U: 'U', W: 'W'}[order] for order in word)] += 1

        clock = second_clock(BASE, q2)
        expected = 2*BASE + 2*at[2] - 4*at[1]
        residue = clock - expected
        assert residue == Affine.constant(residue.coefficients[-1])
        assert residue.coefficients[-1] % 2 == 1
        for end_remainder in (R2_RANGE.start, R2_RANGE[-1]):
            loops = (clock - end_remainder)/72
            nonnegative(loops - 1, ('second clock count', z, remainder, end_remainder))
            last = shift(q2, loops - 1, SECOND_DELTA)
            assert walk(BASE, last, word, ('second last', z, remainder, end_remainder)
                        ) == shift(last, 1, SECOND_DELTA)
            q3 = shift(q2, loops, SECOND_DELTA)
            assert second_clock(BASE, q3) == Affine.constant(end_remainder)
            check_bridge(BASE, q3, (z, remainder, end_remainder))

        fallbacks[f'{z},{remainder}'] = {
            'orders_to_first_paired_block': transition,
            'failed_times_in_this_prefix': failure_times,
            'paired_block': word,
            'second_clock_constant': int(residue.coefficients[-1]),
        }
    return {'base_minimum': BASE_MIN, 'first_remainders_checked': len(fallbacks),
            'second_remainder_endpoints': [R2_RANGE.start, R2_RANGE[-1]],
            'second_remainders_covered': len(R2_RANGE),
            'paired_blocks': dict(section_words),
            'short_prefix_lengths': dict(sorted(switch_lengths.items())),
            'minimum_last_time_offset': str(min(last_offsets)),
            'maximum_last_time_offset': str(max(last_offsets)),
            'residues': fallbacks}


def centered_from(base, q):
    return tuple(6*x-base*v for x, v in zip(q, V))


def physical_from(base, x):
    assert all((base*v+a) % 6 == 0 for a, v in zip(x, V))
    return tuple((base*v+a)//6 for a, v in zip(x, V))


def centered_step(x):
    ordered = [-6]
    for i in range(5):
        pair = (x[i] - 6*(i == 4), -x[4-i] - 6*(i < 4))
        ordered.extend(sorted(pair, reverse=True))
    return tuple(ordered[i] - ordered[10-i] for i in range(5))


def return_to_cycle(x):
    seen = {}
    ticks = 0
    while x not in seen:
        seen[x] = ticks
        x = centered_step(x)
        ticks += 1
    return seen[x], ticks - seen[x], x


def predict(base, source):
    """Evaluate the certified formula on an original difference vector."""
    assert base >= BASE_MIN and base % 6 == 5
    a, b, c, d, e = source
    chart = (a, b, c, 2*(b-e)-base, base-c-d)
    assert all(lo*base <= x <= hi*base for x, (lo, hi) in zip(chart, BOX))
    assert (b-e) % 3 == 1

    q, first_failure, z = incoming(base, source)
    u = first_clock(base, q)
    loops, remainder = divmod(int(u) + 4, 22)
    q = shift(q, loops, FIRST_DELTA)
    q2, transition, states = first_switch(base, q)
    elapsed = first_failure + 20*loops
    last_failure = elapsed + max(
        i for i, (state, word) in enumerate(zip(states, transition))
        if failure_slack(base, state, word) < 0
    )

    second_loops, second_remainder = divmod(second_clock(base, q2) - 601, 72)
    q3 = shift(q2, second_loops, SECOND_DELTA)
    bridge_start = elapsed + len(transition) + 20*second_loops
    assert q3[0] == a-bridge_start
    target_gap = (base+7)//7 + 64
    bridge_steps = target_gap - (base-q3[0])
    state = centered_from(base, q3)
    for _ in range(bridge_steps):
        state = centered_step(state)
    source_time = bridge_start + bridge_steps
    assert source_time == a-base+target_gap
    assert min(full_gaps(base, physical_from(base, state))) == target_gap

    terminal_steps, period, cyclic = return_to_cycle(state)
    cyclic_q = physical_from(base, cyclic)
    return {
        'first_failure': first_failure,
        'last_failure': last_failure,
        'first_overshoot': z,
        'first_remainder': remainder-4,
        'second_remainder': second_remainder+601,
        'second_block_start': elapsed + len(transition),
        'bridge_start': bridge_start,
        'bridge_state': q3,
        'second_block_state': q2,
        'source_time': source_time,
        'source_state': physical_from(base, state),
        'centered_steps': terminal_steps,
        'difference_time': source_time + terminal_steps,
        'word_time': source_time + terminal_steps + 1,
        'period': period,
        'first_cyclic_centered_state': cyclic,
        'first_cyclic_difference': cyclic_q,
        'first_cyclic_word': digits(base, cyclic_q),
        'first_cyclic_words_difference': physical_from(base, centered_step(cyclic)),
    }


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run verification with assertions enabled")
    print(json.dumps(certify(), indent=2))
