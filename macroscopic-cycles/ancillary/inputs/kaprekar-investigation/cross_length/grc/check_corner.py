"""Exact 17-digit d2/d4 return and independent whole-word subtraction."""

from collections import Counter
from fractions import Fraction


LENGTH = 17
BLOCK = 244
START_OFFSET = (0, 0, 0, 0, 1, 0, 0, 0, -2)
DRIFT = (12, 12, 12, -15, 12, -15, -15, -15, 12)


def margin(gaps, row):
    return sum(gaps[-row:]) - sum(gaps[:row])


def translate(state, direction, count):
    return tuple(x + count * dx for x, dx in zip(state, direction))


def rotate(gaps):
    out = [gaps[0]] + [0] * 8
    for j in range(1, 9):
        out[min(2 * j, LENGTH - 2 * j)] = gaps[j]
    return tuple(out)


def paired_step(gaps, secondary=False):
    z = margin(gaps, 1)
    e = (z > 0) - (z < 0)
    h1 = max(0, -1 - z)
    h2 = max(0, -1 - margin(gaps, 2)) if secondary else 0
    h4 = max(0, -1 - margin(gaps, 4)) if secondary else 0
    kicks = (e - h1, -2 * e + h1, e + h1 - h2, -h1 + h2,
             h2, -h2, -h4, h4, 0)
    return translate(rotate(gaps), kicks, 1)


def literal_step(gaps, base):
    # The sorted representative consists of its eight differences and nine zeros.
    word = tuple(sum(gaps[i:]) for i in range(1, 9)) + (0,) * 9
    high = low = 0
    for digit, reflected in zip(word, reversed(word)):
        high = base * high + digit
        low = base * low + reflected
    value = high - low
    digits = []
    for _ in word:
        value, digit = divmod(value, base)
        digits.append(digit)
    assert value == 0
    digits.sort(reverse=True)
    differences = (base,) + tuple(digits[i] - digits[-1 - i] for i in range(8)) + (0,)
    return tuple(a - b for a, b in zip(differences, differences[1:]))


def offset_trace():
    trace = []
    offset = START_OFFSET
    for _ in range(BLOCK):
        trace.append(offset)
        offset = paired_step(offset)
    assert offset == translate(START_OFFSET, DRIFT, 1)
    assert max(abs(x) for state in trace for x in state) == 15
    assert 10 * ((1 - 4 * Fraction(1, 250) - 4 * Fraction(1, 90)) / 9
                 - (Fraction(3, 9) + 15) / 100000) > 1 + Fraction(1, 100000)
    return trace


def first_failures(trace):
    guarded_phases = [(t, margin(gaps, 2), margin(gaps, 4))
                      for t, gaps in enumerate(trace) if t % 4 == 0]
    assert [(t, b) for t, b, _ in guarded_phases if b <= -29] == [(228, -29)]
    assert [(t, b) for t, _, b in guarded_phases if b <= -53] == [(240, -53)]
    result = {}
    for c in range(1, 28):
        for t, b2, b4 in guarded_phases:
            d2, d4 = c + b2, 2 * c + b4 - 4
            failing = tuple(i for i, d in ((2, d2), (4, d4)) if d < -1)
            if failing:
                result[c] = t, failing, d2, d4
                break
        assert c in result
    return result


def check_literal(base, trace, failures):
    histogram = Counter()
    for j in range(27):
        c0 = base // 400
        c = c0 + (j - (c0 - 28)) % 27
        v = base // 95
        u0 = base // 96
        # The integer sum is B = 9A + 4c + 2u + 2v + 3.
        u = u0 + (5 * (base - 4 * c - 2 * v - 3) - u0) % 9
        anchor, remainder = divmod(base - 4 * c - 2 * u - 2 * v - 3, 9)
        assert remainder == 0
        assert Fraction(1, 500) <= Fraction(c, base) <= Fraction(1, 250)
        assert all(Fraction(1, 100) <= Fraction(w, base) <= Fraction(1, 90) for w in (u, v))
        complete = 1 + (c - 28) // 27
        section_c = c - 27 * complete
        assert (j, 2 * j) == ((c - 28) % 27, (2 * c - 56) % 54)
        partial, expected_failure, d2, d4 = failures[section_c]
        stop = BLOCK * complete + partial
        leading = (anchor, anchor, anchor, anchor + c + u + v + 4, anchor,
                   anchor + c + u, anchor + c + v, anchor + c, anchor)
        current = translate(leading, START_OFFSET, 1)
        for t in range(stop + 1):
            block, phase = divmod(t, BLOCK)
            expected_state = translate(leading, DRIFT, block)
            for _ in range(phase % 4):
                expected_state = rotate(expected_state)
            expected_state = translate(expected_state, trace[phase], 1)
            assert current == expected_state
            assert 10 * min(current) > base + 1
            failing = tuple(i for i in (2, 3, 4) if margin(current, i) < -1)
            successor = literal_step(current, base)
            assert successor == paired_step(current, secondary=(t == stop))
            if t < stop:
                assert not failing
            else:
                assert failing == expected_failure
                assert (margin(current, 2), margin(current, 4)) == (d2, d4)
                histogram[failing] += 1
            current = successor
    return histogram


def main():
    trace = offset_trace()
    failures = first_failures(trace)
    histogram = Counter(mark for _, mark, _, _ in failures.values())
    assert histogram == {(2,): 19, (4,): 5, (2, 4): 3}
    print('joint support (j, 2j), 0 <= j < 27; index 27')
    for section, (t, mark, d2, d4) in failures.items():
        print(f'section c={section:2} phase={t:3} walls={mark!s:6} d2={d2:3} d4={d4:3}')
    print('joint weights:', {k: str(Fraction(v, 27)) for k, v in sorted(histogram.items())})
    for base in (100003, 1000003):
        literal = check_literal(base, trace, failures)
        assert literal == histogram
        print('literal base', base, dict(sorted(literal.items())), 'including failure updates', flush=True)


if __name__ == '__main__':
    main()
