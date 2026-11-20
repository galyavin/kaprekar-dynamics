"""Independent whole-word replays of exact affine and centered predictions."""
import argparse
import fcntl
import json
from pathlib import Path
import sys

from certify import predict


def differences_from_sorted(ordered):
    return tuple(ordered[-1-i] - ordered[i] for i in range(5))


def subtract(base, ordered):
    upper = lower = 0
    for lo, hi in zip(ordered, reversed(ordered)):
        upper = upper*base + hi
        lower = lower*base + lo
    value = upper - lower
    out = []
    for _ in range(11):
        value, digit = divmod(value, base)
        out.append(digit)
    assert value == 0
    return tuple(reversed(out))


def failed_guards(base, q):
    if q[-1] == 0:
        raise AssertionError("A positive-inner path hit zero")
    pairs = tuple((q[i] - (i == 4), base - q[4-i] - (i < 4)) for i in range(5))
    return tuple(min(pairs[i])-max(pairs[i+1]) for i in range(4))


def in_source(base, q):
    return min(base-q[0], q[0]-q[1], q[1]-q[2], q[2]-q[3],
               q[3]-q[4], q[4]) >= (base+7)//7 + 64


def check_translation(base, q, slack):
    predicted = differences_from_sorted(sorted(subtract(base, sorted((*q, *[0]*6)))))
    for j in range(4):
        increment = tuple(j*i for i in (5, 4, 3, 2, 1))
        translated = tuple(a+b for a, b in zip(q, increment))
        got = differences_from_sorted(sorted(
            subtract(base+6*j, sorted((*translated, *[0]*6)))
        ))
        expect = tuple(a+b - (min(j, -slack) if i == 3 else 0)
                       for i, (a, b) in enumerate(zip(predicted, increment)))
        assert got == expect, (base, q, slack, j, got, expect)


def replay(base, source, predicted):
    origin = (*source, *[0]*6)
    word = origin
    seen_words = {}
    seen_differences = {}
    first_failure = last_failure = source_time = difference_time = cyclic_q = None
    single_failure_phases = 0
    t = 0
    while word not in seen_words:
        seen_words[word] = t
        order = sorted(word)
        q = differences_from_sorted(order)
        if q in seen_differences and difference_time is None:
            difference_time = seen_differences[q]
            cyclic_q = q
        else:
            seen_differences.setdefault(q, t)
        slacks = failed_guards(base, q)
        if min(slacks) < 0:
            assert slacks[0] >= 0 and slacks[3] >= 0
            assert slacks[1] == slacks[2] and slacks[1] in (-1, -2)
            single_failure_phases += 1
            last_failure = t
            if first_failure is None:
                first_failure = t
                check_translation(base, q, slacks[1])
        if in_source(base, q) and source_time is None:
            source_time = t
            assert q == predicted['source_state']
        if t == predicted['second_block_start']:
            assert q == predicted['second_block_state']
        if t == predicted['bridge_start']:
            assert q == predicted['bridge_state']
        word = subtract(base, order)
        t += 1
    observed = {
        'first_failure': first_failure,
        'last_failure': last_failure,
        'source_time': source_time,
        'difference_time': difference_time,
        'word_time': seen_words[word],
        'period': t-seen_words[word],
        'first_cyclic_difference': cyclic_q,
        'first_cyclic_word': word,
        'first_cyclic_words_difference': differences_from_sorted(sorted(word)),
    }
    assert all(observed[k] == predicted[k] for k in observed), (base, source, observed, predicted)
    return {
        'base': base,
        'initial_word': origin,
        'whole_word_steps': t,
        'failed_difference_phases': single_failure_phases,
        'translation_controls': 4,
        **predicted,
    }


def construct_source(base, z, r1, r2, row):
    b = 73810*base//100000
    n = (base//2000 + 3)//6
    m = (base//400 + 5*n - 1 + z)//2
    e = b - 3*n - (base-3)//2
    approx_c = 57977*base//100000
    approx_a = 94085*base//100000
    u = approx_a + approx_c - 8*n - 10*m - 3*(base-1)//2
    count = (u-r1)//22
    r = b - 4*n - 2*m - 4*count
    p = approx_c - 3*n - 2*count
    twice_offset = r2 - (2*base + 2*p - 4*r + row['second_clock_constant'])
    assert twice_offset % 2 == 0
    p += (twice_offset//2 + 18) % 36 - 18
    c = p + 3*n + 2*count
    a = 3*(base-1)//2 + r1 + 22*count - c + 8*n + 10*m
    d = base-c-1+5*n+z-2*m
    return a, b, c, d, e


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("compute_lock", type=Path)
    args = parser.parse_args()
    rows = json.loads(args.certificate.read_text())['residues']
    controls = (
        (1, -4), (1, 0), (1, 7), (1, 11), (1, 17),
        (2, -4), (2, 0), (2, 6), (2, 12), (2, 17),
    )
    results = []
    with args.compute_lock.open('r+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for index, (z, r1) in enumerate(controls):
            base = 1000001 + 6*(index % 2)
            for r2 in (601, 671):
                source = construct_source(base, z, r1, r2, rows[f'{z},{r1}'])
                prediction = predict(base, source)
                assert prediction['first_overshoot'] == z
                assert prediction['first_remainder'] == r1
                assert prediction['second_remainder'] == r2
                results.append(replay(base, source, prediction))
    total = sum(row['whole_word_steps'] for row in results)
    print(json.dumps({'python': sys.version.split()[0], 'controls': len(results),
                      'whole_word_steps': total, 'results': results}, indent=2))


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run verification with assertions enabled")
    main()
