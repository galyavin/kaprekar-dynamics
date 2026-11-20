"""Reconstruct the archived root cones and strict-profile physical lifts."""
import json
import sys
from math import lcm
from pathlib import Path

from common import (add, complement_period, cycles, difference_rows, insertions, lift, lift_state,
                    literal, physical_period)


def records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


PROFILES = {
    31: {},
    33: {},
    41: {},
    43: {3: 3, 7: 2},
    45: {3: 1, 5: 1, 9: 2, 15: 2},
    51: {3: 2, 5: 1, 9: 3, 17: 0},
    57: {3: 1, 5: 1, 19: -1},
    63: {},
    65: {},
    91: {},
}


def verify(repo):
    evidence = repo/'kaprekar-investigation/cross_length'
    astra = repo/'docs/drafts/research-direction/astra-strategy-16'
    refs = records(evidence/'lc_classification/root-certificates.jsonl')
    refs += records(evidence/'lc2_order12/certificate.jsonl')
    refs.sort(key=lambda record: record['order'])
    assert [ref['order'] for ref in refs] == list(range(3, 14))
    routes = {}
    for ref in refs:
        start, end = ref['cycle']
        routes[ref['order']] = insertions(ref['source'], ref['insertions'], end, start)

    also = records(astra/'period132-root.jsonl')[0]
    start, end = also['cycle']
    special = insertions(also['source'], also['insertions'], end, start)
    standard = routes[12]
    assert special.rays == standard.rays
    assert add(also['source'], refs[12-3]['source'], -1) == (0,0,0,0,0,0,12592,0,0,0,0,0)

    # Order r implies 2**r is either +1 or -1 modulo n. For each
    # candidate, either exhibit a rational strict direction or a zero
    # row on the entire invariant, root-equal profile subspace.
    screen = []
    for r in routes:
        route = routes[r]
        ref = refs[r-3]
        strict = []
        excluded = []
        for n in range(5,2**r+2,2):
            if (2**r-1)%n and (2**r+1)%n:
                continue
            _, labels, *parts = cycles(n)
            if len(labels)!=r:
                continue
            overrides = PROFILES.get(n,{})
            m = (n+1)//2
            rows = [difference_rows(tuple(int(i in part) for i in range(m)))[1:]
                    for part in parts]
            if not parts:
                excluded.append((n,'full'))
                continue
            zero_rows = [j+2 for j,values in enumerate(zip(*rows)) if not any(values)]
            if zero_rows:
                excluded.append((n,zero_rows))
                continue
            assert all(sum(weight*value for weight,value in zip(
                (overrides.get(part[0],1) for part in parts),values))>0 for values in zip(*rows))
            screen.append((n,ref,route,overrides))
            strict.append(n)
        print(f'Complete screen: {r=} {strict=} {excluded=}',flush=True)

    for n,ref,route,overrides in screen:
        r = ref['order']
        mean_scaling = n in (45,91)
        labels, complement, profile, D, delta = lift(n, ref['source'], route, overrides, mean_scaling)
        if n not in PROFILES:
            continue
        # Pick generic offsets, including different residues within each
        # complementary orbit. Check the prefix literally, then test the cycle.
        B = D*1000+7
        free = sum(map(len, complement))
        source = lift_state(B, tuple(ref['source']), labels, complement, profile, tuple(range(free)))
        literal(source, 10)
        rperiod = route.period
        expected = lcm(complement_period(complement), rperiod)
        # Take distinct offsets to realize the full complementary period.
        terminal = lift_state(B, route.terminal, labels, complement, profile,
                              tuple(i*21 for i in range(free)))
        literal(terminal, expected)
        assert physical_period(terminal, expected) == expected
        print(f'  literal: 10 source transitions and terminal {expected} at B={B}', flush=True)
        if r == 12:
            assert complement_period(complement) == 12
            labels, complement, pi, d, _ = lift(n, also['source'], special, overrides, True)
            B = d*1000+9
            terminal = lift_state(B, special.terminal, labels, complement, pi,
                                  tuple(i*21 for i in range(free)))
            assert len({terminal[j] for part in complement for j in part}) == free
            literal(terminal, 132)
            assert physical_period(terminal, 132) == 132
            print(f'  period-132 literal: {B=}; positive limiting lower bound with same delta, divisor {d}', flush=True)


if __name__ == '__main__':
    verify(Path(sys.argv[1]))
