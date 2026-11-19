"""Finite whole-word controls for the exact anchored clocks."""

from fractions import Fraction as F
import json

import certificate as cert
import channel as kap


def route_prediction(base, source, cohort, routes):
    exit_h = {0: -4, 2: -6}[cohort]
    n, rem = divmod(cert.h(base, source)-exit_h, 6)
    assert rem == 0
    alignment = kap.shift(source, n, tuple(-x for x in kap.V))
    alignment = kap.shift(alignment, 1, cert.ALIGN_DELTA[cohort])
    z = 2 - cert.g(base, alignment) % 2
    qF, first = cert.first_failure(cohort, z, source, base)
    clock = int(kap.first_clock(base, qF))
    n_late, r1 = divmod(clock+4, 22)
    r1 -= 4
    qR = kap.shift(qF, n_late, cert.LATE_DELTA)
    route = routes[cohort, z, r1]
    q = qR
    for order in cert.word(route['switch']):
        q = kap.apply_order(base, q, order)
    t2 = first+20*n_late+len(route['switch'])
    last = first+20*n_late+max(route['failures'])
    n_paired, r2 = divmod(kap.second_clock(base,q)-601, 72)
    r2 += 601
    assert r2 in kap.R2_RANGE
    q3 = kap.shift(q,n_paired,cert.PAIRED_DELTA)
    t3 = t2+20*n_paired
    assert q3[0] == source[0]-t3
    entry_gap = (base+7)//7+64
    bridge = entry_gap-(base-q3[0])
    x = kap.centered_from(base,q3)
    for _ in range(bridge):
        x = kap.centered_step(x)
    assert min(kap.full_gaps(base,kap.physical_from(base,x))) == entry_gap
    delay, period, cyclic = kap.return_to_cycle(x)
    arrival = t3+bridge+delay+1
    return dict(
        first_failure=first, last_failure=last, permanent_entry=t3+bridge,
        first_cyclic_word_time=arrival, period=period,
        cyclic_word=kap.digits(base,kap.physical_from(base,cyclic)),
        z=z, r1=r1, r2=r2, paired_block_start=t2, bridge_start=t3, bridge_state=q3,
    )


def whole_word_replay(base, source):
    original = (*source, 0, 0, 0, 0, 0, 0)
    word = original
    seen = {}
    first = last = None
    entry = None
    gap = (base+7)//7+64
    t = 0
    while word not in seen:
        seen[word] = t
        ascending = sorted(word)
        q = tuple(ascending[10-i]-ascending[i] for i in range(5))
        assert q[-1] > 0
        if entry is None and min(kap.full_gaps(base, q)) >= gap:
            entry = t
        if not cert.is_paired(base, q):
            first = t if first is None else first
            last = t
        # Independent integer subtraction and radix decoding.
        hi = lo = 0
        for a,b in zip(ascending,reversed(ascending)):
            hi = base*hi+b
            lo = base*lo+a
        value = hi-lo
        digits = []
        for _ in ascending:
            value,digit = divmod(value,base)
            digits.append(digit)
        assert value == 0
        word = tuple(reversed(digits))
        t += 1
    return dict(first_failure=first, last_failure=last, permanent_entry=entry,
                first_cyclic_word_time=seen[word], period=t-seen[word],cyclic_word=word)


def choose_input(base, cohort, desired_z):
    mids = [int((lo+hi)*base/2) for lo,hi in kap.BOX]
    a,b,c,h,g = mids
    separation = (base+h)//2
    separation += (cohort-separation) % 3
    q = (a,b,c,base-c-g,b-separation)
    n = (cert.h(base,q)-{0:-4,2:-6}[cohort])//6
    aligned = kap.shift(q, n, tuple(-x for x in kap.V))
    aligned = kap.shift(aligned, 1, cert.ALIGN_DELTA[cohort])
    if 2-cert.g(base,aligned)%2 != desired_z:
        q = (*q[:3],q[3]+1,q[4])
    chart = (a,b,c,2*(b-q[4])-base,base-c-q[3])
    assert all(lo*base <= value <= hi*base for (lo,hi),value in zip(kap.BOX,chart))
    return q


def main():
    with open('routes.json') as f:
        routes = {(r['cohort'],r['z'],r['r1']): r for r in json.load(f)}
    outcomes = []
    for base, z in ((1_000_001,1),(1_000_007,2)):
        for cohort in (0,2):
            source = choose_input(base, cohort, z)
            expected = route_prediction(base, source, cohort, routes)
            assert expected['z'] == z
            observed = whole_word_replay(base, source)
            assert all(expected[key] == value for key,value in observed.items())
            t = observed['first_cyclic_word_time']
            ramp = max(0, min(1, 2000*(F(t,base)-F(97,800))))
            outcomes.append(dict(base=base,cohort=cohort,source=source,ramp=str(ramp),
                                 **expected))
    print(json.dumps(outcomes,indent=2))


if __name__ == '__main__':
    main()
