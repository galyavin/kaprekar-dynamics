"""Match the cohort-2 quotient and the remaining weighted whole-box labels."""

from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path

import certificate as cert
import channel as kap
import cohort1 as orbit
import control


def rows_at_switch(cohort, z, r1):
    base, p, r = orbit.B, orbit.P, orbit.R
    pin = {1: -2, 2: 0}[cohort]
    return (3 * (base - 1) / 2 + r1 - p, r, p, base - 1 + z - p,
            r - (base + pin - 1) / 2)


class Cohort2Guards(orbit.BoxGuards):
    """Same leading chart, with the actual cohort-2 intercepts."""

    def __init__(self, first_failure, z, r1):
        count = (kap.first_clock(cert.B, first_failure) - r1) / 22
        state = kap.shift(first_failure, count, cert.LATE_DELTA)
        r, p = state[1:3]
        reference_p, reference_r, _ = orbit.prefix_coordinates(
            cert.B, kap.A, kap.B2, kap.C, kap.H0, kap.G0, z, r1
        )
        assert p.coefficients[:-1] == reference_p.coefficients[:-1]
        assert r.coefficients[:-1] == reference_r.coefficients[:-1]
        self.p0, self.r0 = p.coefficients[-1], r.coefficients[-1]


def build_section(cohort, label, switch, block, failures, guards):
    z, r1, r2 = label
    q2 = orbit.step_many(orbit.B, rows_at_switch(cohort, z, r1), switch)
    clock = kap.second_clock(orbit.B, q2)
    kappa = orbit.integer(clock.eval(0, 0, 0))
    assert kappa % 2 == 1
    assert clock == 2 * orbit.B + 2 * orbit.P - 4 * orbit.R + kappa
    q3 = orbit.translate(q2, (clock - r2) / 72, orbit.SECOND_DRIFT)
    checkpoint, word, delta, from_q3 = orbit.checkpoint_from(q3, block, guards)
    base, p, r, _ = zip(*(row.terms for row in checkpoint))
    columns = (
        tuple(orbit.integer(F(36 * x, orbit.SCALE)) for x in p),
        tuple(orbit.integer(F(2 * x + y, orbit.SCALE)) for x, y in zip(p, r)),
    )
    lift = tuple(orbit.integer(F(36*x, orbit.SCALE) - 6*v) for x,v in zip(base, orbit.V))
    pin = r1 - 11 * (max(failures) // 10)
    return orbit.Section(label, word, checkpoint, delta, columns, lift, kappa,
                         from_q3, pin)


def build_sections(routes, accepted, records):
    first_failures = {z: cert.first_failure(2, z)[0] for z in (1, 2)}
    two = {}
    for record in routes:
        if record['cohort'] != 2:
            continue
        z, r1 = record['z'], record['r1']
        guard = Cohort2Guards(first_failures[z], z, r1)
        for r2 in kap.R2_RANGE:
            label = z, r1, r2
            two[label] = build_section(
                2, label, cert.word(record['switch']), cert.word(record['paired']),
                record['failures'], guard
            )

    one = {}
    residual = {
        tuple(label)
        for row in accepted["matches"]["6"] if row["K_mod_6"] % 2 == 0
        for side in ("left_residual", "right_residual") for label in row[side]
    }
    for label in sorted(residual):
        z, r1, _ = label
        record = records[f'{z},{r1}']
        one[label] = build_section(
            1, label, record['orders_to_first_paired_block'], record['paired_block'],
            record['failed_times_in_this_prefix'], orbit.BoxGuards(z, r1)
        )
    return one, two


def compute_group_drifts(sections):
    drifts = {}
    for section in sections:
        delta, normal = orbit.decompose(section.displacement, section.columns)
        assert not any(normal), (section.label, normal)
        assert orbit.lattice_fraction(orbit.component(11, delta)) == (0, 0)
        assert len({orbit.lattice_fraction(orbit.component(j, delta))
                    for j in range(11)}) == 11
        assert drifts.setdefault(section.group, delta) == delta
        for k in range(7):
            assert orbit.difference(section.origin(k+6), section.origin(k)) == section.lift
    return drifts


def verify_witness(a, b, k):
    s, t, j = orbit.witness(a, b, k, 1)
    assert 0 <= j <= 10 and a.group == b.group
    assert all(
        right - left == s * c1 + t * c2 + j * delta
        for left, right, c1, c2, delta in
        zip(a.origin(k), b.origin(k+1), *a.columns, a.displacement)
    )
    return s, t, j


def match(left, right, k, drifts):
    a = orbit.keys_for(left, k, drifts)
    b = orbit.keys_for(right, k+1, drifts)
    pairs, unpaired_left, unpaired_right = [], [], []
    for key in sorted(a.keys() | b.keys()):
        first, second = a.get(key, ()), b.get(key, ())
        for x, y in zip(first, second):
            pairs.append((x, y, *verify_witness(left[x], right[y], k)))
        unpaired_left.extend(first[len(second):])
        unpaired_right.extend(second[len(first):])
    assert len(pairs) + len(unpaired_left) == len(left)
    assert len(pairs) + len(unpaired_right) == len(right)
    return pairs, sorted(unpaired_left), sorted(unpaired_right)


def compare(one, two, accepted):
    assert len(two) == 1584
    drifts = compute_group_drifts((*one.values(), *two.values()))
    summary = {}
    replay_witnesses = []
    for k in (0, 2, 4):
        pairs, left2, right2 = match(two, two, k, drifts)
        residuals = accepted['matches']['6'][k]
        assert residuals['K_mod_6'] == k and residuals['matched'] == 1520
        left = {(c, *label): one[tuple(label)]
                for label in residuals['left_residual'] for c in (0, 1)}
        right = {(c, *label): one[tuple(label)]
                 for label in residuals['right_residual'] for c in (0, 1)}
        left.update({(2, *label): two[label] for label in left2})
        right.update({(2, *label): two[label] for label in right2})
        joint, pending_left, pending_right = match(left, right, k, drifts)
        assert len(pending_left) == len(pending_right)
        assert len(pending_left) + 2 * 1520 + len(pairs) + len(joint) == 3 * 1584
        summary[str(k)] = {
            'cohort2_matched': len(pairs),
            'cohort2_return_counts': sorted(Counter(row[-1] for row in pairs).items()),
            'cohort2_left_residual': left2, 'cohort2_right_residual': right2,
            'additional_whole_box_matches': len(joint),
            'additional_matches_by_cohort': sorted(Counter((a[0], b[0])
                                                          for a,b,*_ in joint).items()),
            'whole_box_left_residual': pending_left,
            'whole_box_right_residual': pending_right,
            'whole_box_residual_mass_per_arm': str(F(len(pending_left), 3*1584)),
        }
        a,b,s,t,j = max(pairs, key=lambda pair: pair[-1])
        replay_witnesses.append((k, (2,*a), (2,*b), s,t,j))
        if joint:
            replay_witnesses.append((k, *joint[0]))
    return summary, replay_witnesses


def source_from(base, cohort, section, h, r, N, n, m):
    z, r1, r2 = section.label
    p0 = ((r2-section.clock_constant)//2-base) % 36
    p = p0+36*h+2*r
    pin = 0 if cohort == 2 else -2
    qR = (3*(base-1)//2+r1-p, r, p, base-1+z-p, r-(base+pin-1)//2)
    qF = kap.shift(qR,-N,cert.LATE_DELTA)
    qA = kap.shift(qF,-m,cert.PREP_DELTA)
    qE = kap.shift(qA,-1,cert.ALIGN_DELTA.get(cohort,(0,0,0,0,0)))
    source = kap.shift(qE,n,kap.V)
    assert orbit.in_box(base,source) and (source[1]-source[4])%3 == cohort
    return source


def physical_observation(base, cohort, source, section, routes, K, h, r):
    if cohort == 1:
        pred = kap.predict(base,source)
        label = pred['first_overshoot'],pred['first_remainder'],pred['second_remainder']
        word_time, cyclic_word = pred['word_time'],pred['first_cyclic_word']
    else:
        pred = control.route_prediction(base,source,cohort,routes)
        label = pred['z'],pred['r1'],pred['r2']
        word_time, cyclic_word = pred['first_cyclic_word_time'],pred['cyclic_word']
    assert label == section.label, (label,section.label)

    q = pred['bridge_state']
    for _ in range(section.from_q3):
        assert cert.is_paired(base,q)
        order = kap.known_sorted(kap.digits(base,q))
        q = kap.apply_order(base,q,order)
    predicted_q = kap.shift(section.origin(K),h,section.columns[0])
    predicted_q = kap.shift(predicted_q,r,section.columns[1])
    predicted_q = kap.shift(predicted_q,K,orbit.V)
    assert q == predicted_q
    X = kap.centered_from(base,q)
    cyclic_q = (*cyclic_word[:4],cyclic_word[4]+1)
    t4 = pred['bridge_start']+section.from_q3
    return X, t4, word_time, kap.centered_from(base,cyclic_q), pred['period']


def replay(witness, one, two, routes):
    k,a,b,s,t,j = witness
    cohort_a, cohort_b = a[0], b[0]
    sec_a = ({0:one,1:one,2:two}[cohort_a])[a[1:]]
    sec_b = ({0:one,1:one,2:two}[cohort_b])[b[1:]]
    K = 400002+k
    base = 6*K+5
    ref_cohort = 2 if cohort_a == 2 else 0
    z,r1,r2 = sec_a.label
    reference = control.choose_input(base,ref_cohort,z)
    align = kap.shift(reference,(cert.h(base,reference)-{0:-4,2:-6}[ref_cohort])//6,
                      tuple(-x for x in orbit.V))
    align = kap.shift(align,1,cert.ALIGN_DELTA[ref_cohort])
    n = (cert.h(base,reference)-{0:-4,2:-6}[ref_cohort])//6
    m = (cert.g(base,align)+z)//2
    qF = kap.shift(align,m,cert.PREP_DELTA)
    N = (int(kap.first_clock(base,qF))-r1)//22
    qR = kap.shift(qF,N,cert.LATE_DELTA)
    r=qR[1]
    p0 = ((r2-sec_a.clock_constant)//2-base)%36
    h=(qR[2]-p0-2*r)//36

    source_a = source_from(base,cohort_a,sec_a,h,r,N,n,m)
    source_b = source_from(base+6,cohort_b,sec_b,h-s,r-t,N,n,m)
    obs_a = physical_observation(base,cohort_a,source_a,sec_a,routes,K,h,r)
    obs_b = physical_observation(base+6,cohort_b,source_b,sec_b,routes,K+1,h-s,r-t)
    X,t4a,Ta,terminal_a,period_a = obs_a
    Xb,t4b,Tb,terminal_b,period_b = obs_b
    delay=j*len(sec_a.word)
    for _ in range(delay):
        X=kap.centered_step(X)
    assert X == Xb
    assert (terminal_a,period_a) == (terminal_b,period_b)
    assert (Ta-t4a)-(Tb-t4b) == delay
    return dict(K_mod_6=k,bases=(base,base+6), labels=(a,b),
                sources=(source_a,source_b), time_to_checkpoint=(t4a,t4b),
                time_to_first_cyclic_word=(Ta,Tb), time_to_align=delay,
                first_cyclic_centered_difference=terminal_a,period=period_a)


def load(name):
    return json.loads(Path(name).read_text())


def main():
    routes = load('routes.json')
    accepted = load('accepted_quotient.json')
    one, two = build_sections(routes, accepted, load('endpoint_verification.json')['residues'])
    summary, witnesses = compare(one, two, accepted)
    print(json.dumps({
        'cohort2_sections':len(two),
        'cohort1_residual_sections':len(one),
        'cohort2_checkpoint_word_counts':sorted(Counter(
            orbit.display(section.word) for section in two.values()
        ).items()),
        'by_K_mod_6':summary,
        'replay_witnesses': [
            replay(witness, one, two, {(r['cohort'],r['z'],r['r1']):r for r in routes})
            for witness in witnesses
        ],
    }, indent=2))


if __name__ == '__main__':
    main()
