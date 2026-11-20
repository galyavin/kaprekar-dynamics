import json

import channel as kap

B, Q = kap.BASE, kap.SOURCE
VEC = kap.V
AFF_MIN = kap.BASE_MIN
S, U, V, W = kap.S, kap.U, kap.FAILURE, kap.W
C = "5067128934a"
ORDERS = dict(S=S, U=U, V=V, W=W, C=C)
ALIGN_DELTA = {0: (-5, 0, 1, 2, -1), 2: (-5, 2, 3, 4, -1)}
PREP_DELTA = (-10,-2,0,2,-2)
LATE_DELTA = (-20,-4,-2,2,-4)
PAIRED_DELTA = (-20,-4,0,4,-4)
R1 = range(-4, 18)


def constant(row):
    x = row.coefficients[-1]
    assert row == kap.Affine.constant(x), row
    assert x.denominator == 1, x
    return int(x)


def h(base, q):
    return 2*(q[1]-q[4])-base+1


def g(base, q):
    return base-q[2]-q[3]-1


def rows_are_nonnegative(row):
    if isinstance(row, kap.Affine):
        slope = row.normalized_minimum()
        return slope >= 0 and slope * AFF_MIN + row.coefficients[-1] >= 0
    return row >= 0


def pair_rows(base, q):
    d = kap.digits(base, q)
    for i in range(4):
        for a in (i, i+6):
            for b in (i+1, i+7):
                yield d[a]-d[b]


def is_paired(base, q):
    rows = tuple(pair_rows(base, q))
    if all(map(rows_are_nonnegative, rows)):
        return True
    # For each failing phase certify an integer guard uniformly <= -1.
    assert any(rows_are_nonnegative(-1-r) for r in rows)
    return False


def word(s):
    return tuple(ORDERS[x] for x in s)


def encode(orders):
    inv = {w:k for k,w in ORDERS.items()}
    return ''.join(inv[w] for w in orders)


def repeat_affine(q, w, n, delta, failures, label):
    kap.nonnegative(n-1, ('positive repetitions', label))
    failure_rows = {}
    for index, start in enumerate((q, kap.shift(q,n-1,delta))):
        t = start
        for i, order in enumerate(w):
            assert order[0] == '5' and order[-1] == 'a'
            if i in failures:
                negative = {j for j,row in enumerate(pair_rows(B,t))
                            if rows_are_nonnegative(-1-row)}
                if index == 0:
                    failure_rows[i] = negative
                else:
                    negative &= failure_rows[i]
                assert negative, (label,index,i)
            else:
                assert is_paired(B,t), (label,index,i)
            t = kap.walk(B,t,(order,), (label,index,i))
        assert t == kap.shift(start,1,delta), (label, index, t)
    return kap.shift(q,n,delta)


def first_failure(cohort, z, source=Q, base=B):
    exit_h = {0:-4,2:-6}[cohort]
    affine = isinstance(base,kap.Affine)
    n = (h(base,source)-exit_h)/6 if affine else (h(base,source)-exit_h)//6
    q = kap.shift(source,n,tuple(-x for x in VEC))
    align = word('SUSSS')
    if affine:
        q = repeat_affine(q, align, 1, ALIGN_DELTA[cohort], (),
                          ('alignment',cohort))
    else:
        q = kap.walk(base,q,align,'numeric alignment')
    m = (g(base,q)+z)/2 if affine else (g(base,q)+z)//2
    prep = word({0:'SUSSSSUSSS',2:'SSSSSSUSSS'}[cohort])
    if affine:
        q = repeat_affine(q, prep, m, PREP_DELTA, (), ('prep',cohort,z))
    else:
        q = kap.shift(q,m,PREP_DELTA)
    pinned_h, pinned_g = h(base,q), g(base,q)
    if affine:
        pinned_h, pinned_g = constant(pinned_h), constant(pinned_g)
    assert (pinned_h, pinned_g) == ({0:-2,2:0}[cohort], -z)
    assert not is_paired(base,q)
    return q, 5*n+5+10*m


def late_word(cohort):
    return word({0:'VUSSSSUSSS',2:'VSSSSSUSSS'}[cohort])*2


def discover(q):
    """Propose a bounded switch into a paired return with Task 029 displacement."""
    words = []
    states = [q]
    pairings = []
    for _ in range(96):
        # These orders still have to pass exact guards at their points of use.
        q = states[-1]
        pairings.append(is_paired(B,q))
        next_order = kap.known_sorted(kap.digits(B,q))
        words.append(next_order)
        states.append(kap.apply_order(B,q,next_order))
    for i in range(len(words)-40+1):
        a,b = words[i:i+20], words[i+20:i+40]
        if (a == b and all(pairings[i:i+40]) and V not in a
                and states[i+20] == kap.shift(states[i],1,PAIRED_DELTA)):
            return tuple(words[:i]), tuple(a), tuple(j for j in range(i) if not pairings[j])
    raise AssertionError(('no paired return', encode(words)))


def certify():
    leaves = []
    for cohort in (0,2):
        for z in (1,2):
            q, t = first_failure(cohort,z)
            for r1 in R1:
                n = (kap.first_clock(B,q)-r1)/22
                end = repeat_affine(q,late_word(cohort),n,LATE_DELTA,
                                    (0,10) if z == 1 else (0,),
                                    ('late',cohort,z,r1))
                assert constant(kap.first_clock(B,end)) == r1
                switch, paired, failures = discover(end)
                t2 = t + 20*n + len(switch)
                for i, order in enumerate(switch):
                    assert order[0] == '5' and order[-1] == 'a'
                    assert is_paired(B,end) == (i not in failures)
                    end = kap.walk(B,end,(order,),('switch',cohort,z,r1,i))
                assert end[0] == Q[0]-t2
                clock = kap.second_clock(B,end)
                for r2 in (601, 671):
                    loops = (clock-r2)/72
                    bridge = repeat_affine(end,paired,loops,PAIRED_DELTA,(),
                                           ('paired',cohort,z,r1,r2))
                    assert constant(kap.second_clock(B,bridge)) == r2
                    assert bridge[0] == Q[0] - (t2+20*loops)
                    kap.check_bridge(B,bridge,('bridge',cohort,z,r1,r2))
                leaves.append(dict(cohort=cohort,z=z,r1=r1,switch=encode(switch),
                                   failures=failures,paired=encode(paired)))
    return leaves


if __name__ == '__main__':
    print(json.dumps(certify(),indent=2))
