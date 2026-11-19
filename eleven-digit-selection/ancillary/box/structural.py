"""Use the half-return clock to put a twelve-shift pair on one affine section."""
from collections import Counter
from fractions import Fraction as Q
import json

import cohort1 as orbit
from sections import delta_at, prefixes, q3


def primitive_rate(drift, word):
    return tuple(Q(d,len(word)) for d in drift)


def sections():
    for label, guard, switch, block, failed, q2, j2, kappa in prefixes():
        assert len(block) == 20 and tuple(block[:10]) == tuple(block[10:])
        assert all(order[0] == '5' and order[-1] == 'a' for order in (*switch, *block))
        # The ten-step displacement is affine in r2. Verify both extremes of
        # the whole old section interval, plus the first section at q2.
        half_drift = tuple(x//2 for x in orbit.SECOND_DRIFT)
        for rows in (q2, q3(q2,j2,601), q3(q2,j2,671)):
            assert orbit.step_many(orbit.B,rows,block[:10]) == orbit.translate(rows,1,half_drift)
        for r2 in range(601,637,2):
            # Pass from the upper representative to the common lower one.
            upper = q3(q2,j2,r2+36)
            lower = q3(q2,j2,r2)
            assert guard.walk(upper,block[:10]) == lower
            rows, word, drift, delay = orbit.checkpoint_from(lower, block, guard)
            assert all(order[0] == '5' and order[-1] == 'a' for order in word)
            yield (*label,r2), guard, rows, word, drift, delay


def report():
    totals = Counter()
    slope_types = {}
    solution = None
    for label, _, rows, word, drift, _ in sections():
        rate = primitive_rate(drift,word)
        assert rate[0] == -1 and len(word) in (10,20)
        linear = tuple(tuple(Q(x,orbit.SCALE) for x in row.terms[:3]) for row in rows)
        if solution is None:
            a,b,c = linear[1]
            # Centered second coordinate equals the oriented drift times the
            # number of steps needed to align the anchors.
            coefficient = c+rate[1]*Q(10,9)
            intercept = 12*a+18*b-8-rate[1]*Q(80,3)
            solution = -intercept/coefficient
            assert solution == -12
        y = int(solution)
        delta_j = 60-4*y
        assert delta_j % 36 == 0
        checkpoint_clock = 10*delta_j//36
        delay = checkpoint_clock+10
        assert delay % len(word) == 0 and delay//len(word) <= 10
        lift = tuple(delta_at(row,(12,18,y))-2*v for row,v in zip(rows,orbit.V))
        assert lift == tuple(d*delay//len(word) for d in drift), label
        key = (len(word),rate,linear)
        totals[key] += 1
        slope_types[key] = lift
    print(json.dumps({'canonical_sections': sum(totals.values()), 'change_in_p_r': [18,int(solution)],
                      'change_in_J': delta_j, 'checkpoint_clock_change': checkpoint_clock,
                      'source_delay': delay,
                      'slope_types': [{'period': period, 'per_step': list(map(str,rate)),
                                       'coefficients': [[str(x) for x in row] for row in linear],
                                       'lift': slope_types[period,rate,linear], 'sections': n}
                                      for (period,rate,linear),n in sorted(totals.items())]}, indent=2))


if __name__ == '__main__':
    report()
