"""Prefix rows for both independently guarded cohorts."""
from fractions import Fraction as Q
import json
from pathlib import Path

import certificate as prefix
import cohort1 as orbit
import quotient as box


def prefixes():
    one = json.loads(Path('endpoint_verification.json').read_text())['residues']
    two = {(row['z'], row['r1']): row for row in json.loads(Path('routes.json').read_text())
           if row['cohort'] == 2}
    for cohort in (1, 2):
        for z in (1, 2):
            for r1 in range(-4, 18):
                if cohort == 1:
                    row = one[f'{z},{r1}']
                    switch, block = row['orders_to_first_paired_block'], row['paired_block']
                    failed = row['failed_times_in_this_prefix']
                    guard = orbit.BoxGuards(z, r1)
                else:
                    row = two[z, r1]
                    switch, block = prefix.word(row['switch']), prefix.word(row['paired'])
                    failed = row['failures']
                    guard = box.Cohort2Guards(prefix.first_failure(2, z)[0], z, r1)
                q2 = orbit.step_many(orbit.B, box.rows_at_switch(cohort, z, r1), switch)
                a, b, _, d, e = q2
                j = 2*a + 2*b - 2*d + 4*e - 3*orbit.B
                kappa = orbit.integer(j.eval(0, 0, 0))
                assert j == 2*orbit.B + 2*orbit.P - 4*orbit.R + kappa
                yield (cohort, z, r1), guard, switch, block, failed, q2, j, kappa


def q3(q2, j, r2):
    return orbit.translate(q2, (j-r2)/72, orbit.SECOND_DRIFT)


def delta_at(row, delta):
    return orbit.integer(sum(Q(x, orbit.SCALE)*dx for x, dx in zip(row.terms, delta)))
