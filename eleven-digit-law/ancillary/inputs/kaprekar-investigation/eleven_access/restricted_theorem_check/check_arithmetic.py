"""Exact integrals, box transport algebra, and clock constants."""
import json
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb, factorial, lcm, prod
from pathlib import Path

from check_controls import BOX, BOX_VERTICES, P_FORM, R_FORM
from exact_model import CENTRE, combine, leading

ACCESS = Path(__file__).resolve().parents[1]
UNIT = (0, 0, 0, 0, 0, 1)
COORDS = tuple(tuple(int(i == j) for j in range(6)) for i in range(5))


def integrate(rows, box):
    polynomial = {(0,)*5: Q(1)}
    for row in rows:
        new = {}
        for alpha, coef in polynomial.items():
            for j, a in enumerate(row):
                if a:
                    key = tuple(k + int(i == j) for i, k in enumerate(alpha))
                    new[key] = new.get(key, Q(0)) + coef*a
        polynomial = new
    return sum(c * prod((b**(n+1)-a**(n+1))/(n+1) for (a, b), n in zip(box, alpha))
               for alpha, c in polynomial.items())


def gaps(rows):
    return tuple(combine(1, a, -1, b) for a, b in zip((UNIT, *rows), (*rows, (0,)*6)))


def data_checks():
    exact = json.loads((ACCESS/'r2_routes/control-weights.json').read_text())
    a, b, c, d, e = COORDS
    chart = (a, b, c, combine(1, UNIT, -1, combine(1, c, 1, e)),
             combine(1, b, Q(-1, 2), combine(1, UNIT, 1, d)))
    density = gaps(chart)
    assert all(leading(g, x) > 0 for g in density for x in BOX_VERTICES)
    box_mass = Q(factorial(11), 2) * integrate(density, BOX)
    cube = tuple((Q(a, 50000)-Q(1, 100000), Q(a, 50000)+Q(1, 100000))
                 for a in (2938, 3225, 2886, 1913, 387))
    cube_rows = tuple(combine(Q(1, 6), row, centre, UNIT) for row, centre in zip(COORDS, CENTRE))
    cube_mass = Q(factorial(11), 6**5) * integrate(gaps(cube_rows), cube)
    # In gapped source coordinates only the N,n,m columns of the inverse prefix remain.
    prefix_index = abs(22 * (5*0 - (-2)*3))
    marked_index = prefix_index * 36 * 6
    assert prefix_index == 132 and marked_index == 28512
    for key, value in [('box_original_mass', box_mass), ('cube_original_mass', cube_mass),
                       ('each_marked_label_original_mass', box_mass/marked_index)]:
        print(key, value)
        assert value == Q(exact[key])

    # Reconstruct the centered equal-anchor displacement under the box twelve-shift.
    matrices = [
        ((Q(17,18), Q(-14,9), Q(10,9)), (Q(-4,9),Q(-4,9),Q(17,9)),
         (Q(-2,3),Q(1,3),Q(4,3)), (Q(4,9),Q(-14,9),Q(10,9)),
         (Q(-17,18),Q(-4,9),Q(17,9))),
        ((Q(17,18),Q(-14,9),Q(10,9)), (Q(37,18),Q(-31,9),Q(8,9)),
         (Q(3),Q(-3),Q(-1)), (Q(25,9),Q(-20,9),Q(-14,9)),
         (Q(25,18),Q(-10,9),Q(-7,9))),
    ]
    velocities = [(-1,Q(-11,10),Q(-3,5),-1,Q(-11,10)), (-1,Q(-7,5),Q(-3,10),Q(1,5),Q(1,10))]
    hints = []
    for matrix, speed in zip(matrices, velocities):
        def residual(r):
            return tuple(sum(coef*var for coef, var in zip(row, (12,18,r))) - 12*c
                         - (10 + Q(10,36)*(60-4*r))*d for row,c,d in zip(matrix,CENTRE,speed))
        f0, f1 = residual(0)[1], residual(1)[1]
        r = -f0/(f1-f0)
        assert r == -12 and residual(r) == (0,)*5
        hints.append(r)
    delta_j = 60 - 4*hints[0]
    assert delta_j == 108 and delta_j % 72 == 36
    print('whole-box twelve-shift', (0, -12, 18, -6, -18), 'delta J', delta_j,
          'physical clock change', Q(10,36)*delta_j-(10+Q(10,36)*delta_j))

    key = exact['clock_marginal_class5_minus11']
    leaf_forms = [tuple(map(Q, row)) for row in (
        ('-11/324','-133/162','-35/81'), ('-7/408','-57/68','-15/34'),
        ('1/576','-247/288','-65/144'), ('31/96','-19/16','-5/8'),
        ('29/54','-38/27','-20/27'), ('139/144','-133/72','-35/36'), ('863/918','-836/459','-440/459'))]
    counts = tuple(map(Q, (0, -21, -35, 0, -12, '191/3', '13/3')))
    assert sum(counts) == 0
    signed = tuple(sum(weight*a for weight,a in zip(counts, column))/28512 for column in zip(*leaf_forms))
    assert signed == tuple(map(Q,key['signed_table_clock_coefficients']))
    row = combine(1, combine(signed[1], P_FORM, signed[2], R_FORM), signed[0], UNIT)
    extremals = [leading(row, v) for v in BOX_VERTICES]
    assert (min(extremals), max(extremals)) == tuple(map(Q,key['box_pointwise_bounds']))
    effect = Q(factorial(11),2)*integrate((*density,row), BOX)
    assert effect == Q(key['original_contribution'])
    assert effect/box_mass == Q(key['conditional_mean'])
    print('signed box clock contribution', effect)


def borrow_checks():
    counts = dict(strict_paired=0, tied=0, unpaired=0, shorter=0)
    for base in range(2, 15):
        for ascending in combinations_with_replacement(range(base), 5):
            q = tuple(reversed(ascending))
            if not q[0]:
                continue
            height = sum(a > 0 for a in q)
            used = q[:height]
            digits = (*used[:-1], used[-1]-1, *([base-1]*(11-2*height)),
                      *(base-a-1 for a in reversed(used[1:])), base-used[0])
            value = sum(q[i]*(base**(10-i)-base**i) for i in range(5))
            schoolbook = []
            for _ in range(11):
                value, digit = divmod(value, base)
                schoolbook.append(digit)
            assert not value and tuple(reversed(schoolbook)) == digits
            if height < 5:
                counts['shorter'] += 1
                continue
            pairs = tuple((digits[i], digits[6+i]) for i in range(5))
            if any(min(a) < max(b) for a, b in zip(pairs, pairs[1:])):
                counts['unpaired'] += 1
            elif len(set(digits)) < 11:
                counts['tied'] += 1
            else:
                counts['strict_paired'] += 1
    assert all(counts.values())
    print('exact small-base borrowed-column counts', counts)


def offset_stride():
    # A toy nonhomogeneous sheet exercises compatible-base refinements.
    a = (Q(1,2), Q(0), Q(1,3), Q(-5,4), Q(2,5))
    stride = 6*lcm(*(x.denominator for x in a))
    assert stride == 360
    actual = [j for j in range(6, 2*stride+1, 6)
              if all((j*z/6).denominator == 1 for z in a)]
    assert actual == [stride, 2*stride]
    print('normal stride control', stride)


def clock_checks():
    assert comb(509,5) == 279158015851
    assert 1-Q(3,16)*Q(5,48) == Q(251,256)
    assert Q(32,6) == Q(16,3)
    print('terminal constants', comb(509,5), Q(251,256), Q(16,3))


if __name__ == '__main__':
    data_checks()
    borrow_checks()
    offset_stride()
    clock_checks()
