"""Integrate the original density and independently project the marked box leaves."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations, product
import json
from math import factorial, gcd, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def difference(left, right):
    return tuple(a-b for a,b in zip(left,right))


def mass(bounds, coordinates, jacobian, factor=None):
    unit = (1, 0, 0, 0, 0, 0)
    zero = (0,) * 6
    augmented = (unit, *coordinates, zero)
    gaps = [difference(a,b) for a,b in zip(augmented,augmented[1:])]
    if factor is not None:
        gaps.append(factor)
    moments = [[(u**(i+1)-l**(i+1))/(i+1) for i in range(len(gaps)+1)] for l,u in bounds]
    integral = Q(0)
    for monomial in product(*(tuple((j,a) for j,a in enumerate(g) if a) for g in gaps)):
        exponents = [0]*6
        for j,_ in monomial:
            exponents[j] += 1
        integral += prod(a for _,a in monomial)*prod(
            values[exponent] for values,exponent in zip(moments,exponents[1:]))
    return factorial(11)*jacobian*integral


def determinant(columns):
    return sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)) *
               prod(columns[i][p[i]] for i in range(3)) for p in permutations(range(3)))


def box_profiles():
    archive = ROOT/"docs/drafts/research-direction/replay-task071/results/leaves.jsonl"
    table = defaultdict(lambda: defaultdict(Q))
    with archive.open() as leaves:
        for line in leaves:
            leaf = json.loads(line)
            cohort, overshoot, first, pin = leaf["label"]
            if (cohort,overshoot,first,leaf["rho"]) != (0,1,17,2) or pin not in (603,639):
                continue
            step = gcd(*(col[0] for col in leaf["columns"]))
            fibre = Q(abs(determinant(leaf["columns"])),step)
            for base in range(5,72,6):
                if (base-leaf["point"][0])%step==0:
                    table[base,pin][tuple(leaf["profile"])] += 216/fibre
    short = ("-11/324", "-133/162", "-35/81")
    mid = ("29/54", "-38/27", "-20/27")
    long_profiles = {("139/144", "-133/72", "-35/36"): Q(2,3),
                     ("863/918", "-836/459", "-440/459"): Q(1,3)}
    out = {}
    for base in range(5,72,6):
        low = table[base,603]
        high = table[base,639]
        expected_other = {mid if base%12==5 else short: Q(1)}
        expected_low = long_profiles if base%24 in (11,17) else expected_other
        expected_high = expected_other if base%24 in (11,17) else long_profiles
        assert low == expected_low and high == expected_high
        out[base] = {str(pin): [dict(coefficients=profile, weight=str(weight))
                               for profile,weight in sorted(weights.items())]
                     for pin,weights in ((603,low),(639,high))}
    return out


def clock_marginal(bounds, coordinates, box):
    source = HERE.parent/"box_structure/limit_weights.json"
    graphs = json.loads(source.read_text())["graphs"]
    assert all(Q(row["signed_units"])==Q(row["arm_5_units"])-Q(row["arm_11_units"])
               for row in graphs)
    signed = [Q(row["signed_units"])/28512 for row in graphs]
    assert sum(signed)==0
    s0, sp, sr = [sum(w*Q(row["s"][i]) for w,row in zip(signed,graphs)) for i in range(3)]

    def chart(point):
        a,b,c,h,g = point
        n = (a+c-Q(3,2)-Q(11,2)*h-5*g)/22
        return c-h/2-2*n, b-Q(3,2)*h-g-4*n

    def gap(point):
        p,r = chart(point)
        return s0+sp*p+sr*r

    values = [gap(v) for v in product(*bounds)]
    assert min(values)>0
    for v in product(*bounds):
        p,r = chart(v)
        for row in graphs:
            c0,cp,cr = map(Q,row["s"])
            assert 0<v[0]+c0+cp*p+cr*r<1
    constant = gap((Q(0),)*5)
    affine = (constant, *(gap(tuple(Q(i==j) for j in range(5)))-constant for i in range(5)))
    contribution = mass(bounds,coordinates,Q(1,2),affine)
    assert min(values)<=contribution/box<=max(values)
    return dict(signed_table_clock_coefficients=list(map(str,(s0,sp,sr))),
                box_pointwise_bounds=list(map(str,(min(values),max(values)))),
                conditional_mean=str(contribution/box), original_contribution=str(contribution))


def main():
    center = [Q(i,50000) for i in (2938,3225,2886,1913,387)]
    radius = Q(1,100000)
    cube = mass([(x-radius,x+radius) for x in center],
                [tuple([Q(5-i,6)]+[Q(i==j,6) for j in range(5)]) for i in range(5)], Q(1,6**5))
    box_bounds = [(Q(a,100000), Q(b,100000))
                  for a,b in ((94040,94130),(73770,73850),(57940,58015),(20,80),(200,300))]
    box_coordinates = [(0,1,0,0,0,0), (0,0,1,0,0,0), (0,0,0,1,0,0),
                       (1,0,0,-1,0,-1), (Q(-1,2),0,1,0,Q(-1,2),0)]
    box = mass(box_bounds,box_coordinates,Q(1,2))
    print(json.dumps(dict(cube_original_mass=str(cube), box_original_mass=str(box),
                          each_marked_label_original_mass=str(box/28512),
                          both_original_mass=str(cube+box),
                          clock_marginal_class5_minus11=clock_marginal(box_bounds,box_coordinates,box),
                          marked_profiles=box_profiles()),indent=2))


if __name__ == "__main__":
    main()
