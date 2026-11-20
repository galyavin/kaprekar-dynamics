"""Fresh exact integration of the two-step anchors, row 181 geometry, and ledger sums.

Run with Python 3 from the repository root. No repository code is imported.
The historical three-step cell integrals are read as input; their rational
sum, counts, and rounding are checked here.
"""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import json
import re

DATA = Path("kaprekar-investigation/eleven_access")
ZERO = (0,) * 5
ONE = {ZERO: 1}


def variable(index):
    exponent = list(ZERO)
    exponent[index] = 1
    return {tuple(exponent): 1}


def add(left, right, scale=1):
    result = defaultdict(int, left)
    for power, value in right.items():
        result[power] += scale * value
    return {power: value for power, value in result.items() if value}


def multiply(left, right):
    result = defaultdict(int)
    for ex, a in left.items():
        for ey, b in right.items():
            power = tuple(x + y for x, y in zip(ex, ey))
            result[power] += a * b
    return {power: value for power, value in result.items() if value}


def folded_density():
    """Sum products of gap lengths of the chosen lifts of the five ordered lows."""
    total = {}
    for mask in product((False, True), repeat=5):
        high = [add(ONE, variable(i), -1) for i in range(5) if mask[i]]
        low = [variable(i) for i in reversed(range(5)) if not mask[i]]
        points = [ONE, *high, *low, {}]
        weight = ONE
        for a, b in zip(points, points[1:]):
            weight = multiply(weight, add(a, b, -1))
        total = add(total, weight)
    return total


def primitive_at(polynomial, index, other=None, constant=0, slope=0):
    """Evaluate the antiderivative in t[index] at constant + slope*t[other]."""
    result = defaultdict(Q)
    for power, coefficient in polynomial.items():
        n = power[index] + 1
        for j in range(n + 1):
            factor = Q(coefficient, n) * comb(n, j) * constant**(n-j) * slope**j
            if factor:
                exp = list(power)
                exp[index] = 0
                if j:
                    exp[other] += j
                result[tuple(exp)] += factor
    return {power: value for power, value in result.items() if value}


def anchor(U):
    # t1 >= 1-U and t4+t5 <= U are the two parts of the second-step anchor.
    lower = 1 - U
    polynomial = folded_density()
    for i in range(3):
        polynomial = add(
            primitive_at(polynomial, i, i+1, slope=1),
            primitive_at(polynomial, i, constant=lower), -1)
    total = 0
    for lo, hi, c, h in ((lower, U/2, 0, 1), (U/2, Q(1, 2), U, -1)):
        integrated = add(primitive_at(polynomial, 3, 4, c, h),
                         primitive_at(polynomial, 3, constant=lower), -1)
        final = add(primitive_at(integrated, 4, constant=hi),
                    primitive_at(integrated, 4, constant=lo), -1)
        total += final.get(ZERO, 0)
    return factorial(11) * total


def row181(t):
    return [(Q(7, 8), Q(3, 4), Q(5, 8), t, Q(1, 4)),
            (Q(7, 8), Q(5, 8), Q(1, 2), Q(3, 8), Q(5, 8)-t),
            (Q(7, 8), t+Q(1, 4), t, Q(1, 4), Q(1, 8))]


def row_geometry():
    lo, hi = Q(3, 8), Q(1, 2)
    mid = (lo + hi) / 2
    assert row181(Q(13, 32))[0] == tuple(map(Q, json.loads(
        (DATA/"entry_pi/catalogue.json").read_text())[181]["witness"]))
    # The sorting branch at mid stays valid on [lo,hi] by the affine
    # endpoint comparisons, proving the displayed affine images everywhere.
    for j, q in enumerate(row181(mid)):
        digits = (1, *q, *(1-x for x in q))
        order = sorted(range(11), key=digits.__getitem__, reverse=True)
        gaps = [a-b for a, b in zip((1, *q), (*q, 0))]
        gap_order = sorted(range(6), key=gaps.__getitem__)
        # The two moving gaps cross at mid. Use endpoint multiset equality
        # with affine gap forms matched by their slopes, not midpoint ordering.
        actual_gaps, expected_gaps = [], []
        for endpoint in (lo, hi):
            phase = row181(endpoint)[j]
            digits = (1, *phase, *(1-x for x in phase))
            ordered = [digits[i] for i in order]
            assert ordered == sorted(digits, reverse=True)
            image = tuple(ordered[i]-ordered[-1-i] for i in range(5))
            assert image == row181(endpoint)[(j+1)%3]
            actual_gaps.append(tuple(a-b for a,b in zip((1,*phase),(*phase,0))))
            expected_gaps.append((Q(1,8),Q(1,8),Q(1,8),Q(1,4),Q(5,8)-endpoint,endpoint-Q(1,4)))
        assert sorted(zip(*actual_gaps)) == sorted(zip(*expected_gaps))
        assert sum(gaps[i] for i in gap_order[:3]) == Q(3,8)
        assert gaps[gap_order[-1]] == Q(1,4)


def mass_totals():
    log = (DATA/"e1c_row181/verification.log").read_text()
    pat = r'^\([0-5], [0-5], [0-5]\) [0-5]: (\S+) \((\d+) paths, (\d+) domains\)'
    cells = [(Q(w),int(h),int(j)) for w,h,j in re.findall(pat,log,re.M)]
    assert (len(cells),sum(w>0 for w,h,j in cells),
            sum(h for w,h,j in cells),sum(j for w,h,j in cells)) == (25,23,31648,989)
    event = sum(w for w,h,j in cells)
    bound = 1-Q(5124467,80621568)-event
    assert event == Q(re.search(r'^Row 181 upper mass \(exact necessary event\): (\S+)',log,re.M)[1])
    assert bound == Q(re.findall(r'^m\(D\) lower bound \(exact\): (\S+)',log,re.M)[-1])
    assert 10**12*bound//1 == 929921494990
    assert -((-10**12*event)//1) == 6516518719
    print("Three-step event mass",event)
    print("Original entrance lower bound",bound)
    first_four = Q(json.loads((DATA/'en11b/entrance.json').read_text())[4]['by_step_each_class'])
    grid = (DATA/'e1_basins/mass.out').read_text().splitlines()[1:]
    additional = Q(factorial(11)*sum(int(line.split()[2]) for line in grid
                                           if int(line.split()[0]) >= 5),256**5*3072**6)
    assert additional == Q(16984163369908523775525,604462909807314587353088)
    assert 10**12*first_four//1 == 133940217689
    assert 10**12*(first_four+additional)//1 == 162038159009


if __name__ == "__main__":
    for U, expected in ((Q(5,6), Q(5124467,80621568)), (Q(7,8), Q(3882730875,17179869184)),
                        (Q(1),Q(1))):
        value = anchor(U)
        assert value == expected
        print("Two-step anchor",U,value)
    row_geometry()
    mass_totals()
