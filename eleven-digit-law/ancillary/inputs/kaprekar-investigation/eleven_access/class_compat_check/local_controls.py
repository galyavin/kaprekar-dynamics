"""Compute the localized box and cone distinctions from exact leaf geometry."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
import json
from math import gcd, lcm
from pathlib import Path

ACCESS = Path(__file__).resolve().parents[1]
ROOT = ACCESS.parents[1]


def combine(a, source, b, target):
    return tuple(a*x + b*y for x, y in zip(source, target))


def bezout(a, b):
    if b == 0:
        return abs(a), (1 if a >= 0 else -1), 0
    d, u, v = bezout(b, a % b)
    return d, v, u - (a // b) * v


def fibre(columns):
    """Reduce the base row by unimodular column operations, then take its fibre."""
    pivot, *transverse = columns
    kernel = []
    for column in transverse:
        a, b = pivot[0], column[0]
        if a == b == 0:
            kernel.append(tuple(column))
            continue
        d, u, v = bezout(a, b)
        kernel.append(combine(-b // d, pivot, a // d, column))
        pivot = combine(u, pivot, v, column)
    assert pivot[0] and all(v[0] == 0 for v in kernel)
    x, y = kernel
    covolume = abs(x[1]*y[2] - x[2]*y[1])
    assert covolume > 0
    return abs(pivot[0]), covolume


def dot(left, right):
    return sum(a*b for a, b in zip(left, right))


def box_lift():
    slopes = json.loads((ACCESS / "box_structure" / "structure.json").read_text())["slope_types"]
    lifts = []
    for section in slopes:
        matrix = [tuple(map(Q, row)) for row in section["coefficients"]]
        drift = tuple(map(Q, section["per_step"]))
        # Preserve the anchor a by setting dp=18. The centered checkpoint
        # must equal the old one advanced by 10 + 10*(60-4*dr)/36 steps.
        fixed_time = 10 + Q(10*60, 36)
        time_slope = Q(-40, 36)
        dr = (fixed_time*drift[1] - dot(matrix[1], (12, 18, 0)) + 8) / (
            matrix[1][2] - time_slope*drift[1])
        delta_j = 60 - 4*dr
        prefix = 10*delta_j/36
        centered = tuple(dot(row, (12, 18, dr))-2*(5-i) for i, row in enumerate(matrix))
        assert centered == tuple((prefix+10)*d for d in drift)
        assert (prefix+10) % section["period"] == 0
        assert delta_j % 72 == 36
        original_delta = (0, dr, 18, -6, dr-6)
        assert original_delta[1]-original_delta[4] == 6
        assert original_delta[2]+original_delta[3] == 12
        assert original_delta[0]+original_delta[2] == 18
        lifts.append((original_delta, prefix, prefix+10))
    assert all(lift == ((0, -12, 18, -6, -18), 30, 40) for lift in lifts)
    return {"source_delta": list(map(int, lifts[0][0])), "prefix_clock_change": 30,
            "centered_delay": 40, "elapsed_clock_change": -10}


def read_box():
    leaves = ROOT / "docs/drafts/research-direction/replay-task071/results/leaves.jsonl"
    bases = range(5, 72, 6)
    residual = {base: defaultdict(Q) for base in bases}
    primary = defaultdict(lambda: defaultdict(Q))
    marked = {(base, pin): defaultdict(Q) for base in bases for pin in (603, 639)}
    # Determinant of the three nontrivial inverse prefix columns after the
    # integral coordinate change documented in the structural certificate.
    index = abs(22 * (5*0 - (-2)*3))
    assert index == 132 and index*36*6 == 28512
    for line in leaves.open():
        row = json.loads(line)
        step, covolume = fibre(row["columns"])
        assert 72 % step == 0
        for base in bases:
            if (base - row["point"][0]) % step or base % 36 != row["base"]:
                continue
            clock = tuple(map(Q, row["profile"]))
            weight = Q(1, index*covolume)
            residual[base][clock] += weight
            label = (base, *row["label"], row["rho"])
            primary[label][clock] += weight
            if row["label"][:3] == [0, 1, 17] and row["rho"] == 2 and row["label"][3] in (603, 639):
                marked[base, row["label"][3]][clock] += weight
    for (base, cohort, z, first, pin, rho), weights in primary.items():
        assert sum(weights.values()) == Q(1, 28512)
        partner = ((base+12) % 72, cohort, z, first, 601+(pin-601+36) % 72, rho)
        assert primary.get(partner) == weights
    for base in bases:
        assert sum(residual[base].values()) == Q(1, 27)
        assert residual[base] == residual[5 if base % 12 == 5 else 11]
        for pin in (603, 639):
            assert sum(marked[base, pin].values()) == Q(1, 28512)
        partner = (base+12) % 72
        assert marked[base, 603] == marked[partner, 639]
        assert marked[base, 639] == marked[partner, 603]
    short = (Q(-11, 324), Q(-133, 162), Q(-35, 81))
    middle = (Q(29, 54), Q(-38, 27), Q(-20, 27))
    long = {(Q(139, 144), Q(-133, 72), Q(-35, 36)): Q(2, 3),
            (Q(863, 918), Q(-836, 459), Q(-440, 459)): Q(1, 3)}
    for base in bases:
        other = {middle if base % 12 == 5 else short: Q(1)}
        for pin in (603, 639):
            got = {clock: w*28512 for clock, w in marked[base, pin].items()}
            is_long = (base % 24 in (11, 17)) == (pin == 603)
            assert got == (long if is_long else other)
    coeff = tuple(sum((residual[5][key]-residual[11][key])*key[i]
                      for key in set(residual[5]) | set(residual[11])) for i in range(3))
    archived_full = json.loads((ACCESS / "box_structure" / "limit_weights.json").read_text())["graphs"]
    common = {}
    for base in (5, 11):
        common[base] = {tuple(map(Q, row["s"])): Q(row[f"arm_{base}_units"])/28512
                        - residual[base][tuple(map(Q, row["s"]))] for row in archived_full}
        assert min(common[base].values()) >= 0 and sum(common[base].values()) == Q(26, 27)
    assert common[5] == common[11]
    bounds = [[Q(entry, 10**5) for entry in interval] for interval in (
        (94040, 94130), (73770, 73850), (57940, 58015), (20, 80), (200, 300))]
    difference_values = []
    pi_margins = []
    plateau_margins = []
    center = (940800, 738087, 579763, 484, 2491)
    stretch = (446, 365, 342, 280, 449)
    slopes = (590, -307, -447, 127, -296)
    for vertex in product(*bounds):
        a, b, c, h, g = vertex
        pi_margins.append(pi_slack((a, b, c, 1-c-g, b-(1+h)/2)))
        count = (a+c-Q(3, 2)-Q(11, 2)*h-5*g)/22
        p, r = c-h/2-2*count, b-Q(3, 2)*h-g-4*count
        difference_values.append(coeff[0] + coeff[1]*p + coeff[2]*r)
        theta = Q(121049, 10**6) + sum(
            Q(v, s)*(coord-Q(u, 10**6)) for coord, u, s, v in zip(vertex, center, stretch, slopes))
        for clock in common[5]:
            x, y, z = clock
            time = a+x+y*p+z*r
            assert 0 < time < 1
            margin = (time-theta if clock in long else theta-time) - Q(1, 1000)
            plateau_margins.append(margin)
    assert min(difference_values) > 0
    assert min(plateau_margins) > 0

    def table(weights):
        return [[list(map(str, clock)), str(weight)] for clock, weight in sorted(weights.items())]

    return dict(residual={base: table(residual[base]) for base in (5, 11)},
                tested_primary_rows=len(primary), lift=box_lift(),
                marked={f"{base}:{pin}": table(marked[base, pin])
                        for base in (5, 11, 17, 23) for pin in (603, 639)},
                min_pi_slack=str(min(pi_margins)),
                min_difference=str(min(difference_values)),
                max_difference=str(max(difference_values)),
                plateau_slack=str(min(plateau_margins)),
                plateau_difference=str(sum(residual[5][key]-residual[11][key] for key in long)))


def pi_slack(q):
    """Use the five common left-first paired phases and check each required guard."""
    margin = Q(1)
    start = q
    for _ in range(5):
        right = tuple(1-x for x in reversed(q))
        digits = (1, *(value for pair in zip(q, right) for value in pair))
        margin = min(margin, q[-1], 1-q[0],
                     *(a-b for a, b in zip(digits, digits[1:])))
        q = tuple(digits[i]-digits[-1-i] for i in range(5))
    assert margin > 0 and q == start
    return margin


def read_cone():
    ledger = json.loads((ACCESS / "stage_graph" / "cell-lotteries.json").read_text())
    xstar = tuple(Q(i, 50000) for i in (2938, 3225, 2886, 1913, 387))
    ray = [(Q(1),)] + [tuple((x + offset) / xstar[0] for offset in (-Q(1, 10**5), Q(1, 10**5)))
                       for x in xstar[1:]]
    mean = {}
    for root in (3, 5, 9):
        rows = ledger["tables"][ledger["mod12_to_table"][str(root)]]["clock"]
        assert sum(Q(row["weight"]) for row in rows) == 1
        mean[root] = tuple(sum(Q(row["weight"]) * Q(row["numerator"][i], row["denominator"])
                               for row in rows) for i in range(6))
    assert mean[3] == mean[9]
    differences = []
    pi_margins = []
    for scale in (Q(1, 128), Q(1, 8)):
        for direction in product(*ray):
            x = tuple(scale*v for v in direction)
            gaps2 = [a-b for a, b in zip((*x, 0), (0, *x))]
            assert max(gaps2) < Q(1, 7)
            pi_margins.append(pi_slack(tuple((5-i+x[i])/6 for i in range(5))))
            for root in (3, 5):
                for row in ledger["tables"][root]["clock"]:
                    time = sum(Q(a, row["denominator"])*b for a, b in zip(row["numerator"], (*x, 1)))
                    assert 0 < time < 1
            differences.append(sum((a-b)*c for a, b, c in zip(mean[9], mean[5], (*x, 1))))
    assert min(differences) > 0
    speed = {root: vec[1]/2 for root, vec in mean.items()}
    return dict(mean_speeds={root: str(value) for root, value in speed.items()},
                min_pi_slack=str(min(pi_margins)),
                min_difference=str(min(differences)), max_difference=str(max(differences)))


def read_normal_controls():
    rows = json.loads((ACCESS / "class_compat" / "premises.json").read_text())["normal"]
    strides = []
    for record in rows:
        _, period, _, *data = record["sheet"]
        shift, column = data[:5], data[-5:]
        tangent = tuple(x // gcd(*column) for x in column)
        # A primitive integer tangent is the intersection of its real line
        # and the integer lattice. Test all offsets in Z^5 / period Z^5.
        direct = next(6*m for m in range(1, period+1)
                      if any(all((m*s-k*u) % period == 0 for s, u in zip(shift, tangent))
                             for k in range(period)))
        divisor, left_inverse = 0, []
        for coordinate in tangent:
            divisor, old, new = bezout(divisor, coordinate)
            left_inverse = [old*x for x in left_inverse] + [new]
        assert divisor == 1
        projected = [Q(s - u*sum(a*b for a, b in zip(shift, left_inverse)), period)
                     for s, u in zip(shift, tangent)]
        predicted = 6*lcm(*(x.denominator for x in projected))
        assert predicted == direct == record["direct_period"]
        strides.append(direct)
    return strides


def main():
    print(json.dumps({"box": read_box(), "cone": read_cone(),
                      "normal_stride_controls": read_normal_controls()}, indent=2))


if __name__ == "__main__":
    main()
