"""Independently verify a full-rank projected-return route and its lift constants."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

R = 12
FORMS = []
for f in range(2 * R):
    row = [0] * R
    i = f % R
    row[i] = 1
    if f >= R:
        row[(i + 7) % R] = 1
        row[(i - 1) % R] = -1
    FORMS.append(row)


def read(f, x):
    i = f % R
    if f < R:
        return x[i]
    return x[i] + x[(i + 7) % R] - x[(i - 1) % R]


def rank(rows):
    a = [list(map(Fraction, row)) for row in rows]
    j = 0
    for i in range(R):
        p = next((p for p in range(j, len(a)) if a[p][i]), None)
        if p is None:
            continue
        a[j], a[p] = a[p], a[j]
        for t in range(j + 1, len(a)):
            if a[t][i]:
                factor = a[t][i] / a[j][i]
                a[t] = [v - factor * u for u, v in zip(a[j], a[t])]
        j += 1
    return j


def add(x, y, multiplier=1):
    return [u + multiplier * v for u, v in zip(x, y)]


def orbit_order(values, phase=0):
    return tuple(values[(phase - 1 - i) % R] for i in range(R))


def step(x, phase, enabled):
    # Derive by rotating in fold order and subtracting the updated anchor.
    y = x[-1]
    z = x[4] + y - x[0]
    if phase not in enabled:
        y = 1
    if R + phase not in enabled:
        z = -1
    e = int(y > 0) - int(y < 0)
    h1, h2 = max(0, -1-y), max(0, -1-z)
    corrections = [-e+h1] * R
    corrections[0] = 2*h1-3*e
    corrections[1] = 2*h1-h2
    corrections[2] += h2
    corrections[5] = -e+h2
    return tuple(a + b for a, b in zip((x[-1], *x[:-1]), corrections)), h2


def label_order(x, phase=0):
    return [x[(phase - 1 - i) % R] for i in range(R)]


def segment(state, length, tracked):
    x = orbit_order(state)
    minimum = [None] * (2 * R)
    largest = max(map(abs, x))
    export = 0
    for t in range(length):
        phase = t % R
        for f, v in ((phase, x[-1]), (R+phase, x[4]+x[-1]-x[0])):
            old = minimum[f]
            minimum[f] = v if old is None else min(v, old)
        x, h = step(x, phase, tracked)
        export += h
        largest = max(largest, max(map(abs, x)))
    assert length % R == 0
    return label_order(x), minimum, largest, export


def certify(route_path):
    lines = route_path.read_text().splitlines()
    header = lines[0].split()
    assert header[0] == "found"
    claimed_entry, claimed_period = int(header[-3]), int(header[-1])
    path = [list(map(int,line.split())) for line in lines[1:]]
    assert len(path) == R
    state = [0] * R
    tracked = set()
    chosen = []
    stage_data = []
    max_nominal = 0
    determinant_forms = 1
    for j, (f, prefix, period, *shift) in enumerate(path):
        assert len(shift)==R and f not in tracked and prefix>=0 and period>0
        assert all(read(t,shift)==0 for t in tracked)
        start, pre_min, pre_max, pre_export = segment(state, prefix, tracked)
        end, block_min, block_max, block_export = segment(start, period, tracked)
        ray = add(start, end, -1)
        assert all(read(t,ray)==0 for t in tracked)
        chosen.append(f)
        rows = [FORMS[t] for t in chosen]
        assert rank(rows) == j+1
        bigger = {t for t, q in enumerate(FORMS) if rank(rows+[q])==j+1}
        new = bigger-tracked
        for t in new:
            assert read(t, ray)>0
            assert block_min[t]+read(t, shift)>=(1 if t<R else -1)
        determinant_forms *= read(f, ray)
        stage_data.append({
            "added": f, "fresh": sorted(new), "prefix": prefix, "period": period,
            "ray": ray, "shift": shift, "prefix_export": pre_export, "word_export": block_export,
            "prefix_minimum": pre_min, "word_minimum": block_min,
            "largest_nominal": max(pre_max,block_max),
        })
        max_nominal = max(max_nominal, pre_max, block_max)
        tracked = bigger
        state = add(end, shift)
    assert len(tracked)==2*R
    assert determinant_forms > 0

    terminal_source = orbit_order(state)
    x = terminal_source
    seen = {}
    clock, terminal_export = 0, 0
    while x not in seen:
        seen[x] = clock
        next_x, h2 = step(x, clock%R, tracked)
        terminal_export += h2
        max_nominal = max(max_nominal, max(map(abs,x)))
        clock += 1
        x = next_x
        assert clock <= claimed_entry + claimed_period
    assert seen[x] == claimed_entry and clock-seen[x] == claimed_period
    cycle = []
    for t in range(claimed_period):
        cycle.append(x)
        x,h2 = step(x,t%R,tracked)
        assert h2 == 0
    assert x==cycle[0]

    # Coefficient floors by ownership. For each newly fixed read f, verify
    # that earlier stages and the owner's prefix are feasible for arbitrary
    # counts in the ranges [floor_j*L,(floor_j+1)*L).
    floors = []
    for i, stage in enumerate(stage_data):
        minimum_floor = 0
        for f in stage["fresh"]:
            threshold = 1 if f<R else -1
            slopes = [read(f,s["ray"]) for s in stage_data[:i]]
            earlier_shifts = [read(f,s["shift"]) for s in stage_data[:i+1]]
            required = 0
            for j in range(i+1):
                minima = [stage_data[j]["prefix_minimum"][f]]
                if j<i:
                    minima += [stage_data[j]["word_minimum"][f]]
                read_floor = min((a for a in minima if a is not None), default=None)
                if read_floor is None:
                    continue
                constant = read_floor+sum(earlier_shifts[j:])
                negative = sum(min(0,slope)*(floors[t]+1) for t,slope in enumerate(slopes) if t>=j)
                required = max(required, max(0,threshold-constant)-negative)
            slope = read(f,stage["ray"])
            minimum_floor = max(minimum_floor, -(-required//slope))
        floors.append(minimum_floor)
    norm = max_nominal
    loss_budget = terminal_export
    time_budget = clock
    for floor, stage in zip(floors,stage_data):
        norm += max(map(abs,stage["shift"])) + (floor+1)*max(map(abs,stage["ray"]))
        loss_budget += stage["prefix_export"] + (floor+2)*stage["word_export"]
        time_budget += stage["prefix"] + (floor+2)*stage["period"]
    for i,stage in enumerate(stage_data):
        stage["coefficient_floor"] = floors[i]

    C = max(2*norm,loss_budget)
    n=35
    m=18
    k=17
    J={0,1,2,4,8,16,3,6,12,11,13,9,17}
    pi=[Fraction(3,55) if j in J else Fraction(16,275) for j in range(m)]
    assert sum(pi)==1
    margins=[sum(pi[m-i:])-sum(pi[:i]) for i in range(1,m//2+1)]
    delta_pair=min(pi)-Fraction(1,m+1)
    radius=max(abs(sum(pi[i:])-Fraction(m-i,m)) for i in range(1,m))
    gamma=min(margins[2:])
    assert delta_pair>0 and gamma>0 and radius>0
    E=m+2*C+1
    D=math.ceil(4*E/min(delta_pair,gamma/m,radius/m))
    B0=max(2*D,math.ceil(Fraction(4,3*(m+1)*delta_pair)),2*(m+1))
    return {
        "rank": R,
        "determinant_in_read_coordinates": determinant_forms,
        "stages": stage_data,
        "terminal_source": terminal_source,
        "terminal_entry_from_release": claimed_entry,
        "terminal_period": claimed_period,
        "terminal_cycle": cycle,
        "terminal_export_to_entry": terminal_export,
        "time_through_a_period_coefficient": time_budget,
        "spread_coefficient": 2*norm,
        "export_coefficient": loss_budget,
        "common_coefficient": C,
        "lift_profile": [str(p) for p in pi],
        "lift_rows": [str(a) for a in margins],
        "lift_delta_pair": str(delta_pair),
        "lift_radius": str(radius),
        "lift_cut": str(radius/2),
        "lift_D": D,
        "certified_base_threshold": B0,
        "original_probability_lower_bound_numerator": math.factorial(n),
        "original_probability_lower_bound_denominator": (2*D)**k*(2*(m+1))**m,
    }


if __name__=="__main__":
    print(json.dumps(certify(Path(sys.argv[1])),indent=2))
