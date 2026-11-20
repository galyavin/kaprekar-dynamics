"""Compose the certified first two stages and integrate a signed terminal cone."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent
DIM = 5
ZERO = (0,) * DIM


def add(a, b):
    out = a.copy()
    for key, value in b.items():
        out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}


def scale(a, k):
    return {key: value * k for key, value in a.items() if value * k}


def subtract(a, b):
    return add(a, scale(b, -1))


def multiply(a, b):
    out = {}
    for left, v in a.items():
        for right, w in b.items():
            key = tuple(x+y for x, y in zip(left, right))
            out[key] = out.get(key, Q(0)) + v*w
    return {key: value for key, value in out.items() if value}


def constant(c):
    return {ZERO: Q(c)} if c else {}


def coordinate(j):
    return {tuple(int(i==j) for i in range(DIM)): Q(1)}


def evaluate(poly, point):
    return sum(v*product_powers(point, key) for key, v in poly.items())


def product_powers(point, key):
    result = Q(1)
    for x, power in zip(point, key):
        if power:
            result *= x**power
    return result


def integrate(poly, bounds):
    moments = [[(hi**(j+1)-lo**(j+1))/(j+1) for j in range(12)] for lo, hi in bounds]
    total = Q(0)
    for key, value in poly.items():
        for m, i in zip(moments, key):
            value *= m[i]
        total += value
    return total


def gaps(q):
    return [subtract(a,b) for a,b in zip([constant(1), *q], [*q, {}])]


def density(q):
    out = constant(factorial(11))
    for gap in gaps(q):
        out = multiply(out, gap)
    return out


def in_pi(q):
    for _ in range(5):
        if not (0 < q[4] and q[0] < 1 and all(q[i]+q[4-i] > 1 for i in range(5))
                and all(q[i+1]+q[4-i] < 1 for i in range(4))):
            return False
        q = (q[0],q[0]-q[4],q[1]-q[4],q[1]-q[3],q[2]-q[3])
    return True


def check_prefix():
    source = json.loads((INPUT/'stage_graph/sections.json').read_text())
    result = []
    for f in source['faces']:
        tables = f['stage2_by_B_mod_12']
        signatures = set()
        for b, lanes in tables.items():
            first_rem = set(map(int, lanes))
            expected = ({6,18,30} if int(b)%2==0 else {0,12,24}) if f['first_pair']==3 else set(range(0,36,6))
            assert first_rem == expected
            for lane in lanes.values():
                period = lane['period']
                assert lane['short_steps']==0 and lane['failures']==0
                drift = tuple(Q(d,period) for d in lane['drift'])
                guards = []
                for winner in lane['stage2_winners']:
                    if not Q(winner['conditional_volume']):
                        continue
                    row = tuple(period*Q(x) for x in winner['gap_time_row'])
                    guards.append((row,Q(winner['conditional_volume']),winner['leads_to_interlevel_contact']))
                    assert sum(len(form['exit_remainders'])*Q(form['each_lane']) for form in winner['forms'])==1
                signatures.add((drift,tuple(sorted(guards))))
        assert len(signatures)==1
        drift, guards = signatures.pop()
        result.append(dict(first_phase=f['first_phase'], first_pair=f['first_pair'],
                           drift_per_update=drift, clock_cells=guards))
    assert len(result)==10
    return result


def exact_cone():
    center = tuple(Q(i,50000) for i in (2938,3225,2886,1913,387))
    radius = Q(1,100000)
    bounds = [(Q(1,128),Q(1,8)), *[((x-radius)/center[0],(x+radius)/center[0]) for x in center[1:]]]
    t = coordinate(0)
    x = (t, *(multiply(t,coordinate(j)) for j in range(1,5)))
    q = [scale(add(constant(5-j), x[j]),Q(1,6)) for j in range(5)]
    verts = list(product(*bounds))
    for v in verts:
        qv = tuple(evaluate(expr,v) for expr in q)
        assert in_pi(qv)
        g = [evaluate(expr,v) for expr in gaps(q)]
        assert min(g)>Q(1,7), (v,g) # Equivalent to R(X/B)<1/7.
    d = add(add(scale(x[0],Q(-4,5)),scale(x[1],2)),subtract(x[3],scale(x[2],2)))
    assert min(evaluate(d,v) for v in verts)>0
    forms = json.loads((INPUT/'stage_graph/cell-lotteries.json').read_text())['tables']
    speeds = []
    stage_times = []
    for table in forms:
        weights = {}
        assert not table['gap']
        for leaf in table['clock']:
            w = Q(leaf['numerator'][3],leaf['denominator'])
            assert tuple(Q(a,leaf['denominator']) for a in leaf['numerator']) == (
                Q(1,6)-Q(4,5)*w, 2*w, -2*w, w, 0, 0)
            assert w not in weights
            weights[w] = Q(leaf['weight'])
            assert all(0 < evaluate(add(scale(t,Q(1,6)),scale(d,w)),v) < 1 for v in verts)
        assert sum(weights.values())==1
        speeds.append(sum(w*p for w,p in weights.items()))
        coeff = defaultdict(lambda: [Q(0)]*6)
        for edge in table['edge']:
            assert edge['failures']==0
            weight = Q(edge['weight'])*Q(edge['period'],edge['decrement']*edge['distance_denominator'])
            for i,v in enumerate(edge['distance_numerator']):
                coeff[edge['stage']][i] += weight*v
        expected = [Q(1,6)-Q(4,5)*speeds[-1], 2*speeds[-1], -2*speeds[-1], speeds[-1], 0, 0]
        assert [sum(c[i] for c in coeff.values()) for i in range(6)]==expected
        stage_times.append(coeff)
    for k in range(1,6):
        assert all(coeff[k]==stage_times[0][k] for coeff in stage_times)
    signed_times = {k: [a-b for a,b in zip(stage_times[3][k],stage_times[5][k])] for k in (6,7)}
    time_contrasts = {}
    for k,c in signed_times.items():
        out = constant(c[5])
        for expr,value in zip(x,c):
            out = add(out,scale(expr,value))
        time_contrasts[k] = out
    assert all(evaluate(time_contrasts[6],v)<0<evaluate(time_contrasts[7],v) for v in verts)

    g = gaps(q)
    low = [g[0],g[1],add(g[0],g[2]),add(g[1],g[3]),add(add(g[0],g[2]),g[4])]
    for v in verts:
        values = [evaluate(expr,v) for expr in low]
        assert 0 < values[0] < values[1] < values[2] < values[3] < values[4] < Q(1,2)
    inverse = []
    midpoint = tuple((lo+hi)/2 for lo, hi in bounds)
    inside = []
    for mask in product((0,1), repeat=5):
        pre = [subtract(constant(1),p) for p,b in zip(low,mask) if b]
        pre += [p for p,b in reversed(list(zip(low,mask))) if not b]
        if in_pi(tuple(evaluate(p,midpoint) for p in pre)):
            inside.append(mask)
        else:
            inverse.append(pre)
    assert len(inverse)==31 and len(inside)==1
    jac = {(4,0,0,0,0):Q(1,6**5)}
    rows = {}
    for name, domains in [('direct',[q]),('first_entry_one',inverse)]:
        numerator = {}
        for p in domains:
            numerator = add(numerator, density(p))
        weight = multiply(jac,numerator)
        mass, d_total = integrate(weight,bounds), integrate(multiply(weight,d),bounds)
        contrasts = {f'{i}-{j}':(speeds[i]-speeds[j])*d_total for i,j in combinations(range(6),2)}
        rows[name] = dict(mass=mass, distance_integral=d_total, signed_contributions=contrasts,
                          signed_3_minus_5_by_stage={k: integrate(multiply(weight,v),bounds)
                                                     for k,v in time_contrasts.items()})
    mass = sum(r['mass'] for r in rows.values())
    d_total = sum(r['distance_integral'] for r in rows.values())
    return dict(t_and_slopes=bounds, one_step_inside_mask=inside[0], mean_speeds=speeds,
                cohorts=rows, mass=mass, distance_integral=d_total, signed_3_minus_5_stage_rows=signed_times,
                conditional_3_minus_5=(speeds[3]-speeds[5])*d_total/mass,
                original_3_minus_5=(speeds[3]-speeds[5])*d_total)


def boxed_threshold():
    center = tuple(Q(i,50000) for i in (2938,3225,2886,1913,387))
    epsilon = Q(1,200000)
    low, high = Q(10383,1000000), Q(10386,1000000)
    law = json.loads((INPUT/'stage_graph/cell-lotteries.json').read_text())['tables']
    out = {}
    for b in (3,5):
        value = Q(0)
        endpoints = []
        for row in law[b]['clock']:
            c = [Q(a,row['denominator']) for a in row['numerator']]
            middle = sum(x*y for x,y in zip(c,center))+c[5]
            radius = sum(abs(x)*epsilon for x in c[:5])
            if middle-radius > high:
                value += Q(row['weight'])
            else:
                assert middle+radius < low, (b,row,middle+radius,low)
            endpoints.append([middle-radius,middle+radius])
        out[b] = dict(expectation=value, support_intervals=endpoints)
    assert out[3]['expectation']==Q(1093,7776) and out[5]['expectation']==0
    q = [scale(add(constant(5-j),coordinate(j)),Q(1,6)) for j in range(5)]
    mass = integrate(scale(density(q),Q(1,6**5)),[(c-epsilon,c+epsilon) for c in center])
    return dict(ramp=(low,high), pointwise_expectations=out, original_direct_mass=mass,
                original_3_minus_5=mass*Q(1093,7776))


def common_collar(cone_mass):
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    early = json.loads((INPUT/'en11/ledger.json').read_text())
    later = json.loads((INPUT/'en11b/applications.json').read_text())
    weights = {}
    for eta in ('1/8','1/16','1/32'):
        weight = sum(Q(c['covered_at_radius'][eta]) for f in early['faces'] for c in f['exits'])
        weight += sum(Q(row['second_exit'][eta]) for row in later)
        # Round down in exact arithmetic to keep the certificate compact.
        weights[eta] = Q(weight*10**12//1,10**12)
    eta = Q(1,8)
    width = Q(16,3)*eta
    bound = 1-cone_mass-(1-width)*(weights['1/8']-cone_mass)
    return dict(mass_lower_bounds=weights, remainder_upper_bound=bound)


def main():
    work = subprocess.check_output([sys.executable,'-B',str(INPUT/'r2_routes/control_weights.py')])
    assert work == (INPUT/'r2_routes/control-weights.json').read_bytes()
    cone = exact_cone()
    result = dict(box=json.loads(work)['clock_marginal_class5_minus11'],
                  prefix=check_prefix(), cone=cone, small_cube_ramp=boxed_threshold(),
                  common_second_exit_collar=common_collar(cone['mass']))
    print(json.dumps(result,indent=2,default=str))


if __name__ == '__main__':
    main()
