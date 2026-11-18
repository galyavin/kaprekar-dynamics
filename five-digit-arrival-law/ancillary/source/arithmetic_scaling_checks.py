"""Independent finite checks of the uniform and square-root scaling proofs.

Standard library only. The all-base claims are proved in arithmetic_scaling.md;
this script corroborates guards, clock bridges, exact weights and convergence.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import isqrt, sqrt
from pathlib import Path
import json

import four_uniform as four
import five_boundary_layer as layer
from doubling_checks import multiplicity as weight4
from five_digit import multiplicity as weight5


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def raw_subtract(base, digits):
    ds = sorted(digits)
    high = low = 0
    for d in reversed(ds):
        high = base*high+d
    for d in ds:
        low = base*low+d
    value, out = high-low, [0]*len(ds)
    for i in range(len(ds)-1, -1, -1):
        value, out[i] = divmod(value, base)
    return tuple(out)


def raw_pair_step(base, pair, width):
    ds = raw_subtract(base, (*pair, *((0,)*(width-len(pair)))))
    ordered = sorted(ds, reverse=True)
    return tuple(ordered[j]-ordered[-1-j] for j in range(len(pair)))


def graph_distances(nodes, step):
    distances, cycle_states = {}, set()
    for start in nodes:
        if start in distances:
            continue
        path, positions, q = [], {}, start
        while q not in distances and q not in positions:
            positions[q] = len(path)
            path.append(q)
            q = step(q)
        if q in positions:
            cut = positions[q]
            cycle_states.update(path[cut:])
            for p in path[cut:]:
                distances[p] = 0
            path = path[:cut]
        for p in reversed(path):
            distances[p] = distances[step(p)]+1
    return distances, cycle_states


def four_checks():
    rows, pair_total, generic_total = [], 0, 0
    for m in (5, 9, 31, 53):
        for n in range(5):
            base = m*2**n
            if base > 900:
                continue
            pairs = [(u,v) for u in range(1,base) for v in range(u+1)]
            distances, cycles = graph_distances(pairs, lambda q: raw_pair_step(base,q,4))
            hist, excluded = Counter(), 0
            for q, h in distances.items():
                w = weight4(base,q)
                hist[h+1] += w
                if four.regular_pair(base,q):
                    generic_total += 1
                    require(four.generic_pair_arrival(base,q) == h, ('valuation clock',base,q,h))
                    require(four.regular_pair(base,raw_pair_step(base,q,4)), ('invariance',base,q))
                else:
                    excluded += w
            hist[1] -= len(cycles)
            hist[0] += len(cycles)
            total = base**4-base
            require(sum(hist.values()) == total, ('weight',base))
            require(Q(excluded,total) <= Q(9,m)/(1-Q(1,base**3)), ('exceptional mass',base))
            law = four.limiting_law(n)
            tv = sum(abs(Q(hist[t],total)-law.get(t,0)) for t in set(hist)|set(law))/2
            require(tv <= four.error_bound(base), ('TV bound',base,tv))
            require(four.proof_error_budget(base) < Q(26,m), ('proof budget',base))
            pair_total += len(pairs)
            rows.append({'base':base,'odd_part':m,'exponent':n,'pairs':len(pairs),
                         'tv':str(tv),'tv_float':float(tv),'excluded_probability':str(Q(excluded,total))})
    raw_inputs = 0
    for base in (5,7,9,10,14,18):
        nodes = [ds for ds in product(range(base),repeat=4) if len(set(ds))>1]
        distances, _ = graph_distances(nodes,lambda ds:raw_subtract(base,ds))
        for ds,h in distances.items():
            a,b,c,d = sorted(ds,reverse=True)
            if four.regular_pair(base,(a-d,b-c)):
                require(four.generic_actual_arrival(base,ds)==h, ('actual clock',base,ds,h))
                raw_inputs += 1
    return {'pair_states':pair_total,'regular_pairs':generic_total,
            'independent_actual_regular_inputs':raw_inputs,'finite_distributions':rows}


def five_checks():
    rows, starts, high_steps = [], 0, 0
    for base in (199,200,401,402,997,1000):
        # Full pair graph uses raw numeric subtraction, not the high/core formula.
        pairs = [(u,v) for u in range(1,base) for v in range(u+1)]
        distances, _ = graph_distances(pairs,lambda q:raw_pair_step(base,q,5))
        count, max_error = 0,Q(0)
        for r in range(max(2,isqrt(base)//2),isqrt(base)+1):
            for v in range((3*base+4)//5,4*base//5+1):
                q = base-r,v
                try:
                    data = layer.exit_data(base,q)
                except ValueError:
                    continue  # Small bases need not meet the asymptotic band guard.
                h = distances[q]
                require(data['pair_lower'] <= h <= data['pair_upper'], ('layer clock',base,q,h,data))
                require(abs(Q(h+1)-data['center_clock']) <= 10, ('uniform O(1)',base,q,h,data))
                # Validate every high transition for a deterministic subcollection.
                if v % 17 == 0:
                    p = q
                    for _ in range(data['high_steps']):
                        require(raw_pair_step(base,p,5)==(p[0]-1,p[0]+p[1]-base+1),('high branch',base,p))
                        p = raw_pair_step(base,p,5)
                        high_steps += 1
                max_error = max(max_error,abs(Q(h+1)-data['center_clock']))
                count += 1
        starts += count
        rows.append({'base':base,'checked_starts':count,'max_actual_center_error':str(max_error)})
    # Weighted overshoot sums on complete and truncated bands. These are finite
    # convergence observations, not a proof or independent replay of the core.
    bands = []
    a,b,c,d = Q(1,2),Q(1),Q(3,5),Q(4,5)
    mass = layer.band_mass_constant(a,b,c,d)
    for base in (10000,40000,160000):
        root = isqrt(base)
        total_weight, theta_sum, theta2_sum = 0,0.0,0.0
        # Aggregate complete overshoot blocks, summing their polynomial weights
        # analytically; endpoint blocks use the same power-sum identities.
        for r in range(root//2,root+1):
            vmin,vmax = (3*base+4)//5,4*base//5
            v = vmin
            while v <= vmax:
                data = layer.exit_data(base,(base-r,v))
                n,R = data['high_steps'],data['R']
                H = (base+1)//2
                D = n*(2*r+n-3)//2
                end = min(vmax,H+D)
                lo,hi = H+D-end,H+D-v
                sums = [sum_powers(hi,k)-sum_powers(lo-1,k) for k in range(5)]
                A = H+D
                coefficients = [r*(120*A*(base-r-A)-20),120*r*(2*A-base+r),-120*r]
                W = sum(coefficients[k]*sums[k] for k in range(3))
                P = sum(coefficients[k]*sums[k+1] for k in range(3))
                P2 = sum(coefficients[k]*sums[k+2] for k in range(3))
                total_weight += W
                theta_sum += P/(R-2)
                theta2_sum += P2/(R-2)**2
                v = end+1
        bands.append({'base':base,'B_times_band_probability':float(Q(base*total_weight,base**5-base)),
                      'target_mass':str(mass),'phase_mean':theta_sum/total_weight,
                      'phase_second_moment':theta2_sum/total_weight})
    return {'exact_clock_starts':starts,'raw_high_steps':high_steps,'clock_checks':rows,
            'weighted_phase_checks':bands}


def sum_powers(n,k):
    if n < 0:
        return 0
    if k==0:
        return n+1
    if k==1:
        return n*(n+1)//2
    if k==2:
        return n*(n+1)*(2*n+1)//6
    if k==3:
        return (n*(n+1)//2)**2
    if k==4:
        return n*(n+1)*(2*n+1)*(3*n*n+3*n-1)//30
    raise ValueError(k)


def main():
    result = {'four':four_checks(),'five':five_checks(),
              'scope':'Finite direct-digit corroboration; all-base proofs are analytic.'}
    root = Path(__file__).resolve().parent
    (root/'arithmetic_scaling_results.json').write_text(json.dumps(result,indent=2)+'\n')
    lines = ['Uniform arithmetic and second-scale checks passed.',
             f"Four-slot pair states: {result['four']['pair_states']}; regular: {result['four']['regular_pairs']}.",
             f"Independent actual regular four-slot inputs: {result['four']['independent_actual_regular_inputs']}.",
             f"Five-slot exact clock starts: {result['five']['exact_clock_starts']}; raw high steps: {result['five']['raw_high_steps']}.",
             'All finite claims are corroboration; no asymptotic claim is inferred from observed errors.']
    (root/'arithmetic_scaling_verification.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    print(json.dumps(result['five']['weighted_phase_checks'],indent=2))


if __name__=='__main__':
    main()
