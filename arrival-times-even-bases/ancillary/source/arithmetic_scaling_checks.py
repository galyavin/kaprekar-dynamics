from collections import Counter
from fractions import Fraction as Q
import four_uniform as four
from doubling_checks import multiplicity as weight4


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
    return {'pair_states':pair_total,'regular_pairs':generic_total,'finite_distributions':rows}
