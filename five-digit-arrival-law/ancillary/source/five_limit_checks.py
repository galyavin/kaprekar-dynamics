"""Independent finite, geometric, and held-out checks of the five-digit law."""
import argparse
import json
from fractions import Fraction as Q
from math import ceil, isqrt
from pathlib import Path
import subprocess

from five_clock import core_region, core_clock, margins, margin_step, limit_time
from five_digit import finite_graph, subtract, digits, difference_pair, classical_cycles
from five_limit_law import (RHO, cdf, cdf_bounds, integrate, integration_cells,
                            missing_mass, moment_bounds)

ROOT = Path(__file__).resolve().parent
HOLDOUTS = (*range(313, 319), *range(1553, 1559))


def raw_step(base, pair):
    u, v = pair
    representative = u*base**4+v*base**3
    return difference_pair(digits(base, subtract(base, representative)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--holdout-executable')
    args = parser.parse_args()
    pairs = core = 0
    for base in range(30, 91):
        cycles, fate = finite_graph(base, raw_step)
        assert set(cycles) == set(classical_cycles(base).values())
        for pair, (_, h) in fate.items():
            pairs += 1
            assert Q(h+1) <= Q(base, 3)+isqrt(base)+1+20
            if core_region(base, pair):
                core += 1
                following = raw_step(base, pair)
                assert core_region(base, following)
                assert margins(base, following) == margin_step(margins(base, pair))
                assert -3 <= h-core_clock(base, pair) <= 6
    print(f'Direct subtraction: {pairs} pairs; {core} core states; identities and bounds agree.', flush=True)

    # Compare the closed clock against repeated real digit sorting, not its
    # own fold table. Closed drift boundaries are harmless for this identity.
    compared = 0
    for i in range(2, 97):
        for j in range(1, i):
            x, y = Q(i, 97), Q(j, 97)
            u, v = x, y
            for _ in range(200):
                if v <= Q(1, 2) and u-v <= Q(1, 2) and u+v >= 1 and 2*u-v >= 1:
                    break
                ordered = sorted((u, v, Q(1), 1-v, 1-u))
                u, v = ordered[4]-ordered[0], ordered[3]-ordered[1]
            else:
                raise AssertionError(('leading entry failed', x, y))
            assert limit_time(x, y) == Q(2, 3)*max(u+v-1, 2*u-v-1)
            assert 0 <= limit_time(x, y) <= Q(1, 3)
            compared += 1
    print(f'Closed limiting clock equals direct leading digit sorting at {compared} rational points.', flush=True)

    triangle = ((Q(0), Q(0)), (Q(1), Q(0)), (Q(1), Q(1)))
    for p in range(5):
        for q in range(5-p):
            assert integrate(triangle, {(p, q): Q(1)}) == Q(1, (q+1)*(p+q+2))
    assert integrate(triangle, RHO) == 1
    for n in (1, 3, 10, 30):
        assert sum(integrate(poly, RHO) for poly, _, _ in integration_cells(n)) == 1-missing_mass(n)
    thresholds = (Q(1, 25), Q(1, 12), Q(1, 6), Q(1, 4), Q(3, 10), Q(8, 25))
    exact_cdf = {str(t): cdf(t) for t in thresholds}
    for t in thresholds:
        for n in (1, 3, 10):
            lo, hi = cdf_bounds(t, n)
            assert lo <= exact_cdf[str(t)] <= hi
    assert cdf(Q(1, 25)) == Q(713903, 39062500)
    moments = {str(k): moment_bounds(k, 300) for k in (1, 2)}
    for k in (1, 2):
        coarse = moment_bounds(k, 30)
        assert coarse[0] <= moments[str(k)][0] <= moments[str(k)][1] <= coarse[1]
    print('Exact polynomial integrals, normalization, CDF enclosures, and moment refinement agree.', flush=True)

    comparisons = []
    for base in HOLDOUTS:
        path = ROOT/f'five_limit_holdout_{base}.json'
        if args.holdout_executable and not path.exists():
            raw = subprocess.check_output([args.holdout_executable, str(base)], text=True)
            d = json.loads(raw)
            assert d['base'] == base and d['total_input_weight'] == base**5-base
            path.write_text(json.dumps(d, indent=2)+'\n')
        if not path.exists():
            raise FileNotFoundError(f'{path.name}: run with --holdout-executable first')
        d = json.loads(path.read_text())
        hist = {}
        for cycle in d['cycles']:
            for time, weight in cycle['input_arrival_counts'].items():
                hist[int(time)] = hist.get(int(time), 0)+weight
        total = base**5-base
        assert sum(hist.values()) == total == d['total_input_weight']
        assert max(hist) <= Q(base, 3)+ceil(base**0.5)+20
        gaps = {str(t): float(Q(sum(w for h, w in hist.items() if Q(h, base) <= t), total)-exact_cdf[str(t)]) for t in thresholds}
        row = {'base': base, 'mean_over_base': float(Q(d['first_moment_sum'], total*base)),
               'second_moment_over_base_squared': float(Q(d['second_moment_sum'], total*base**2)),
               'maximum_over_base': max(hist)/base,
               'cdf_minus_limit_at_thresholds': gaps}
        comparisons.append(row)
        print(f"Holdout {base}: mean/B={row['mean_over_base']:.9f}, max CDF discrepancy on stated grid={max(map(abs,gaps.values())):.6f}", flush=True)
    output = {'direct_subtraction_pairs': pairs, 'core_states': core,
              'leading_clock_points': compared, 'cdf': {t: str(v) for t, v in exact_cdf.items()},
              'moment_cutoff': 300,
              'moment_intervals': {k: {'rational': [str(v) for v in bounds], 'decimal': [float(v) for v in bounds]} for k, bounds in moments.items()},
              'holdouts': comparisons,
              'scope': 'Held-out bases selected after the clock and law formulas; finite comparisons are corroboration, not convergence-rate bounds.'}
    (ROOT/'five_limit_results.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    main()
