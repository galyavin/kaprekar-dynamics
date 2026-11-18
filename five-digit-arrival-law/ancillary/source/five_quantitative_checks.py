"""Check the all-base witness construction against actual digit subtraction."""
from fractions import Fraction as Q
import json
from math import isqrt, sqrt
from pathlib import Path

from five_clock import core_region, margins
from five_limit_checks import raw_step
from five_quantitative import maximum_witness, finite_cdf_error_bound


def main():
    bases = [*range(30, 251), *range(1000, 1011), *range(1973, 1979),
             100003, 100004, 100005]
    rows = []
    for base in bases:
        start, n = maximum_witness(base)
        q = start
        for j in range(n):
            assert q == (base-1-j, start[1]-j*(j-1)//2)
            assert q[0] > q[1] and 2*q[1] > base+1
            q = raw_step(base, q)
        assert q == (base-1-n, (base-1)//2)
        assert core_region(base, q)
        assert max(margins(base, q)) == q[1]-1-n
        seen = {}
        q = start
        while q not in seen:
            seen[q] = len(seen)
            q = raw_step(base, q)
        actual = seen[q]+1  # The starting pair is in the high region.
        assert actual >= base/3+sqrt(base)/3-5
        assert actual <= base/3+sqrt(base)/3+60
        assert n >= sqrt(base)-2
        rows.append({'base': base, 'pair': start, 'high_steps': n,
                     'actual_time': actual,
                     'centered_time': actual-base/3-sqrt(base)/3})
    # Complete older references corroborate the maximum upper/lower theorem.
    root = Path(__file__).resolve().parent
    references = []
    for p in sorted(root.glob('five_joint_holdout_*.json')):
        d = json.loads(p.read_text()); b = d['base']
        maximum = max(int(t) for c in d['cycles'] for t in c['input_arrival_counts'])
        assert b/3+sqrt(b)/3-5 <= maximum <= b/3+sqrt(b)/3+60
        references.append({'base': b, 'maximum': maximum})
    rates = {str(2**e): str(finite_cdf_error_bound(2**e)) for e in (40, 80, 120)}
    assert Q(rates[str(2**120)]) < Q(rates[str(2**80)]) < 1
    result = {'witness_checks': len(rows), 'witnesses': rows,
              'complete_reference_maxima': references,
              'conservative_cdf_bounds': rates,
              'scope': 'Witness and finite-reference checks corroborate the analytic '
                       'all-base maximum theorem; the CDF rate is proved geometrically.'}
    (root/'five_quantitative_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f'Passed {len(rows)} direct-subtraction witness checks and '
          f'{len(references)} complete-reference maxima; conservative rate tends down.')


if __name__ == '__main__':
    main()
