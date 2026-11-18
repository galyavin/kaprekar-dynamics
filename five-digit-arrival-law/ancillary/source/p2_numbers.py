"""Exact manuscript constants and independent small-base input multiplicities."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path
from five_digit import multiplicity
from five_limit_law import cdf, moment_bounds
counts=0
for B in range(2,9):
    hist=Counter()
    for digits in product(range(B),repeat=5):
        ds=sorted(digits)
        if ds[0] != ds[-1]:
            hist[(ds[4]-ds[0],ds[3]-ds[1])] += 1
            counts+=1
    assert sum(hist.values()) == B**5-B
    assert all(value==multiplicity(B,pair) for pair,value in hist.items())
F=cdf(Q(1,25))
assert F==Q(713903,39062500)
assert 1-F==Q(38348597,39062500)==Q('0.9817240832')
lo,hi=moment_bounds(1,300)
assert Q('0.1827011465') < lo <= hi < Q('0.1827011669')
# Constants in the geometric finite-base CDF bound, after factoring out 1/B.
for N in range(1,1001):
    e=2**(N+5)
    guard=27000*(N+1)*(N+4)*e
    clock=22*96*(N+1)*(N+10+2*e)
    assert guard+clock <= 10**5*(N+4)**2*(e+N+10)
result={'cdf':str(F),'tail':str(1-F),'tail_decimal':'0.9817240832',
        'mean_bounds':[str(lo),str(hi)],'mean_decimal_bounds':[float(lo),float(hi)],
        'ordered_inputs_checked':counts,'geometric_constant_checks':1000}
Path('p2_numbers_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
