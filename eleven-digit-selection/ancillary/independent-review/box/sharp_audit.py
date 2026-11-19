"""Replay the residual routes, verify guard extrema, weigh fibers and audit capture.

Run from any directory with Python 3.13. Only temporary snapshots are written.
"""
from collections import Counter, defaultdict
from contextlib import redirect_stdout
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import product
from math import ceil, floor, gcd, lcm
from pathlib import Path
import fcntl
import io
import json
import os
import shutil
import sys
import tempfile

BOX_INPUTS = Path(__file__).resolve().parents[2] / 'box'
INPUTS = (
    'verify.py', 'certificate.py', 'control.py', 'quotient.py', 'routes.json',
    'quotient.json', 'channel.py', 'endpoint_verification.json', 'cohort1.py',
    'accepted_quotient.json',
)
V = (5,4,3,2,1)
ALIGN = ((-5,0,1,2,-1), (0,0,0,0,0), (-5,2,3,4,-1))
BOX = ((94040,94130),(73770,73850),(57940,58015),(20,80),(200,300))


def chart(B, x, c0, z, r1):
    a,b,c,h,g = x
    n = (h + 1 - (-4,-2,-6)[c0]) / Q(6)
    al = ALIGN[c0]
    m = (g + 5*n - al[2] - al[3] - 1 + z) / 2
    N = (a+c-8*n+al[0]+al[2]-10*m-3*(B-1)/Q(2)-r1) / 22
    return c-3*n+al[2]-2*N, b-4*n+al[1]-2*m-4*N


@cache
def first_cycle(point):
    # A bounded centered trajectory can be realized at this reference base.
    B = 6*10000+5
    def decode(x):
        assert all((z+B*v)%6==0 for z,v in zip(x,V))
        return tuple((z+B*v)//6 for z,v in zip(x,V))
    q = decode(point)
    seen = {}
    while q not in seen:
        a,b,c,d,e = q
        assert B>a>b>c>d>e>0
        digits = (a,b,c,d,e-1,B-1,B-e-1,B-d-1,B-c-1,B-b-1,B-a)
        groups = ((digits[5],), *((digits[i],digits[i+6]) for i in range(5)))
        assert all(min(l)>=max(r) for l,r in zip(groups,groups[1:]))
        rows = sorted(digits, reverse=True)
        seen[q] = len(seen)
        q = tuple(rows[i]-rows[-1-i] for i in range(5))
    return seen[q]


@cache
def clock_corrections(center, drift, length):
    # Parameterize the entire physical segment by any nonconstant coordinate.
    pivot = max(range(5), key=lambda i:abs(drift[i]))
    low, high = sorted((center[pivot],center[pivot]-drift[pivot]))
    corrections = []
    for xi in range(ceil(low),floor(high)+1):
        f = (center[pivot]-xi)/drift[pivot]
        q = tuple(z-f*d for z,d in zip(center,drift))
        if not 0<=f<1 or any(z.denominator!=1 or (z+5*v)%6 for z,v in zip(q,V)):
            continue
        point = tuple(map(int,q))
        for back in range(5042):
            before = tuple(z-back*d for z,d in zip(point,drift))
            mu = first_cycle(before)
            if mu:
                corrections.append((1+mu-length*(f+back),back))
                break
        else:
            raise AssertionError(('no transient predecessor', point))
    assert corrections
    return min(c[0] for c in corrections), max(c[0] for c in corrections), max(c[1] for c in corrections)


def run():
    import verify as cert
    # Same route generator, independently derived chart and full-vertex minima.
    vertices = []
    for x in product(*BOX):
        coordinates = [Q(y,100000) for y in x]
        p,r = chart(1, coordinates, 1, 1, 0)
        p0,r0 = chart(0, (0,)*5, 1, 1, 0)
        vertices.append((p-p0,r-r0))
    original_extrema = cert.coefficient_bounds
    @cache
    def extrema(linear):
        a,b,c = linear
        vals = [a+b*p+c*r for p,r in vertices]
        result = (min(vals)/cert.SCALE, max(vals)/cert.SCALE)
        assert result == original_extrema(linear)
        return result
    cert.coefficient_bounds = extrema
    original_sections = cert.initial_sections
    def sections():
        for label, guards, *rest in original_sections():
            expected = chart(0,(0,)*5,*label[:3])
            assert expected == (guards.p0,guards.r0)
            yield (label,guards,*rest)
    cert.initial_sections = sections
    # Recover leaf lattice fibers and tighten the first-entrance bound.
    indices = Counter()
    slices = defaultdict(Counter)
    terminal = [None,None,0]
    original_terminal = cert.terminal_time_bounds
    def terminal_bounds(state,word,drift,checks,limits,control,clock,guards):
        res = original_terminal(state,word,drift,checks,limits,control,clock,guards)
        if res is None:
            return res
        end,offset,_ = res
        center = tuple(row.constant() for row in cert.shift(state,control-offset,drift))
        low,high,back = clock_corrections(center,tuple(drift),len(word))
        intercept = guards.parts(end)[2]
        low,high = intercept+low,intercept+high
        terminal[0] = min(terminal[0],low) if terminal[0] is not None else low
        terminal[1] = max(terminal[1],high) if terminal[1] is not None else high
        terminal[2] = max(terminal[2],back)
        return res
    cert.terminal_time_bounds = terminal_bounds
    original_branches = cert.branches
    def branches(entry,clock,section,lattice,guards):
        label,k,arm,rho = guards.current
        for leaf in original_branches(entry,clock,section,lattice,guards):
            end,_,_,lat = leaf
            g = gcd(*(col[0] for col in lat.columns))
            D = abs(cert.det(*lat.columns))
            assert D%g==0
            indices[g] += 1
            # In the (p,r) plane the primary class has determinant 216.
            key = (lat.point[0]%g,g,tuple(end.leading().record()[:3]))
            slices[(label,k,arm,rho)][key] += Q(216*g,D)
            yield leaf
    cert.branches = branches
    output = io.StringIO()
    with redirect_stdout(output):
        cert.main()
    data = output.getvalue().encode()
    assert data == (BOX_INPUTS/'result.json').read_bytes()
    result = json.loads(data)
    modulus = lcm(*indices)
    rows = {}
    per_arm = {}
    hist = Counter()
    for (label,k,arm,rho), weights in slices.items():
        for B in range(5+6*(k+arm),modulus,36):
            row = Counter()
            for (b,g,p), weight in weights.items():
                if B%g == b:
                    row[p] += weight
            assert sum(row.values())==1
            rows[B,label,rho] = row
            hist.update(row.values())
            per_arm.setdefault(B,Counter()).update(row)
    assert all(sum(v.values())==1056 for v in per_arm.values())
    comparisons = Counter()
    for (B,label,rho),row in rows.items():
        c,z,r1,r2 = label
        other = (c,z,r1,601+(r2-601+36)%72)
        comparisons[row==rows[(B+12)%modulus,other,rho]] += 1
    equal_by_arm = all(per_arm[B]==per_arm[(B+12)%modulus] for B in per_arm)
    return dict(pinned_result_sha256=sha256(data).hexdigest(),entries=result['completed'],
        leaves=result['leaves'],guard_cutoff=result['minimum_valid_base'],
        distinct_guard_linear_parts=extrema.cache_info().currsize,
        fractional_clock_errors=[str(terminal[0]),str(terminal[1])],
        maximum_backtracked_returns=terminal[2], bounded_tail_states_tested=first_cycle.cache_info().currsize,
        fiber_base_steps=sorted(indices.items()), fiber_modulus=modulus,
        fiber_profile_weight_histogram=sorted((str(k),v) for k,v in hist.items()),
        weighted_r2_shift_matches=comparisons[True], weighted_r2_shift_mismatches=comparisons[False],
        weighted_arm_totals_period12=equal_by_arm,
        arm_profile_masses_in_primary_units={str(B):[[list(p),str(w)] for p,w in sorted(per_arm[B].items())]
                                          for B in (5,11)})


if __name__=='__main__':
    with tempfile.TemporaryDirectory(prefix='box-audit-') as temp:
        for name in INPUTS:
            shutil.copy2(BOX_INPUTS / name, Path(temp) / name)
        os.chdir(temp)
        sys.path.insert(0,temp)
        lock_path = Path(os.environ.get('KAPREKAR_COMPUTE_LOCK', str(Path(temp) / 'compute.lock')))
        with lock_path.open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            summary = run()
    print(json.dumps(summary,indent=2))
