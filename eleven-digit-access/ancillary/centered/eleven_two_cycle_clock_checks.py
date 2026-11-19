"""Independent algebra, guarded raw subtraction, and sharp-scale controls."""
from fractions import Fraction
from itertools import product
from pathlib import Path
from random import Random
import json
import subprocess
import sys
import sympy as sp
import eleven_two_cycle_clock as m
from odd_base_covariance import centered_step
from odd_gap import subtract_pair

ROOT = Path(__file__).resolve().parent
require = m.require


def algebra():
    M = sp.Matrix(m.M)
    c = sp.Matrix(m.C)
    u = sp.Matrix(sp.symbols('p q r'))
    center = sp.Matrix([sp.Rational(-7, 2), -9, 5])
    require(M**3 == -sp.eye(3) and M**6 == sp.eye(3), 'Finite-order phase matrix')
    require(M*center+c == center, 'Rational affine center')
    require(sum((M**j for j in range(6)), sp.zeros(3)) == sp.zeros(3), 'Zero mean of centered phase')
    v = M*u+c
    da, dd = -u[0]-u[2], u[0]-u[1]-1
    require(sp.expand(da+sp.Rational(3, 2)-(v[0]-u[0])/2) == 0, 'Exact a coboundary')
    require(sp.expand(dd-sp.Rational(9, 2)+3*(v[0]-u[0])/2-(v[1]-u[1])) == 0,
            'Exact d coboundary')
    require(sp.expand(-2*da+v[0]-u[0]) == 3, 'Clock ignores the alternating phase')
    forms = []
    v = u
    drift = sp.zeros(2, 1)
    for _ in range(6):
        p, q, r = v
        forms.extend((r-1, p-q-2, p+r+2, p+2*r-1))
        drift += sp.Matrix([-p-r, p-q-1])
        v = M*v+c
    require(sp.simplify(v-u) == sp.zeros(3, 1), 'Affine sixth iterate, including constants')
    require(sp.simplify(drift) == sp.Matrix([-9, 27]), 'Universal longitudinal displacement')
    forms = set(map(sp.expand, forms))
    p, q, r = u
    require({p+9, 2-p, r-1, 9-r, p-q-2, -p+q+9} <= forms, 'Derived finite polytope box')
    require(all(f.subs(dict(zip(u, center))) > 0 for f in forms), 'Nonempty real interior')
    return dict(distinct_affine_guards=len(forms), matrix_power=6,
                center=list(map(str, center)), twelve_phase_return_drift=[-9, 27])


def phases():
    count = 0
    for u in product(range(-10, 11), repeat=3):
        v, (da, dd) = m.phase_step(u)
        p, q, r = u
        require(v[0] <= 2, 'Prepared phase bound')
        require(dd+3*da-(v[1]-q) == max(q-p-max(-1, -r)-1, 0), 'Transverse defect identity')
        a = 10000
        d = -60000
        x = m.offsets(a, d, u)
        y = centered_step(x)
        require(x[1]+x[3] <= -6 and x[2] <= -3 and y[1]+y[3] >= -6 and y[2] >= -3,
                'Independent primitive sign guards')
        require(centered_step(y) == m.offsets(a+6*da, d+6*dd, v), 'Independent two-step primitive')
        s = m.second_shift(u, 4)
        require(s >= -2, 'Four-return lower bound')
        if p <= 2 and s <= 0:
            require(-4 <= p <= 2 and -10 <= q <= 4 and -4 <= r <= 1, 'Small prepared exception box')
        require(m.second_shift(u, 16) >= 1, 'Universal clock control')
        count += 1
    pts = [u for u in product(range(-9, 3), range(-18, 1), range(1, 10)) if m.in_phase_polytope(u)]
    for u in pts:
        v, drift = m.phase_prefix(u, 12)
        require(v == u and drift == (-9, 27), 'All physical points of the entire polytope')
        require(m.in_phase_polytope(m.affine_phase(u)), 'Polytope invariance')
        require(m.second_shift(u, 2) == 3, 'Exact four-update clock')
    rng = Random(5801)
    for _ in range(1000):
        u = tuple(rng.randrange(-10**12, 10**12) for _ in range(3))
        require(m.second_shift(u, 16) >= 1, 'Large independent integer phases')
    for _ in range(100):
        a, b = rng.sample(pts, 2)
        t = Fraction(rng.randrange(101), 100)
        u = tuple(t*x+(1-t)*y for x, y in zip(a, b))
        require(m.in_phase_polytope(u), 'Exact rational convex controls')
        require(m.phase_prefix(u, 12) == (u, (-9, 27)), 'Real affine-region displacement')
    return dict(integer_phase_sources=count, primitive_updates=2*count,
                complete_polytope_lattice_points=len(pts), large_integer_sources=1000,
                rational_polytope_sources=100)


def whole_step(B, value):
    ds = []
    for _ in range(11):
        value, digit = divmod(value, B)
        ds.append(digit)
    require(value == 0, 'Eleven-slot whole-number control')
    ds.sort()
    lo = hi = 0
    for z in ds: lo = B*lo+z
    for z in reversed(ds): hi = B*hi+z
    out = hi-lo
    v = out
    digits = []
    for _ in range(11):
        v, digit = divmod(v, B)
        digits.append(digit)
    digits.sort(reverse=True)
    q = tuple(digits[i]-digits[-1-i] for i in range(5))
    return out, q


def actual():
    rng = Random(5802)
    count = steps = 0
    for beta in range(6):
        for scale in (1, 10**6, 10**30):
            for _ in range(8):
                u = tuple(rng.randrange(-100*scale, 100*scale+1) for _ in range(3))
                A = 10**6*scale
                B = 6000*A+beta
                a = 6*A+(5*beta)%6
                d = -36*A+(-2*beta)%6
                q = m.pair(B, a, d, u)
                require(m.coordinates(B, q) == (a, d, u), 'Physical coordinates in every base class')
                initial = q
                value = 0
                for digit in (*q, 0, 0, 0, 0, 0, 0): value = B*value+digit
                for _ in range(16):
                    require(m.in_regime(B, q), 'Actual paired and sign guards for both updates')
                    predicted = m.two_step(B, q)
                    for _ in range(2):
                        raw = subtract_pair(B, q)
                        value, full = whole_step(B, value)
                        require(raw == full, 'Whole-number subtraction versus independent digit columns')
                        q = raw
                        steps += 1
                    require(q == predicted, 'Phase return versus actual whole-integer result')
                require(q[1]-initial[1] == m.second_shift(u, 16) >= 1, 'Physical thirty-two-update clock')
                count += 1
    return dict(sources=count, whole_integer_updates=steps, digit_column_updates=steps,
                maximum_base_digits=len(str(B)), base_classes=6)


def nonnegative_on_family(f, t, j):
    f = sp.Poly(sp.expand(f), t, j)
    require(f.total_degree() <= 1, 'Affine family guard')
    a = f.coeff_monomial(t)
    b = f.coeff_monomial(j)
    c = f.coeff_monomial(1)
    # The domain t>=2, 0<=j<=t-1 has vertices (2,0),(2,1)
    # and recession rays (1,0),(1,1).
    require(min(a, a+b, 2*a+c, 2*a+b+c) >= 0,
            'Guard on every t>=2, 0<=j<=t-1: '+str(f.as_expr()))


def long_family():
    t, j = sp.symbols('t j', integer=True)
    guards = identities = 0
    for beta in range(6):
        B = 600*t+beta
        a = 18*t+(5*beta)%6-18*j
        d = -108*t+(-2*beta)%6+54*j
        q = tuple(sp.expand((B*v+x)/6) for v, x in zip(m.V, m.offsets(a, d, (-4, -10, 5))))
        initial = q
        for h in range(8):
            valid = (B-1-q[0], q[-1]-1, *(a-b for a, b in zip(q, q[1:])))
            bs = tuple((q[i]-int(i == 4), B-q[4-i]-int(i < 4)) for i in range(5))
            paired = tuple(a-b for left, right in zip(bs, bs[1:]) for a in left for b in right)
            x = tuple(6*z-B*v for z, v in zip(q, m.V))
            signs = (-6-x[1]-x[3], -3-x[2]) if h % 2 == 0 else (x[1]+x[3]+6, x[2]+3)
            ds = (*q[:-1], q[-1]-1, B-1, *(B-z-1 for z in reversed(q[1:])), B-q[0])
            order = sorted(range(11), key=lambda i: ds[i].subs({t: 10, j: 3}), reverse=True)
            sorted_ds = tuple(ds[i] for i in order)
            sort_guards = tuple(a-b for a, b in zip(sorted_ds, sorted_ds[1:]))
            for f in (*valid, *paired, *signs, *sort_guards):
                nonnegative_on_family(f, t, j)
                guards += 1
            q = tuple(sp.expand(sorted_ds[i]-sorted_ds[-1-i]) for i in range(5))
        expected = tuple(f.subs(j, j+1) for f in initial)
        for a, b in zip(q, expected):
            require(sp.expand(a-b) == 0, 'Eight actual updates advance the family by one block')
            identities += 1
    steps = 0
    for n in (2, 3, 7, 31, 101):
        base, q = m.long_family(n)
        for k in range(n):
            require((base, q) == m.long_family(n, k), 'Long family block source')
            for h in range(8):
                if h % 2 == 0: require(m.in_regime(base, q), 'Entire long family remains in the regime')
                q = subtract_pair(base, q)
                steps += 1
        require((base, q) == m.long_family(n, n), 'Long family endpoint')
    for base in (*range(1200, 1206), *range(1794, 1806)):
        q = m.long_family_base(base)
        for k in range(base//600):
            require(q == m.long_family_base(base, k), 'All-base lifted family source')
            for h in range(8):
                if h % 2 == 0: require(m.in_regime(base, q), 'Every base-class lift remains guarded')
                q = subtract_pair(base, q)
                steps += 1
        require(q == m.long_family_base(base, base//600), 'All-base lifted family endpoint')
    return dict(unbounded_affine_guards=guards, return_polynomial_identities=identities,
                additional_raw_updates=steps, longest_guarded_primitive_prefix=808,
                scale='8 floor(B/600) updates in every B>=1200')


def negative_controls():
    rng = Random(5803)
    counter = None
    for _ in range(1000):
        B = 1000
        q = tuple(sorted((rng.randrange(1, B) for _ in range(5)), reverse=True))
        if m.in_regime(B, q): continue
        a, d, u = m.coordinates(B, q)
        v, (da, dd) = m.phase_step(u)
        predicted = m.pair(B, a+6*da, d+6*dd, v)
        actual = subtract_pair(B, subtract_pair(B, q))
        if predicted != actual:
            counter = dict(base=B, source=q, unguarded_prediction=predicted, actual=actual)
            try: m.two_step(B, q)
            except ValueError: pass
            else: raise ArithmeticError('Unguarded actual return was accepted')
            break
    require(counter is not None, 'Omitted guards have an actual counterexample')
    require(m.second_shift((6, 5, -6), 12) == -2, 'A shorter universal clock is false')
    B = 6000000
    source = m.pair(B, 6000, -36000, (6, 5, -6))
    q = source
    for _ in range(12):
        require(m.in_regime(B, q), 'Short-clock counterexample has actual guards')
        q = subtract_pair(B, subtract_pair(B, q))
    require(q[1] == source[1]-2, 'Actual twenty-four-update clock counterexample')
    require(not m.in_regime(600, (497, 390, 281, 182, 108)), 'The unrounded t=1 family fails a sign guard')
    rejected = 0
    for name in ('eleven_two_cycle_clock_audit.py', 'eleven_two_cycle_clock_checks.py'):
        r = subprocess.run([sys.executable, '-O', str(ROOT/name)], text=True, capture_output=True)
        require(r.returncode != 0 and 'Optimized Python is not accepted' in r.stderr, 'Disabled assertions rejected')
        rejected += 1
    return dict(omitted_guard_counterexample=counter, false_twelve_return_increment=-2,
                actual_short_clock_counterexample=dict(base=B, source=source, image=q),
                rejected_small_family_base=600,
                optimized_programs_rejected=rejected)


def main():
    if not __debug__: raise RuntimeError('Optimized Python is not accepted')
    out = {}
    for f in (algebra, phases, actual, long_family, negative_controls):
        out[f.__name__] = f()
        print(f.__name__, out[f.__name__], flush=True)
    import eleven_two_cycle_clock_audit as audit
    records = json.loads((ROOT/'eleven_two_cycle_clock_audit_results.json').read_text())
    audit.accept(records['records'])
    out['accepted_unbounded_lemmas'] = len(records['records'])
    out['status'] = 'all independent controls passed; unbounded audit accepted'
    (ROOT/'eleven_two_cycle_clock_results.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
