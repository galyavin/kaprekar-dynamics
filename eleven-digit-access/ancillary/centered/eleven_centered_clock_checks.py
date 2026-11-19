"""Independent replay, symbolic finite-base bridge, and sharp-scale checks."""
from itertools import product
from math import comb, factorial
from pathlib import Path
from random import Random
import hashlib
import json
import subprocess
import sys
import sympy as sp
import eleven_centered_clock as m
import eleven_centered_clock_audit as audit
from eleven_two_cycle_clock import offsets
from eleven_two_cycle_clock_checks import whole_step, nonnegative_on_family
from odd_base_covariance import in_region, centered
from odd_gap import subtract_pair

ROOT = Path(__file__).resolve().parent
require = m.require


def branch(x):
    s = x[0]+x[4]
    return [0 if s <= -7 else 1 if s <= 5 else 2,
            int(x[1]+x[3] >= -6), int(x[2] >= -3)]


def witness_replay(cert):
    count = steps = 0
    depths = []
    for case in audit.inventory(cert):
        nodes = case['nodes']
        paths = [()]
        for ni, node in enumerate(nodes[1:], 1):
            path = paths[node['parent']]+(ni,)
            paths.append(path)
            x = tuple(node['witness'])
            beta = m.base_class(x)
            minimum = min(m.gaps(x))
            require(minimum < -m.THRESHOLD and m.gaps(x)[case['index']] == minimum,
                    'Recorded source realizes its source minimum')
            for depth, ci in enumerate(path, 1):
                current = nodes[ci]
                require(branch(x) == current['branch'], 'Witness branch replay')
                x = m.centered_step(x)
                require(m.base_class(x) == beta, 'Preserved physical base class')
                gs = m.gaps(x)
                require(gs.index(min(gs)) == current['minimum'], 'Earliest minimum convention')
                require(minimum-m.GROWTH*depth <= min(gs) <= minimum,
                        'Surviving path and one-step growth bounds')
                steps += 1
            count += 1
        depths.append(max(map(len, paths)))
    return dict(physical_source_witnesses=count, primitive_replay_steps=steps,
                maximum_depths=depths)


def algebra():
    # Every affine sorting cell preserves the same physical coset.
    z = sp.symbols('z0:5', integer=True)
    identities = 0
    for beta in range(6):
        x = tuple(6*a-beta*v for a, v in zip(z, m.V))
        for label in product(range(3), range(2), range(2)):
            forms, _ = audit.successor(audit.INITIAL, list(label), 0)
            y = [sp.expand(sum(a*b for a, b in zip(f[:5], x))+f[5]) for f in forms]
            for i, f in enumerate(y, 1):
                require(all(int(c) % 6 == 0 for c in sp.Poly(f-i*beta, *z).coeffs()),
                        'Universal physical-coset preservation on each affine cell')
                identities += 1
    # Write all actual gaps as t+h_i, where h_i>=0. The sufficient
    # paired condition 7t>=B+1 becomes t>=sum(h)+1.
    h = sp.symbols('h0:6', nonnegative=True)
    u = sp.Symbol('u', nonnegative=True)
    t = sum(h)+1+u
    gs = [t+a for a in h]
    B = sum(gs)
    q = [sum(gs[i+1:]) for i in range(5)]
    pairs = [(q[i]-int(i == 4), B-q[4-i]-int(i < 4)) for i in range(5)]
    guards = [B-1-q[0], q[-1]-1]
    guards += [a-b for left, right in zip(pairs, pairs[1:]) for a in left for b in right]
    for f in guards:
        require(all(c >= 0 for c in sp.Poly(sp.expand(f), *h, u).coeffs()),
                'Sufficient paired-region condition for every real gap vector')
        identities += 1
    # Coercivity and state count: centered gaps sum to zero and all have
    # residue -beta. They determine the five coordinates uniquely.
    totals = []
    for beta in range(6):
        lower = (-m.ABSORBING+beta+5)//6
        total = beta-6*lower
        require(total >= 0 and comb(total+5, 5) <= m.STATE_BOUND, 'Finite state count in every coset')
        totals.append(total)
    require(totals == [504, 499, 500, 501, 502, 503], 'Exact translated simplex sizes')
    require(m.ABSORBING == 120+32*12 == 504, 'Absorbing envelope')
    require(3534 == 7*m.ABSORBING+6, 'Finite-base permanence threshold')
    require(m.STATE_BOUND == 279158015851, 'Explicit, deliberately coarse finite-state constant')
    return dict(affine_identities_and_guards=identities, coset_simplex_totals=totals,
                bounded_states_per_coset=m.STATE_BOUND,
                sufficient_guard='7 min(actual gaps) >= B+1')


def sharp_family():
    t, j = sp.symbols('t j', integer=True)
    guards = identities = 0
    residuals = []
    for beta in range(6):
        x = offsets(18*t+(5*beta)%6-18*j, -108*t+(-2*beta)%6+54*j, (-4, -10, 5))
        initial = x
        rs = []
        for phase in range(8):
            gs = m.gaps(x)
            mi = min(range(6), key=lambda i: gs[i].subs({t: 20, j: 3}))
            r = sp.expand(-gs[mi])
            rs.append(str(r))
            ds = [-6]
            for i in range(5):
                pair = (x[i]-6*int(i == 4), -x[4-i]-6*int(i < 4))
                pair = sorted(pair, key=lambda a: sp.sympify(a).subs({t: 20, j: 3}), reverse=True)
                nonnegative_on_family(pair[0]-pair[1], t, j)
                guards += 1
                ds.extend(pair)
            for g in gs:
                nonnegative_on_family(g-gs[mi], t, j)
                guards += 1
            # At every phase and every 0<=j<t, R>=54t+22.
            nonnegative_on_family(r-54*t-22, t, j)
            guards += 1
            x = tuple(sp.expand(ds[i]-ds[10-i]) for i in range(5))
        for a, b in zip(x, initial):
            require(sp.expand(a-b.subs(j, j+1)) == 0, 'Universal eight-step return')
            identities += 1
        R0 = 126*t+30+(5*beta)%6-(-2*beta)%6
        # B>=6000t+beta; stronger source central guard implies integer
        # min gap >=ceil((B+1)/7)+64. Increasing B only helps.
        nonnegative_on_family(6000*t+beta-7*R0-2694, t, j)
        guards += 1
        residuals.append(rs)
    return dict(unbounded_family_guards=guards, polynomial_return_identities=identities,
                residuals_by_coset_and_phase=residuals,
                lower_capture='8 floor(B/6000), every B>=54000',
                residual_lower_bound='54t+22>504 for t>=9')


def source_weight():
    # The inherited original-input law is Dirichlet(2,...,2) on six
    # gaps. On min(g)>=1/7, write g_i=(1+h_i)/7, sum(h_i)=1.
    symmetric = sp.Rational(factorial(11), 7**11)*sum(
        sp.Rational(comb(6, r), factorial(5+r)) for r in range(7))
    h = sp.symbols('h0:5')
    polynomial = sp.Poly(sp.expand(sp.prod(1+a for a in h)*(2-sum(h))), *h)
    expanded = sp.Rational(factorial(11), 7**11)*sum(
        c*sp.Rational(sp.prod(factorial(a) for a in mon), factorial(5+sum(mon)))
        for mon, c in polynomial.terms())
    require(symmetric == expanded == sp.Rational(805597, 1977326743),
            'Two exact integrations of the original source mass')
    return dict(limiting_original_input_mass=str(symmetric),
                percent=float(100*symmetric), integrations=2,
                scope='initial central region only; no claimed full basin weight')


def lookahead_rank():
    rng = Random(6003)
    sources = 0
    for beta in range(6):
        for scale in (1, 100, 10**12):
            for _ in range(40):
                x = tuple(6*rng.randrange(-scale, scale+1)-beta*v for v in m.V)
                r = m.residual(x)
                w = m.lyapunov(x)
                require(0 <= w <= 32*r, 'Nonnegative finite-lookahead rank bound')
                if r > 120:
                    require(m.lyapunov(m.centered_step(x)) <= w-6, 'One-step rank through arbitrary switching')
                sources += 1
    counter = (270, -576, 600, 294, -324)
    require(m.lyapunov(counter) < 32*m.residual(counter)-11718,
            'Forward growth does not imply a reverse growth bound')
    return dict(physical_sources=sources, base_classes=6,
                rank='min_{0<=j<32}(32 R(Psi^j X)+6j)',
                rejected_reverse_growth_source=counter,
                theorem_dependency='uniform lower-hit assertion; finite controls alone are insufficient')


def actual_controls():
    rng = Random(6001)
    sources = raw = whole = maxcapture = 0
    terminal_r = 0
    for B in (*range(3534, 3540), *range(6000, 6006), *range(54000, 54006)):
        low = (B+7)//7+64
        slack = B-6*low
        for _ in range(12):
            # Six positive gaps with arbitrary integral excesses.
            cuts = sorted(rng.randrange(slack+1) for _ in range(5))
            extras = [b-a for a, b in zip((0, *cuts), (*cuts, slack))]
            gs = [low+a for a in extras]
            q = tuple(sum(gs[i+1:]) for i in range(5))
            require(m.central_source(B, q), 'Declared full-dimensional source region')
            value = 0
            for digit in (*q, 0, 0, 0, 0, 0, 0): value = B*value+digit
            seen = {}
            qs = []
            R0 = m.residual(centered(B, q))
            first_small = None
            while q not in seen:
                seen[q] = len(qs)
                qs.append(q)
                require(in_region(B, q), 'Paired permanence on complete actual trajectory')
                x = centered(B, q)
                r = m.residual(x)
                if r <= m.THRESHOLD and first_small is None: first_small = len(qs)-1
                require(r <= max(R0+384, 504), 'Global envelope')
                if first_small is not None: require(r <= 504, 'Post-threshold envelope')
                predicted = m.centered_step(x)
                q = subtract_pair(B, q)
                value, full = whole_step(B, value)
                require(q == full and centered(B, q) == predicted, 'Whole-integer versus column and centered maps')
                raw += 1
                whole += 1
                require(len(qs) < 100000, 'Control resource cap, never theorem evidence')
            require(first_small is not None and first_small <= m.preparation_bound(centered(B, qs[0])),
                    'Adaptive clock preparation bound')
            capture = seen[q]
            maxcapture = max(maxcapture, capture)
            require(capture <= m.actual_bound(B, qs[0]), 'Declared actual first-pair-capture bound')
            terminal_r = max(terminal_r, *(m.residual(centered(B, p)) for p in qs[capture:]))
            require(terminal_r <= 504, 'Every observed cycle lies in the proved envelope')
            sources += 1
    long_steps = 0
    for B in (*range(54000, 54006), 59999, 60000, 60001):
        q = m.long_family(B)
        require(m.central_source(B, q), 'Long family starts in the central region')
        for j in range(B//6000):
            require(q == m.long_family(B, j), 'Every all-base return source')
            for _ in range(8):
                require(m.residual(centered(B, q)) > 504, 'Every phase is noncyclic by the universal envelope')
                q = subtract_pair(B, q)
                long_steps += 1
        require(q == m.long_family(B, B//6000), 'All-base linear prefix endpoint')
    for B in (10**24+beta for beta in range(6)):
        for j in (0, B//12000, B//6000-1):
            q = m.long_family(B, j)
            for _ in range(8):
                q = subtract_pair(B, q)
                long_steps += 1
            require(q == m.long_family(B, j+1), 'Huge-base return blocks')
    return dict(complete_actual_sources=sources, raw_column_updates=raw,
                independent_whole_integer_updates=whole, maximum_pair_capture=maxcapture,
                maximum_terminal_residual=terminal_r, all_base_family_updates=long_steps,
                maximum_base_digits=len(str(B)))


def negative_controls(cert):
    rejected = 0
    for name, mutate in (
            ('missing_case', lambda c: c['cases'].pop()),
            ('missing_edge', lambda c: c['cases'][0]['nodes'][0]['children'].pop()),
            ('bad_branch', lambda c: c['cases'][0]['nodes'][1]['branch'].__setitem__(0, 3)),
            ('false_horizon', lambda c: c.__setitem__('horizon', 1))):
        bad = json.loads(json.dumps(cert))
        mutate(bad)
        try:
            for case in audit.inventory(bad): audit.expand(case)
        except ArithmeticError: rejected += 1
        else: raise ArithmeticError('Corrupt certificate accepted: '+name)
    try: m.base_class((-64, -132, 410, 352, 70))
    except ValueError: rejected += 1
    else: raise ArithmeticError('Unphysical integral centered vector accepted')
    for script in ('eleven_centered_clock_audit.py', 'eleven_centered_clock_checks.py'):
        p = subprocess.run([sys.executable, '-O', str(ROOT/script)], capture_output=True, text=True)
        require(p.returncode != 0 and 'Optimized Python' in p.stderr, 'Disabled assertions rejected')
        rejected += 1
    return dict(rejected_controls=rejected)


def main():
    if not __debug__: raise RuntimeError('Optimized Python is not accepted')
    cert = json.loads(audit.CERT.read_text())
    out = dict(certificate_sha256=hashlib.sha256(audit.CERT.read_bytes()).hexdigest())
    for name, fn in (('algebra', algebra), ('sharp_family', sharp_family),
                     ('source_weight', source_weight),
                     ('lookahead_rank', lookahead_rank),
                     ('actual_controls', actual_controls),
                     ('witness_replay', lambda: witness_replay(cert)),
                     ('negative_controls', lambda: negative_controls(cert))):
        out[name] = fn()
        print(name, out[name], flush=True)
    # This check report never substitutes for the separate complete audit.
    out['status'] = 'all symbolic and independent controls passed; coverage acceptance is separate'
    (ROOT/'eleven_centered_clock_check_results.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
