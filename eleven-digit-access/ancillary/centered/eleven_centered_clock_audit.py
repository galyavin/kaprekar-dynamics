"""Independent Z3/Yices coverage audit of physical centered minimum trees.

Reconstructs each affine successor by sorting the five original offset
pairs. It does not import the discovery branch formulas or accept saved
discovery UNSAT results. Complete source inventory and all coverage nodes
are mandatory; per-batch output is only a checkpoint.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import time
import z3
import eleven_centered_clock as m

ROOT = Path(__file__).resolve().parent
CERT = ROOT/'eleven_centered_clock_certificates.json'
CHECKPOINT = ROOT/'audit-checkpoints'
ZERO = (0,)*6
ONE = (0, 0, 0, 0, 0, 1)
INITIAL = tuple(tuple(int(i == j) for j in range(6)) for i in range(5))
require = m.require
_worker_certificate = None
_worker_case = None


def add(*xs): return tuple(map(sum, zip(*xs)))
def mul(c, x): return tuple(c*a for a in x)
def sub(a, b): return add(a, mul(-1, b))
def affine_gaps(x): return (mul(-1, x[0]), *[sub(a, b) for a, b in zip(x, x[1:])], x[-1])


def successor(x, branch, minimum):
    require(isinstance(branch, list) and len(branch) == 3 and
            all(type(i) is int for i in branch) and branch[0] in (0, 1, 2) and
            branch[1] in (0, 1) and branch[2] in (0, 1), 'Exact affine branch label')
    require(type(minimum) is int and 0 <= minimum < 6, 'Exact minimum label')
    outer, inner, middle = branch
    a = add(x[0], x[4])
    b = add(x[1], x[3])
    guards = ([sub(mul(-7, ONE), a)] if outer == 0 else
              [add(a, mul(6, ONE)), sub(mul(5, ONE), a)] if outer == 1 else
              [sub(a, mul(6, ONE))])
    guards.append(add(b, mul(6, ONE)) if inner else sub(mul(-7, ONE), b))
    guards.append(add(x[2], mul(3, ONE)) if middle else sub(mul(-4, ONE), x[2]))
    flags = (outer != 0, bool(inner), bool(middle), bool(inner), outer == 2)
    digits = [mul(-6, ONE)]
    for i, flag in enumerate(flags):
        left = sub(x[i], mul(6*int(i == 4), ONE))
        right = sub(mul(-1, x[4-i]), mul(6*int(i < 4), ONE))
        digits.extend((left, right) if flag else (right, left))
    y = tuple(sub(digits[i], digits[10-i]) for i in range(5))
    gs = affine_gaps(y)
    guards.extend(sub(sub(v, gs[minimum]), ONE if j < minimum else ZERO)
                  for j, v in enumerate(gs) if j != minimum)
    return y, guards


def inventory(cert):
    require(cert['format'] == 'eleven-centered-minimum-trees-v1', 'Certificate version')
    require((cert['horizon'], cert['threshold'], cert['one_step_growth']) ==
            (m.HORIZON, m.THRESHOLD, m.GROWTH), 'Declared universal assertion')
    cases = cert['cases']
    require(tuple(c['index'] for c in cases) == tuple(range(6)), 'All six possible source minima')
    return cases


def expand(case):
    index = case['index']
    require(type(index) is int and 0 <= index < 6, 'Source minimum index')
    gs = affine_gaps(INITIAL)
    minimum = gs[index]
    root = [sub(g, minimum) for g in gs]+[sub(mul(-m.THRESHOLD-1, ONE), minimum)]
    nodes = case['nodes']
    require(isinstance(nodes, list) and nodes and nodes[0]['parent'] is None, 'Nonempty rooted tree')
    states = [(INITIAL, root, 0)]
    incoming = [0]*len(nodes)
    for i, node in enumerate(nodes):
        if i:
            parent = node['parent']
            require(type(parent) is int and 0 <= parent < i, 'Earlier parent')
            x, guards, depth = states[parent]
            y, sorting = successor(x, node['branch'], node['minimum'])
            gy = affine_gaps(y)
            extra = sorting+[sub(minimum, gy[node['minimum']])]
            extra += [sub(g, sub(minimum, mul(m.GROWTH*(depth+1), ONE))) for g in gy]
            states.append((y, guards+extra, depth+1))
        require(states[i][2] < m.HORIZON, 'Every node precedes the declared horizon')
        children = node['children']
        require(isinstance(children, list) and len(children) == len(set(children)), 'Distinct children')
        for child in children:
            require(type(child) is int and i < child < len(nodes) and nodes[child]['parent'] == i,
                    'Correct forward edge and parent')
            incoming[child] += 1
    require(incoming == [0]+[1]*(len(nodes)-1), 'One connected tree with complete edges')
    return minimum, states


def primitive(x):
    digits = [-6]
    for i in range(5):
        a = x[i]-6*int(i == 4)
        b = -x[4-i]-6*int(i < 4)
        digits.extend((z3.If(a >= b, a, b), z3.If(a <= b, a, b)))
    return tuple(digits[i]-digits[10-i] for i in range(5))


def numeric_gaps(x): return (-x[0], *[a-b for a, b in zip(x, x[1:])], x[-1])


def query(case, expanded, ni, omit_child=None, drop_lattice=False):
    ctx = z3.Context()
    x0 = z3.Ints('x0 x1 x2 x3 x4', ctx=ctx)
    def expr(v): return sum(a*b for a, b in zip(v[:5], x0))+v[5]
    minimum, states = expanded
    x, guards, depth = states[ni]
    s = z3.SolverFor('QF_LIA', ctx=ctx)
    s.add(*[expr(g) >= 0 for g in guards])
    if not drop_lattice:
        beta = z3.Int('beta', ctx=ctx)
        z = z3.Ints('z0 z1 z2 z3 z4', ctx=ctx)
        s.add(beta >= 0, beta < 6,
              *[a == 6*b-beta*v for a, b, v in zip(x0, z, m.V)])
    y = primitive(tuple(map(expr, x)))
    gy = numeric_gaps(y)
    s.add(z3.Or(*[g <= expr(minimum) for g in gy]),
          *[g >= expr(minimum)-m.GROWTH*(depth+1) for g in gy])
    for ci in case['nodes'][ni]['children']:
        if ci == omit_child: continue
        child = case['nodes'][ci]
        _, sorting = successor(x, child['branch'], child['minimum'])
        s.add(z3.Not(z3.And(*[expr(g) >= 0 for g in sorting])))
    return s.to_smt2()


def dependency():
    ctx = z3.Context()
    x = z3.Reals('x0 x1 x2 x3 x4', ctx=ctx)
    minimum = z3.Real('minimum', ctx=ctx)
    s = z3.SolverFor('QF_LRA', ctx=ctx)
    s.add(*[g >= minimum for g in numeric_gaps(x)],
          z3.Or(*[g < minimum-m.GROWTH for g in numeric_gaps(primitive(x))]))
    return s.to_smt2()


def solve(smt, logic='QF_LIA'):
    ctx = z3.Context()
    s = z3.SolverFor(logic, ctx=ctx)
    s.set(timeout=60000)
    s.add(z3.parse_smt2_string(smt.replace('(check-sat)', ''), ctx=ctx))
    start = time.monotonic()
    zs = str(s.check())
    zt = time.monotonic()-start
    exe = shutil.which('yices-smt2')
    require(exe is not None, 'External yices-smt2 is required')
    start = time.monotonic()
    r = subprocess.run([exe, '--timeout=2'], input=f'(set-logic {logic})\n'+smt,
                       text=True, capture_output=True, timeout=7)
    ys = r.stdout.strip()
    can_split = logic == 'QF_LIA' and '(declare-fun beta () Int)' in smt
    empty_exit = r.returncode if r.returncode and not (ys or r.stderr.strip()) else None
    if empty_exit is not None and can_split:
        # An empty combined result settles nothing. Require fresh decisions in all six classes.
        ys = 'unknown'
    else:
        require(r.returncode == 0 and ys in ('sat', 'unsat', 'unknown'),
                f'Unexpected Yices response (exit {r.returncode}): {r.stdout}{r.stderr}')
    initial_ys = ys
    split = []
    if ys == 'unknown' and can_split:
        # The source already asserts 0<=beta<6. These six disjoint
        # cases exhaust it exactly; UNKNOWN is never accepted itself.
        for beta in range(6):
            part = smt.replace('(check-sat)', f'(assert (= beta {beta}))\n(check-sat)')
            r = subprocess.run([exe, '--timeout=60'], input=f'(set-logic {logic})\n'+part,
                               text=True, capture_output=True, timeout=65)
            status = r.stdout.strip()
            require(r.returncode == 0 and status in ('sat', 'unsat', 'unknown'), 'Unexpected split Yices response')
            split.append(dict(beta=beta, status=status, sha256=hashlib.sha256(part.encode()).hexdigest()))
        ys = 'unsat' if all(r['status'] == 'unsat' for r in split) else 'sat' if any(r['status'] == 'sat' for r in split) else 'unknown'
    out = dict(z3=zs, yices=ys, z3_seconds=zt, yices_seconds=time.monotonic()-start,
               sha256=hashlib.sha256(smt.encode()).hexdigest(),
               yices_initial=initial_ys, yices_calls=1+len(split))
    if empty_exit is not None:
        out['yices_initial_exit_code'] = empty_exit
    if split: out['yices_split'] = split
    if zs == 'sat':
        out['source'] = [str(s.model().eval(z3.Int(f'x{i}', ctx=ctx), model_completion=True)) for i in range(5)]
    return out


def accepted_record(r):
    if r['z3'] != 'unsat' or r['yices'] != 'unsat': return False
    if r.get('yices_initial') == 'unknown':
        split = r.get('yices_split', [])
        return (r.get('yices_calls') == 7 and tuple(s['beta'] for s in split) == tuple(range(6))
                and all(s['status'] == 'unsat' for s in split))
    return not r.get('yices_split')


def batch(item):
    global _worker_certificate, _worker_case
    index, lo, hi = item
    if _worker_certificate is None:
        raw = CERT.read_bytes()
        _worker_certificate = (json.loads(raw), hashlib.sha256(raw).hexdigest())
    cert, digest = _worker_certificate
    if _worker_case is None or _worker_case[0] != index:
        case = inventory(cert)[index]
        _worker_case = (index, case, expand(case))
    _, case, expanded = _worker_case
    records = []
    for ni in range(lo, hi):
        smt = query(case, expanded, ni)
        r = solve(smt)
        records.append(dict(node=ni, **r))
        if r['z3'] != 'unsat' or r['yices'] != 'unsat':
            (CHECKPOINT/f'failed_{index}_{ni}.smt2').write_text(smt)
            break
    out = dict(index=index, start=lo, stop=hi,
               certificate_sha256=digest, records=records)
    (CHECKPOINT/f'batch_{index}_{lo}_{hi}.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    require(len(records) == hi-lo and all(accepted_record(r) for r in records),
            f'Incomplete independent coverage in case {index}, batch {lo}:{hi}')
    return dict(index=index, start=lo, stop=hi, nodes=len(records),
                z3_seconds=sum(r['z3_seconds'] for r in records),
                yices_seconds=sum(r['yices_seconds'] for r in records))


def accept(cert, records):
    cases = inventory(cert)
    require(tuple(r['index'] for r in records) == tuple(range(6)), 'Complete ordered case records')
    for case, out in zip(cases, records):
        require(tuple(r['node'] for r in out['records']) == tuple(range(len(case['nodes']))), 'Every node audited once')
        require(all(accepted_record(r) for r in out['records']), 'Both engines accept every coverage node')


def assemble(cert, batches):
    digest = hashlib.sha256(CERT.read_bytes()).hexdigest()
    records = []
    for case in inventory(cert):
        rs = []
        for index, lo, hi in batches:
            if index != case['index']: continue
            item = json.loads((CHECKPOINT/f'batch_{index}_{lo}_{hi}.json').read_text())
            require(item['certificate_sha256'] == digest and (item['index'], item['start'], item['stop']) == (index, lo, hi),
                    'Checkpoint belongs to this exact certificate and interval')
            rs.extend(item['records'])
        records.append(dict(index=case['index'], records=rs))
    accept(cert, records)
    return records


def controls(cert):
    # Every child has a recorded physical witness. Removing one root
    # child must expose at least that independent successor cell.
    case = inventory(cert)[0]
    expanded = expand(case)
    child = case['nodes'][0]['children'][0]
    smt = query(case, expanded, 0, omit_child=child)
    hole = solve(smt)
    require(hole['z3'] == hole['yices'] == 'sat', 'A deliberately omitted child must be detected')
    x = tuple(map(int, hole['source']))
    m.base_class(x)
    require(m.residual(x) > m.THRESHOLD and m.residual(m.centered_step(x)) >= m.residual(x),
            'Independent replay of missing-child witness')
    return [dict(kind='omitted_child', case=0, child=child, **hole)]


def main():
    if not __debug__: raise RuntimeError('Optimized Python is not accepted')
    pa = argparse.ArgumentParser()
    pa.add_argument('--resume', action='store_true')
    pa.add_argument('--workers', type=int, default=4)
    args = pa.parse_args()
    CHECKPOINT.mkdir(exist_ok=True)
    cert = json.loads(CERT.read_text())
    cases = inventory(cert)
    for c in cases: expand(c)
    dep = solve(dependency(), 'QF_LRA')
    require(dep['z3'] == dep['yices'] == 'unsat', 'One-step growth dependency')
    print('dependency', dep, flush=True)
    batches = [(c['index'], lo, min(lo+500, len(c['nodes']))) for c in cases for lo in range(0, len(c['nodes']), 500)]
    digest = hashlib.sha256(CERT.read_bytes()).hexdigest()
    todo = []
    for item in batches:
        index, lo, hi = item
        path = CHECKPOINT/f'batch_{index}_{lo}_{hi}.json'
        if args.resume and path.exists():
            r = json.loads(path.read_text())
            if (r['certificate_sha256'] == digest and tuple(v['node'] for v in r['records']) == tuple(range(lo, hi))
                    and all(accepted_record(v) for v in r['records'])):
                continue
        todo.append(item)
    print(dict(total_nodes=sum(len(c['nodes']) for c in cases), batches=len(batches), pending_batches=len(todo)), flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(batch, item) for item in todo]
        for future in as_completed(futures): print('batch', future.result(), flush=True)
    records = assemble(cert, batches)
    negative = controls(cert)
    rejected = 0
    for kind in ('missing_case', 'missing_node', 'sat', 'unknown'):
        bad = json.loads(json.dumps(records))
        if kind == 'missing_case': bad.pop()
        elif kind == 'missing_node': bad[0]['records'].pop()
        else: bad[0]['records'][0]['yices'] = kind
        try: accept(cert, bad)
        except ArithmeticError: rejected += 1
        else: raise ArithmeticError('Corrupt acceptance record passed')
    # A split with a missing base class cannot certify its original node.
    for broken in ([], [dict(beta=i, status='unsat') for i in range(5)]):
        require(not accepted_record(dict(z3='unsat', yices='unsat', yices_initial='unknown',
                                         yices_calls=7, yices_split=broken)), 'Incomplete coset split rejected')
        rejected += 1
    out = dict(status='complete: all physical source cases and nodes accepted',
               certificate_sha256=digest, dependency=dep, cases=records, controls=negative,
               rejected_corruptions=rejected,
               solver_calls=2+sum(1+r.get('yices_calls', 1) for c in records for r in c['records'])+
                            sum(1+r.get('yices_calls', 1) for r in negative),
               yices_coset_splits=sum(bool(r.get('yices_split')) for c in records for r in c['records']),
               z3_version=z3.get_version_string(),
               yices_version=subprocess.run(['yices-smt2', '--version'], text=True, capture_output=True).stdout.strip())
    (ROOT/'eleven_centered_clock_audit_results.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print({k: v for k, v in out.items() if k not in ('cases', 'dependency', 'controls')}, flush=True)


if __name__ == '__main__': main()
