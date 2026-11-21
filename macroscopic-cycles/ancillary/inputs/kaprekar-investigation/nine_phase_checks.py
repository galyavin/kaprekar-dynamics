"""Solver-free affine proof replay and independent nine-slot first arrivals.

Run --generate once to regenerate the saved sorting itineraries. Normal runs
reconstruct every guard and transition, prove the parameter-box bounds with
rational arithmetic, and trace complete numeric subtraction independently.
"""
import argparse
import copy
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
from math import factorial

from nine_phase import (MIN_BASE,BOX,TERMINAL,BASE_ROW,ONE,SEED,STAGES,
                        evaluate,stage,source,exact_arrival,source_mass)


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def shift(matrix, delta):
    return tuple(tuple(x+(d if j == len(r)-1 else 0) for j,x in enumerate(r))
                 for r,d in zip(matrix,delta))


def set_count(matrix, k, value):
    return tuple(tuple(r[-1]+value*r[k] if j == len(r)-1 else 0 if j == k else x
                       for j,x in enumerate(r)) for r in matrix)


def compose(row,matrix):
    return tuple(sum(c*r[j] for c,r in zip(row,matrix)) for j in range(len(matrix[0])))


def step(matrix, base, one, order):
    """Reconstruct from nine borrowed digits; no affine-word generator used."""
    require(sorted(order) == list(range(9)), 'sorting permutation')
    ds = (*matrix[:-1],sub(matrix[-1],one),sub(base,one),
          *[sub(sub(base,r),one) for r in reversed(matrix[1:])],sub(base,matrix[0]))
    guards = [sub(matrix[-1],one),sub(sub(base,one),matrix[0]),
              *[sub(a,b) for a,b in zip(matrix,matrix[1:])]]
    guards += [sub(sub(ds[i],ds[j]),one) if i > j else sub(ds[i],ds[j])
               for i,j in zip(order,order[1:])]
    result = tuple(sub(ds[order[i]],ds[order[-1-i]]) for i in range(4))
    return result,guards


def replay(record):
    matrix = tuple(tuple(int(i == j) for j in range(6)) for i in range(4))
    guards, offsets = [],[]
    for order in record['orders']:
        offsets.append(matrix)
        matrix,extra = step(matrix,BASE_ROW,ONE,order)
        guards.extend(extra)
    require(matrix == tuple(map(tuple,record['matrix'])), 'record matrix')
    return matrix,guards,offsets


def box_bound(row):
    slope = row[4]+sum(c*(lo if c >= 0 else hi) for c,(lo,hi) in zip(row,BOX))
    require(slope >= 0 and MIN_BASE*slope+row[5] >= 0, ('box guard',row,slope))


def edge_bound(row):
    # Variables (V,B,1), 4B/5+256 <= V <= 403B/500.
    v,b,c = row
    slope = b+v*(Q(4,5) if v >= 0 else Q(403,500))
    constant = c+(256*v if v >= 0 else 0)
    require(slope >= 0 and MIN_BASE*slope+constant >= 0, ('edge guard',row))


def specs():
    for slow in (False,True):
        for k,length,delta in STAGES:
            start = stage(k,slow)
            yield f'{"slow" if slow else "fast"}:{k}',start,length,shift(start,delta)
    start = set_count(stage(1,True),1,1)
    target = ((1,0,0,0,0,-1),(2,0,0,0,-1,-2),
              (3,0,0,0,-2,-3),(-1,0,0,0,1,-1))
    yield 'connector',start,40,target
    yield 'fixed',TERMINAL,1,TERMINAL


def terminal_generate(offset):
    # Rows (coefficient of b, constant), B=5b, b>=200000.
    matrix = ((4,offset),(3,2*offset),(2,3*offset),(1,-offset-2))
    seen,orders = {},[]
    while matrix not in seen:
        seen[matrix] = len(orders)
        base,one = (5,0),(0,1)
        ds = (*matrix[:-1],sub(matrix[-1],one),sub(base,one),
              *[sub(sub(base,r),one) for r in reversed(matrix[1:])],sub(base,matrix[0]))
        order = sorted(range(9),key=lambda i:(-evaluate(ds[i],(200000,1)),i))
        orders.append(order)
        matrix,_ = step(matrix,base,one,order)
        require(len(orders) <= 10000,'terminal generation cap')
    return {'offset':offset,'arrival':seen[matrix],
            'period':len(orders)-seen[matrix],'orders':orders}


def generate(path):
    from odd_affine import word
    records = []
    for key,start,length,target in specs():
        q = tuple(evaluate(r,SEED) for r in start)
        certificate = word(MIN_BASE,q,length)
        records.append({'key':key,'orders':certificate['orders'],'matrix':certificate['matrix']})
    certificate = word(MIN_BASE,(803611,607222,410833,196387),32)
    records.append({'key':'edge','orders':certificate['orders'],'matrix':certificate['matrix']})
    data = {'min_base':MIN_BASE,'box':[[str(x) for x in pair] for pair in BOX],
            'records':records,'terminals':[terminal_generate(e) for e in (256,257)]}
    path.write_text(json.dumps(data,separators=(',',':'))+'\n')


def verify(data):
    require(data['min_base'] == MIN_BASE,'base domain')
    require(data['box'] == [[str(x) for x in pair] for pair in BOX],'box domain')
    expected = list(specs())
    require([r['key'] for r in data['records']] == [s[0] for s in expected]+['edge'], 'record inventory')
    guard_count,capture_count = 0,0
    for record,(key,start,length,target) in zip(data['records'],expected):
        require(len(record['orders']) == length,('word length',key))
        matrix,guards,offsets = replay(record)
        require(tuple(compose(r,(*start,BASE_ROW,ONE)) for r in matrix) == target,('return identity',key))
        anchors = [start]
        if ':' in key:
            k = int(key[-1])
            anchors.append(set_count(start,k,2 if key == 'slow:1' else 1))
        for anchor in anchors:
            for guard in guards:
                box_bound(compose(guard,(*anchor,BASE_ROW,ONE)))
                guard_count += 1
            if key.startswith('slow') or key == 'connector':
                # The eventual 21-cycle has first coordinate <= 4B/5+3.
                for offset in offsets:
                    row = list(compose(offset[0],(*anchor,BASE_ROW,ONE)))
                    row[4] -= Q(4,5);row[5] -= 4
                    box_bound(row);capture_count += 1
        if key == 'fast:1':
            last = set_count(start,1,1)
            hits = []
            for j,offset in enumerate([*offsets,matrix]):
                rows = tuple(compose(r,(*last,BASE_ROW,ONE)) for r in offset)
                differences = tuple(sub(r,t) for r,t in zip(rows,TERMINAL))
                require(all(not any(r[:-1]) for r in differences),'constant final-word differences')
                if all(r[-1] == 0 for r in differences):
                    hits.append(j)
            require(hits == [42],'first fixed capture')
    record = data['records'][-1]
    require(len(record['orders']) == 32,'edge length')
    matrix,guards,offsets = replay(record)
    edge = ((1,0,0),(2,-1,0),(3,-2,0),(-1,1,-2))
    base,one = (0,1,0),(0,0,1)
    require(tuple(compose(r,(*edge,base,one)) for r in matrix)
            == shift(edge,(-2,-4,-6,2)),'edge displacement')
    for guard in guards:
        edge_bound(compose(guard,(*edge,base,one)));guard_count += 1
    for offset in offsets:
        r = compose(offset[0],(*edge,base,one))
        require(r[:2] == (1,0) and -2 <= r[2] <= 2,'edge safe capture')
        capture_count += 1
    require([r['offset'] for r in data['terminals']] == [256,257], 'terminal inventory')
    terminals = []
    terminal_guards = 0
    for record in data['terminals']:
        e = record['offset']
        matrix = ((4,e),(3,2*e),(2,3*e),(1,-e-2))
        path,seen = [],{}
        base,one = (5,0),(0,1)
        for order in record['orders']:
            require(matrix not in seen,'premature symbolic repetition')
            seen[matrix] = len(path);path.append(matrix)
            matrix,guards = step(matrix,base,one,order)
            for a,c in guards:
                require(a >= 0 and 200000*a+c >= 0,'terminal affine guard')
                terminal_guards += 1
        require(matrix in seen,'terminal closure')
        h,period = seen[matrix],len(path)-seen[matrix]
        require((h,period) == (record['arrival'],record['period']) == (16*e+(18 if e%2 == 0 else 34),21),'terminal time')
        cycle = path[h:]
        require(all(r[0][0] == 4 and -1 <= r[0][1] <= 3 for r in cycle),'cycle extent')
        # Different affine states must stay different at EVERY b >= 200000.
        # A repeat of numeric values could precede the symbolic repeat.
        # They all have the same leading tuple after the first fixed number
        # of steps here: certify it directly, so distinct constants suffice.
        require(all(tuple(r[0] for r in state) == (4,3,2,1) for state in path),
                'terminal common leading tuple')
        terminals.append({'offset':e,'pair_time':h,'period':period})
    return {'return_guard_checks':guard_count,'safe_capture_checks':capture_count,
            'terminal_guard_checks':terminal_guards,'terminal_families':terminals}


def raw_step(base,q):
    digits = sorted((*q,0,0,0,0,0))
    low=high=0
    for d in digits:low=base*low+d
    for d in reversed(digits):high=base*high+d
    value,out=high-low,[]
    for _ in range(9):
        value,d=divmod(value,base);out.append(d)
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(4))


def raw_arrival(base,q):
    seen={}
    while q not in seen:
        seen[q]=len(seen);q=raw_step(base,q)
    return seen[q],len(seen)-seen[q]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true')
    args=parser.parse_args();root=Path(__file__).resolve().parent
    path=root/'nine_phase_certificates.json'
    if args.generate:generate(path)
    data=json.loads(path.read_text());result=verify(data)
    controls=0
    for mutation in ('delete','matrix','domain'):
        broken=copy.deepcopy(data)
        if mutation == 'delete':broken['records'][0]['orders'].pop()
        elif mutation == 'matrix':broken['records'][0]['matrix'][0][-1]+=1
        else:broken['box'][0][0]='0'
        try:verify(broken)
        except AssertionError:controls+=1
        else:raise AssertionError('Corrupt certificate accepted.')
    # Rational elimination independently checks rank and lattice separation.
    A=[list(map(Q,r[:4])) for r in stage(3)]
    augmented=[r+[Q(-4 if i == 2 else 0)] for i,r in enumerate(A)]
    determinant=Q(1)
    for i in range(4):
        j=next(j for j in range(i,4) if augmented[j][i])
        if j != i:augmented[i],augmented[j]=augmented[j],augmented[i];determinant*=-1
        pivot=augmented[i][i];determinant*=pivot
        augmented[i]=[x/pivot for x in augmented[i]]
        for j in range(4):
            if j != i:
                c=augmented[j][i];augmented[j]=[x-c*y for x,y in zip(augmented[j],augmented[i])]
    translation=[row[-1] for row in augmented]
    require(abs(determinant) == 1440 and any(x.denominator != 1 for x in translation),'distinct lattice cosets')
    # Positivity of every leading gap on all box vertices.
    for params in product(*BOX):
        q=[evaluate(r,(*params,1,0)) for r in stage(3)]
        require(1>q[0]>q[1]>q[2]>q[3]>0,'strict density support')
    mass=source_mass();require(mass>0,'positive mass')
    # An independent tensor Newton-Cotes rule is exact for degree <=5 in
    # each parameter; it does not expand or integrate polynomial coefficients.
    weights=(19,75,50,50,75,19);quadrature=Q(0)
    for indices in product(range(6),repeat=4):
        params=tuple(lo+(hi-lo)*Q(i,5) for i,(lo,hi) in zip(indices,BOX))
        q=tuple(evaluate(r,(*params,1,0)) for r in stage(3))
        value=Q(factorial(9))
        for gap in (1-q[0],*[a-b for a,b in zip(q,q[1:])],q[-1]):value*=gap
        for i,(lo,hi) in zip(indices,BOX):value*=Q(weights[i],288)*(hi-lo)
        quadrature+=value
    require(quadrature==mass,'independent exact source integration')
    cost_rows=((3,2,1,1),(6,12,6,4),(10,20,30,20),(2,4,6,-2))
    lambda_coefficients=(Q(4,15),Q(19,60),Q(1,6),Q(7,5))
    all_coefficients=(Q(52,45),Q(-23,180),Q(67,18),Q(-59,15))
    require(tuple(evaluate(r,lambda_coefficients) for r in cost_rows)==(3,12,42,0),
            'source cost functional and edge failure')
    require(tuple(evaluate(r,all_coefficients) for r in cost_rows)==(3,12,42,32),
            'compatible four-return cost functional')
    require(evaluate((1,2,3,-1),all_coefficients)==16,'varying terminal value')
    rng=Random(942032);checks=[];steps=0
    for B in (1000000,1000005,2000000,2000005):
        params=tuple(rng.randrange(int(lo*B)+1,int(hi*B)) for lo,hi in BOX)
        pair=[]
        for slow in (False,True):
            expected=exact_arrival(B,params,slow)
            actual,period=raw_arrival(B,expected['source'])
            require((actual,period)==(expected['pair_time'],expected['period']),('direct total clock',B,params,slow,actual,expected))
            steps+=actual+period;pair.append(expected)
        checks.append({'base':B,'parameters':params,'fast':pair[0],'slow':pair[1]})
    result.update({'negative_controls':controls,'determinant':str(determinant),
                   'parameter_translation':[str(x) for x in translation],
                   'source_mass_each':str(mass),'source_mass_each_decimal':float(mass),
                   'certificate_sha256':sha256(path.read_bytes()).hexdigest(),
                   'independent_numeric_steps':steps,'independent_full_clocks':checks})
    (root/'nine_phase_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f"PASS: {result['return_guard_checks']} return guards, "
            f"{result['terminal_guard_checks']} terminal guards, "
            f"{result['safe_capture_checks']} safe-capture checks; {controls} corruption controls.\n"
            f"PASS: 8 independent full clocks, {steps} numeric subtraction steps.\n"
            f"Each phase family has exact limiting input mass {mass} = {float(mass):.12g}.\n"
            "All affine guards hold for unbounded bases in the stated domain; no solver is used.\n")
    (root/'nine_phase_verification.txt').write_text(report);print(report,end='')


if __name__ == '__main__':main()
