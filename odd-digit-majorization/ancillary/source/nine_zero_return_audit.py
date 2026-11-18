"""Unbounded Z3/CVC5 audit of sparse nine-slot zero-return sources.

Also certifies G(F(q)) <= G(q)+8 for every positive source, including
zero outputs. No multi-step dispersion descent is asserted here.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
import hashlib
import json
import time
import z3
import cvc5
from odd_gap import subtract_pair
from nine_clock_reduction import dispersion, return_cover, require

ROOT = Path(__file__).resolve().parent
NETWORK = ((0,1),(2,3),(4,5),(6,7),(0,2),(1,3),(4,6),(5,7),
           (1,2),(5,6),(0,4),(3,7),(1,5),(2,6),(1,4),(3,6),
           (2,4),(3,5),(3,4))
KEYS = ('dispersion_increase',)+tuple(f'partial_sum_{r}' for r in range(1,5))+tuple(
    f'return_{beta}_{rank}_{depth}'
    for beta in (0,1) for rank in (2,3) for depth in range(1,5))


def absolute(x):
    return z3.If(x >= 0, x, -x)


def symbolic_cover(B, q, rank):
    u,v,w,_ = q
    if rank == 2:
        return z3.Or(absolute(u-v)<=3, absolute(u+v-B)<=3,
                     absolute(2*u-B)<=4, absolute(2*v-B)<=4, B-u<=3)
    return z3.Or(
        z3.And(v==w, z3.Or(absolute(u-v)<=2, absolute(u+v-B)<=2,
                          absolute(2*v-B)<=3)),
        z3.And(u==v, z3.Or(absolute(u+w-B)<=2, absolute(2*u-B)<=3)),
        z3.And(absolute(2*u-B)<=3, absolute(2*v-B)<=3),
        z3.And(B-u<=3, z3.Or(absolute(2*v-B)<=3, v==w,
                            absolute(v+w-B)<=3)),
        z3.And(absolute(4*u-3*B)<=8,
               z3.Or(*[absolute(4*x-j*B)<=12
                       for x in (v,w) for j in (1,2,3)])))


def symbolic_dispersion(B, q):
    u,v,w,x = q
    gs = (B-u,u-v,v-w,w-x,x)
    return sum(absolute(a-b) for i,a in enumerate(gs) for b in gs[i+1:])


def symbolic_top_sum(B,q,r):
    u,v,w,x = q
    gs = (B-u,u-v,v-w,w-x,x)
    sums = [sum(gs[i] for i in ids) for ids in combinations(range(5),r)]
    value = sums[0]
    for x in sums[1:]:
        value = z3.If(value>=x,value,x)
    return value


def prefix(depth, beta=None, rank=None):
    B = z3.Int('B') if beta is None else 2*z3.Int('a')+beta
    q = z3.Ints('q0_0 q0_1 q0_2 q0_3')
    D = [B>=2, B>q[0], q[0]>=1, q[0]>=q[1], q[1]>=q[2],
         q[2]>=q[3], q[3]>=0]
    if rank is None:
        D.append(q[3]>=1)
    else:
        D.extend((z3.Int('a')>=20, q[rank-1]>=1))
        D.extend(x==0 for x in q[rank:])
    qs = [q]
    for t in range(depth):
        # Exact borrowed digits, with the central maximal digit B-1 removed.
        hi = [z3.If(q[i+1]>0,q[i],z3.If(q[i]>0,q[i]-1,B-1))
              for i in range(3)]+[z3.If(q[3]>0,q[3]-1,B-1)]
        ds = hi+[B-q[3]-1,B-q[2]-1,B-q[1]-1,B-q[0]]
        for k,(i,j) in enumerate(NETWORK):
            a,b = ds[i],ds[j]
            h,l = z3.Ints(f'h{t}_{k} l{t}_{k}')
            D.extend((h==z3.If(a>=b,a,b),l==a+b-h))
            ds[i],ds[j] = h,l
        q = z3.Ints(f'q{t+1}_0 q{t+1}_1 q{t+1}_2 q{t+1}_3')
        D.extend((q[0]==B-1-ds[7], q[1]==ds[0]-ds[6],
                  q[2]==ds[1]-ds[5], q[3]==ds[2]-ds[4],
                  B>q[0], q[0]>=1, q[0]>=q[1], q[1]>=q[2],
                  q[2]>=q[3], q[3]>=0))
        qs.append(q)
    return B,qs,D


def obligation(key):
    if key == 'dispersion_increase':
        B,qs,D = prefix(1)
        D.append(symbolic_dispersion(B,qs[1])>symbolic_dispersion(B,qs[0])+8)
    elif key.startswith('partial_sum_'):
        r = int(key.rsplit('_',1)[1])
        B,qs,D = prefix(1)
        D.append(symbolic_top_sum(B,qs[1],r)>symbolic_top_sum(B,qs[0],r)+1)
    else:
        beta,rank,depth = map(int,key.split('_')[1:])
        B,qs,D = prefix(depth,beta,rank)
        D.extend(q[3]>=1 for q in qs[1:-1])
        D.extend((qs[-1][3]==0,z3.Not(symbolic_cover(B,qs[0],rank))))
    s = z3.SolverFor('QF_LIA')
    s.add(*D)
    return B,qs,s.to_smt2()


def cvc_status(smt):
    s = cvc5.Solver()
    s.setLogic('QF_LIA')
    s.setOption('tlimit-per','90000')
    parser = cvc5.InputParser(s)
    parser.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6,
                         smt.replace('(check-sat)',''),'nine-zero-return')
    sm = parser.getSymbolManager()
    while True:
        c = parser.nextCommand()
        if c.isNull():
            break
        c.invoke(s,sm)
    t = time.monotonic()
    return str(s.checkSat()),time.monotonic()-t


def solve(item):
    key,smt = item
    s = z3.SolverFor('QF_LIA')
    s.set(timeout=90000)
    s.add(z3.parse_smt2_string(smt.replace('(check-sat)','')))
    t = time.monotonic()
    zs = str(s.check())
    zt = time.monotonic()-t
    cs,ct = cvc_status(smt)
    return dict(key=key,z3=zs,cvc5=cs,z3_seconds=zt,cvc5_seconds=ct,
                sha256=hashlib.sha256(smt.encode()).hexdigest())


def accept(records):
    require(tuple(r['key'] for r in records)==KEYS,'Complete ordered inventory')
    require(all(r['z3']==r['cvc5']=='unsat' for r in records),
            'Every obligation must be UNSAT in both solvers')


def controls():
    out = []
    for key in ('false_increase_four','false_rank_two_empty_cover'):
        if key == 'false_increase_four':
            B,qs,D = prefix(1)
            D.append(symbolic_dispersion(B,qs[1])>symbolic_dispersion(B,qs[0])+4)
        else:
            B,qs,D = prefix(1,1,2)
            D.append(qs[1][3]==0)
        s = z3.SolverFor('QF_LIA')
        s.set(timeout=90000)
        s.add(*D)
        smt = s.to_smt2()
        zs = str(s.check())
        cs,_ = cvc_status(smt)
        require(zs==cs=='sat','False strengthening must be SAT in both solvers')
        m = s.model()
        base = m.eval(B).as_long()
        path = [tuple(m.eval(x).as_long() for x in q) for q in qs]
        require(subtract_pair(base,path[0])==path[1],'Independent witness replay')
        if key == 'false_increase_four':
            require(dispersion(base,path[1])>dispersion(base,path[0])+4,
                    'Numeric false potential bound')
        else:
            require(path[0][2:]==(0,0) and path[1][3]==0
                    and return_cover(base,path[0]),'Numeric zero return control')
        out.append(dict(key=key,z3=zs,cvc5=cs,base=base,path=path,smtlib=smt))
    return out


def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is not accepted.')
    statements = [dict(key=k,negated_lemma_smtlib=obligation(k)[2]) for k in KEYS]
    (ROOT/'nine_zero_return_obligations.json').write_text(
        json.dumps(statements,indent=2)+'\n')
    records = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        fs = [pool.submit(solve,(r['key'],r['negated_lemma_smtlib'])) for r in statements]
        for f in as_completed(fs):
            r = f.result()
            records.append(r)
            print(r,flush=True)
    records.sort(key=lambda r:KEYS.index(r['key']))
    accept(records)
    negative = controls()
    rejected = 0
    for kind in ('missing','sat','unknown'):
        bad = [dict(r) for r in records]
        if kind == 'missing':
            bad.pop()
        else:
            bad[0]['cvc5'] = kind
        try:
            accept(bad)
        except ArithmeticError:
            rejected += 1
        else:
            raise ArithmeticError('Corrupted result accepted')
    data = dict(records=records,negative_controls=negative,
                corruptions_rejected=rejected,actual_solver_calls=46,
                z3_version=z3.get_version_string(),cvc5_version=cvc5.__version__,
                scope='Sparse zero-return cover for B>=40; every ordered gap partial sum increases by at most one, and dispersion by at most eight. No multi-step descent or global nine-slot capture-time bound.')
    (ROOT/'nine_zero_return_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: twenty-one universal lemmas in Z3/CVC5; two SAT controls; three rejected corruptions.',flush=True)


if __name__ == '__main__':
    main()
