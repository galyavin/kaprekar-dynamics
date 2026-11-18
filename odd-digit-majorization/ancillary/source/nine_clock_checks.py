"""Independent arithmetic replay for the global nine-slot clock theorem."""
from pathlib import Path
from random import Random
from copy import deepcopy
import hashlib,json
import sympy as sp
from odd_gap import subtract_pair,gaps
from nine_clock import dispersion,terminal_residual,pair_clock_bound,number_clock_bound
from nine_clock_audit import inventory,expand_case,accept,require

ROOT=Path(__file__).resolve().parent


def evaluate(row,point):
    return sum(x*y for x,y in zip(row[:5],point))+row[5]


def verify_witnesses(cases):
    count=0
    prefixes={}
    updates=0
    for case in cases:
        _,states=expand_case(case)
        order=tuple(map(int,case['order']))
        for node,(q,guards,depth) in zip(case['nodes'],states):
            if 'source_witness' not in node:
                continue
            point=node['source_witness']
            require(len(point)==5 and all(type(x) is int for x in point),'Integer witness')
            B=point[0]
            source=tuple(point[1:])
            require(B>=40 and B>source[0]>=source[1]>=source[2]>=source[3]>=1,
                    'Valid original difference source')
            require(all(evaluate(row,point)>=0 for row in guards),'Every affine path guard')
            g=gaps(B,source)
            require(all(g[i]<=g[j] for i,j in zip(order,order[1:])),'Source gap order')
            key=tuple(point)
            if key not in prefixes:
                prefixes[key]=[source]
            path=prefixes[key]
            while len(path)<=depth:
                path.append(subtract_pair(B,path[-1]))
                updates+=1
            require(tuple(evaluate(row,point) for row in q)==path[depth],
                    'Affine certificate state versus independent digit columns')
            require(all(v[-1]>=1 and dispersion(B,v)>=dispersion(B,source)
                        for v in path[:depth+1]),'Nondecreasing positive prefix')
            require(terminal_residual(B,source)>100,'Outside exceptional tube')
            count+=1
    return dict(branch_witnesses=count,distinct_sources=len(prefixes),raw_updates=updates)


def eight_step_family():
    t=sp.Symbol('t')
    B=34*t-9
    path=[(26*t-12,20*t-9,14*t-6,8*t-3),
          (26*t-6,18*t-10,12*t-6,6*t-3),
          (28*t-6,20*t-4,14*t,6*t-4),
          (28*t-5,22*t-3,14*t,6*t-4),
          (28*t-5,22*t-1,16*t+1,8*t-3),
          (28*t-6,20*t-1,14*t+2,6*t-2),
          (28*t-7,22*t-3,14*t+1,6*t-3),
          (28*t-6,22*t-5,16*t,8*t-4),
          (28*t-7,20*t-1,14*t-1,6*t-5),
          (28*t-4,22*t-3,14*t+2,6*t)]
    path=[tuple(map(sp.expand,q)) for q in path]
    guards=0
    identities=0
    values=[]
    def nonnegative(v):
        nonlocal guards
        v=sp.expand(v)
        require(sp.Poly(v,t).degree()<=1 and v.coeff(t)>=0 and v.subs(t,6)>=0,
                'Exact affine inequality for every t>=6')
        guards+=1
    def ordered(values):
        values=sorted(values,key=lambda x:(x.coeff(t),x.subs(t,0)),reverse=True)
        for a,b in zip(values,values[1:]):
            nonnegative(a-b)
        return values
    for i,q in enumerate(path):
        for v in (B-1-q[0],q[0]-q[1],q[1]-q[2],q[2]-q[3],q[3]-1):
            nonnegative(v)
        gs=ordered([sp.expand(x) for x in gaps(B,q)])
        values.append(sp.expand(sum((4-2*j)*x for j,x in enumerate(gs))))
        if i==9:
            break
        u,v,w,x=q
        ds=[u,v,w,x-1,B-1,B-x-1,B-w-1,B-v-1,B-u]
        descending=sum(q[j]*B**(8-j) for j in range(4))
        ascending=sum(q[j]*B**j for j in range(4))
        require(sp.expand(descending-ascending-sum(d*B**(8-j) for j,d in enumerate(ds)))==0,
                'Whole-number polynomial subtraction identity')
        identities+=1
        ds=ordered([sp.expand(x) for x in ds])
        require(all(sp.expand(ds[j]-ds[8-j]-path[i+1][j])==0 for j in range(4)),
                'Unbounded actual difference update')
        identities+=4
    expected=[12*t+c for c in (24,32,36,34,34,24,30,34,28,20)]
    require(values==expected,'Exact dispersion path')
    for v in values[:9]:
        nonnegative(v-values[0])
    nonnegative(values[0]-values[9]-1)
    m=6*t-3
    M=8*t+3
    nonnegative(2*m-M)
    require(sp.expand(B-4*m-M)==2*t,'Unbounded terminal residual')
    raw=0
    for v in (6,7,10,51,100,1000,10**6,10**12):
        base=int(B.subs(t,v))
        points=[tuple(int(x.subs(t,v)) for x in q) for q in path]
        for a,b in zip(points,points[1:]):
            require(subtract_pair(base,a)==b,'Large-base family replay')
            raw+=1
    return dict(affine_guards=guards,polynomial_identities=identities,raw_updates=raw,
                base=str(B),parameter_range='integer t>=6',residual='2*t',
                path=[[str(v) for v in q] for q in path],
                dispersion=[str(v) for v in values])


def integer_digits(B,n):
    ds=[]
    for _ in range(9):
        n,d=divmod(n,B)
        ds.append(d)
    require(n==0,'Nine-slot integer range')
    return ds


def number_step(B,n):
    ds=sorted(integer_digits(B,n))
    a=b=0
    for x,y in zip(ds,reversed(ds)):
        a=B*a+x
        b=B*b+y
    return b-a


def number_differences(B,n):
    ds=sorted(integer_digits(B,n),reverse=True)
    return tuple(ds[i]-ds[-1-i] for i in range(4))


def complete_captures():
    rng=Random(560211)
    records=[]
    number_updates=pair_updates=0
    for B in (195,196,197,198,199,1000,1001,1002,1003,1004):
        for _ in range(12):
            ds=[rng.randrange(B) for i in range(9)]
            n=0
            for d in ds:
                n=B*n+d
            source=number_differences(B,n)
            require(source[0]>0,'Nonconstant full-number input')
            seen={}
            cur=n
            while cur not in seen:
                seen[cur]=len(seen)
                nxt=number_step(B,cur)
                require(number_differences(B,nxt)==subtract_pair(B,number_differences(B,cur)),
                        'Whole-number and independent column transitions agree')
                cur=nxt
                number_updates+=1
                require(len(seen)<100000,'Finite-check safety cap')
            tn=seen[cur]
            pairs={}
            q=source
            while q not in pairs:
                pairs[q]=len(pairs)
                q=subtract_pair(B,q)
                pair_updates+=1
                require(len(pairs)<100000,'Pair safety cap')
            tq=pairs[q]
            require(tq<=tn<=tq+1,'Actual first-capture bridge')
            require(tq<=pair_clock_bound(B) and tn<=number_clock_bound(B),'Global bound')
            records.append(dict(base=B,source=source,pair_capture=tq,number_capture=tn))
    return dict(captures=len(records),number_updates=number_updates,pair_updates=pair_updates,
                records=records)


def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is not accepted.')
    raw=(ROOT/'nine_clock_certificates.json').read_bytes()
    certificate=json.loads(raw)
    audit=json.loads((ROOT/'nine_clock_audit_results.json').read_text())
    cases=inventory(certificate)
    require(audit['accepted'] and audit['certificate_sha256']==hashlib.sha256(raw).hexdigest(),
            'Accepted audit matches this certificate')
    accept(audit['records'],cases)
    require(audit['dependency']['z3']==audit['dependency']['yices']=='unsat','Hint dependency')
    rejected=0
    for kind in ('missing_case','wrong_horizon','duplicate_sort_label','self_parent','missing_edge'):
        try:
            if kind in ('missing_case','wrong_horizon'):
                bad=dict(certificate)
                if kind=='missing_case':
                    bad['cases']=cases[:-1]
                else:
                    bad['horizon']=8
                inventory(bad)
            else:
                bad=deepcopy(next(c for c in cases if len(c['nodes'])>1))
                if kind=='duplicate_sort_label':
                    bad['nodes'][1]['digits_order'][1]=bad['nodes'][1]['digits_order'][0]
                elif kind=='self_parent':
                    bad['nodes'][1]['parent']=1
                else:
                    bad['nodes'][0]['children'].pop()
                expand_case(bad)
        except ArithmeticError:
            rejected+=1
        else:
            raise ArithmeticError('Corrupted certificate accepted')
    witnesses=verify_witnesses(cases)
    print(witnesses,flush=True)
    family=eight_step_family()
    print({k:v for k,v in family.items() if k not in ('path','dispersion')},flush=True)
    captures=complete_captures()
    print({k:v for k,v in captures.items() if k!='records'},flush=True)
    # Reuse the unchanged complete finite graphs and zero-face audit from chapter 55.
    earlier=json.loads((ROOT/'nine_clock_reduction_results.json').read_text())
    graphs=[r for r in earlier['holdouts'] if r['mode']=='graph']
    require(tuple(r['base'] for r in graphs)==(128,129)
            and all(r['latest_descent']<=9 and r['maximum_one_step_increase']<=8 for r in graphs),
            'Retained independent graph controls')
    result=dict(witnesses=witnesses,eight_step_family=family,captures=captures,
                certificate_corruptions_rejected=rejected,
                reused_complete_graph_states=sum(r['states'] for r in graphs),
                scope='Independent arithmetic and symbolic controls for the accepted unbounded tree audit; finite experiments do not replace it.')
    (ROOT/'nine_clock_check_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: certificate branch witnesses, infinite eight-step obstruction and complete original-number captures.',flush=True)


if __name__=='__main__':
    main()
