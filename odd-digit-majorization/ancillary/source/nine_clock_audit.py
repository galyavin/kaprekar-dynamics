"""Verify the complete nine-step dispersion certificate, over unbounded bases.

Each tree node has a single actual update and linear ancestral guards.
Children are disjoint digit/gap sorting cells. Regenerated coverage
queries must be UNSAT in both Z3 and Yices. No saved SMT is trusted.
"""
from concurrent.futures import ProcessPoolExecutor,as_completed
from itertools import permutations
from pathlib import Path
import hashlib,json,subprocess,time
import z3
from odd_gap import subtract_pair
from nine_clock import HORIZON,THRESHOLD,BASE_MINIMUM,dispersion,terminal_residual
from nine_zero_return_audit import NETWORK,obligation as dependency_obligation,prefix as full_prefix
from seven_clock_audit import yices_status

ROOT=Path(__file__).resolve().parent
KEYS=tuple((''.join(map(str,p)),case) for p in permutations(range(5)) for case in (1,2))
VARS=z3.Ints('B u v w x')
ZERO=(0,0,0,0,0,0)
ONE=(0,0,0,0,0,1)
BASE=(1,0,0,0,0,0)
INITIAL=tuple(tuple(int(i==j) for j in range(6)) for i in range(1,5))


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def add(*xs):
    return tuple(map(sum,zip(*xs)))


def mul(c,x):
    return tuple(c*a for a in x)


def sub(a,b):
    return add(a,mul(-1,b))


def expr(a):
    return sum(x*y for x,y in zip(a[:5],VARS))+a[5]


def gaps(q):
    return (sub(BASE,q[0]),*(sub(a,b) for a,b in zip(q,q[1:])),q[-1])


def sorting_guards(values,order):
    require(isinstance(order,list) and all(type(i) is int for i in order)
            and sorted(order)==list(range(len(values))),'Complete sorting permutation')
    # Descending values, with increasing original index at equality.
    return [sub(sub(values[i],values[j]),ONE if i>j else ZERO)
            for i,j in zip(order,order[1:])]


def affine_child(q,od,og,g0,depth):
    d=[q[0],q[1],q[2],sub(q[3],ONE),sub(sub(BASE,q[3]),ONE),
       sub(sub(BASE,q[2]),ONE),sub(sub(BASE,q[1]),ONE),sub(BASE,q[0])]
    dg=sorting_guards(d,od)
    d=[d[i] for i in od]
    nq=(sub(sub(BASE,ONE),d[7]),sub(d[0],d[6]),sub(d[1],d[5]),sub(d[2],d[4]))
    gs=gaps(nq)
    gg=sorting_guards(gs,og)
    G=add(*[mul(4-2*i,gs[j]) for i,j in enumerate(og)])
    extra=dg+gg+[sub(nq[3],ONE),sub(G,g0),sub(add(g0,mul(8*(depth+1),ONE)),G)]
    return nq,extra,dg+gg


def root_domain(order,case):
    require(isinstance(order,str) and len(order)==5 and sorted(order)==list('01234'),
            'Full source gap order')
    require(type(case) is int and case in (1,2),'Residual branch')
    order=tuple(map(int,order))
    gs=gaps(INITIAL)
    g0=add(*[mul(2*i-4,gs[j]) for i,j in enumerate(order)])
    guards=[sub(BASE,mul(BASE_MINIMUM,ONE)),sub(sub(BASE,INITIAL[0]),ONE),
            *(sub(a,b) for a,b in zip(INITIAL,INITIAL[1:])),sub(INITIAL[-1],ONE)]
    guards.extend(sub(gs[j],gs[i]) for i,j in zip(order,order[1:]))
    residual=(sub(sub(BASE,gs[order[-1]]),mul(4,gs[order[0]]))
              if case==1 else sub(BASE,mul(6,gs[order[0]])))
    guards.append(sub(residual,mul(THRESHOLD+1,ONE)))
    return g0,guards


def expand_case(case):
    g0,guards=root_domain(case['order'],case['residual_case'])
    nodes=case['nodes']
    require(isinstance(nodes,list) and len(nodes)>0,'Nonempty tree')
    require(nodes[0]['parent'] is None,'Root parent')
    states=[(INITIAL,guards,0)]
    incoming=[0]*len(nodes)
    for i,node in enumerate(nodes):
        if i:
            parent=node['parent']
            require(type(parent) is int and 0<=parent<i,'Earlier parent')
            q,pg,depth=states[parent]
            nq,extra,_=affine_child(q,node['digits_order'],node['gaps_order'],g0,depth)
            states.append((nq,pg+extra,depth+1))
        require(states[i][2]<HORIZON,'Every node precedes the ninth update')
        children=node['children']
        require(isinstance(children,list) and len(children)==len(set(children)),
                'Unique child list')
        for child in children:
            require(type(child) is int and i<child<len(nodes)
                    and nodes[child]['parent']==i,'Valid forward tree edge')
            incoming[child]+=1
    require(incoming==[0]+[1]*(len(nodes)-1),'Complete connected tree')
    return g0,states


def coverage_query(case,expanded,ni,omit_child=None):
    g0,states=expanded
    q,guards,depth=states[ni]
    s=z3.SolverFor('QF_LIA')
    s.add(*[expr(g)>=0 for g in guards])
    B=VARS[0]
    u,v,w,x=map(expr,q)
    ds=[u,v,w,x-1,B-x-1,B-w-1,B-v-1,B-u]
    for k,(i,j) in enumerate(NETWORK):
        a,b=ds[i],ds[j]
        h,l=z3.Ints(f'h{k} l{k}')
        s.add(h==z3.If(a>=b,a,b),l==a+b-h)
        ds[i],ds[j]=h,l
    nq=(B-1-ds[7],ds[0]-ds[6],ds[1]-ds[5],ds[2]-ds[4])
    gs=(B-nq[0],nq[0]-nq[1],nq[1]-nq[2],nq[2]-nq[3],nq[3])
    G=sum(z3.If(a>=b,a-b,b-a) for i,a in enumerate(gs) for b in gs[i+1:])
    s.add(nq[3]>=1,G>=expr(g0),G<=expr(g0)+8*(depth+1))
    for ci in case['nodes'][ni]['children']:
        if ci==omit_child:
            continue
        c=case['nodes'][ci]
        _,_,sorting=affine_child(q,c['digits_order'],c['gaps_order'],g0,depth)
        s.add(z3.Not(z3.And(*[expr(v)>=0 for v in sorting])))
    return s


def inventory(certificate):
    require(certificate['format']=='nine-clock-affine-trees-v1','Certificate format')
    require((certificate['base_minimum'],certificate['threshold'],certificate['horizon'])
            ==(BASE_MINIMUM,THRESHOLD,HORIZON),'Declared universal assertion')
    cases=certificate['cases']
    require(tuple((c['order'],c['residual_case']) for c in cases)==KEYS,
            'Complete ordered source-case inventory')
    return cases


def solve_smt(smt):
    s=z3.SolverFor('QF_LIA')
    s.set(timeout=60000)
    s.add(z3.parse_smt2_string(smt.replace('(check-sat)','')))
    t=time.monotonic()
    zs=str(s.check())
    zt=time.monotonic()-t
    ys,yt=yices_status(smt,60)
    return dict(z3=zs,yices=ys,z3_seconds=zt,yices_seconds=yt,
                sha256=hashlib.sha256(smt.encode()).hexdigest())


def audit_case(case):
    expanded=expand_case(case)
    records=[]
    for ni in range(len(case['nodes'])):
        r=solve_smt(coverage_query(case,expanded,ni).to_smt2())
        records.append(dict(node=ni,**r))
    return dict(order=case['order'],residual_case=case['residual_case'],
                nodes=len(records),records=records)


def accept(records,cases):
    require(tuple((r['order'],r['residual_case']) for r in records)==KEYS,
            'Complete ordered audit inventory')
    for r,c in zip(records,cases):
        require(r['nodes']==len(c['nodes']) and
                tuple(x['node'] for x in r['records'])==tuple(range(len(c['nodes']))),
                'Every node independently checked')
        require(all(x['z3']==x['yices']=='unsat' for x in r['records']),
                'Every coverage obligation UNSAT in both engines')


def controls(cases):
    # Omit one nonempty root cell: coverage must become SAT.
    c=next(c for c in cases if c['order']=='12340' and c['residual_case']==1)
    expanded=expand_case(c)
    require(c['nodes'][0]['children'],'Nontrivial control root')
    ci=c['nodes'][0]['children'][0]
    solver=coverage_query(c,expanded,0,omit_child=ci)
    smt=solver.to_smt2()
    r=solve_smt(smt)
    require(r['z3']==r['yices']=='sat','Missing-cell control must be SAT')
    require(str(solver.check())=='sat','Recover coverage-hole witness')
    m=solver.model()
    B,u,v,w,x=[m.eval(v,model_completion=True).as_long() for v in VARS]
    q=(u,v,w,x)
    target=subtract_pair(B,q)
    require(q[-1]>=1 and target[-1]>=1 and terminal_residual(B,q)>THRESHOLD
            and dispersion(B,target)>=dispersion(B,q),'Independent coverage-hole replay')
    # Eight updates do not suffice: a fixed member of the infinite exact family.
    t=100
    base=34*t-9
    q=(26*t-12,20*t-9,14*t-6,8*t-3)
    path=[q]
    for _ in range(8):
        path.append(subtract_pair(base,path[-1]))
    require(terminal_residual(base,q)>THRESHOLD and
            all(v[-1]>=1 and dispersion(base,v)>=dispersion(base,q) for v in path),
            'Independent eight-step control replay')
    Bv,qs,D=full_prefix(8)
    s=z3.SolverFor('QF_LIA')
    s.add(*D,Bv==base,*[a==b for a,b in zip(qs[0],q)])
    from nine_zero_return_audit import symbolic_dispersion
    s.add(*[z[3]>=1 for z in qs],
          *[symbolic_dispersion(Bv,z)>=symbolic_dispersion(Bv,qs[0]) for z in qs[1:]])
    short=solve_smt(s.to_smt2())
    require(short['z3']==short['yices']=='sat','False eight-step claim must be SAT')
    return [dict(key='missing_root_cell',base=B,source=(u,v,w,x),target=target,**r),
            dict(key='false_eight_step_horizon',base=base,source=q,path=path,**short)]


def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is not accepted.')
    raw=(ROOT/'nine_clock_certificates.json').read_bytes()
    certificate=json.loads(raw)
    cases=inventory(certificate)
    # Establish the only dynamical hint before any tree obligation.
    dependency=solve_smt(dependency_obligation('dispersion_increase')[2])
    require(dependency['z3']==dependency['yices']=='unsat','One-step +8 dependency')
    records=[]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(audit_case,c) for c in cases]):
            r=f.result()
            records.append(r)
            print(dict(done=len(records),order=r['order'],residual_case=r['residual_case'],
                       nodes=r['nodes'],failed=sum(x['z3']!='unsat' or x['yices']!='unsat'
                                                for x in r['records']),
                       seconds=time.monotonic()-start),flush=True)
            (ROOT/'nine_clock_audit_checkpoint.json').write_text(json.dumps(records,indent=2)+'\n')
    records.sort(key=lambda r:KEYS.index((r['order'],r['residual_case'])))
    accept(records,cases)
    negative=controls(cases)
    rejected=0
    for kind in ('missing_case','missing_node','sat','unknown'):
        bad=json.loads(json.dumps(records))
        if kind=='missing_case':
            bad.pop()
        elif kind=='missing_node':
            bad[0]['records'].pop()
        else:
            bad[0]['records'][0]['yices']=kind
        try:
            accept(bad,cases)
        except ArithmeticError:
            rejected+=1
        else:
            raise ArithmeticError('Corrupt acceptance record admitted')
    n=sum(r['nodes'] for r in records)
    result=dict(accepted=True,cases=len(cases),nodes=n,records=records,
                dependency=dependency,negative_controls=negative,
                corruptions_rejected=rejected,solver_calls=2*n+7,
                certificate_sha256=hashlib.sha256(raw).hexdigest(),
                z3_version=z3.get_version_string(),
                yices_version=subprocess.run(['yices-smt2','--version'],capture_output=True,
                                             text=True,check=True).stdout.splitlines()[0],
                scope='For B>=40, every ten-state positive prefix with terminal residual >100 has a strict dispersion decrease among updates one through nine.')
    (ROOT/'nine_clock_audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'nine_clock_audit_checkpoint.json').unlink()
    print(f'PASS: {len(cases)} source cases and {n} affine coverage nodes in both solvers; two SAT controls; four rejected corruptions.',flush=True)


if __name__=='__main__':
    main()
