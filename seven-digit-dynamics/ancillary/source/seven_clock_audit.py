"""Z3/Yices audit of actual seven-slot gap-potential descent.

The +8 one-step bound is a separately audited dependency of the partitioned
nine-step queries. Every state is an actual positive borrowed prefix. Source
gap sign masks partition the entire domain, including equality boundaries.
"""
from concurrent.futures import ProcessPoolExecutor,as_completed
from itertools import product
from pathlib import Path
import hashlib,json,time,subprocess,shutil
import z3,cvc5
from odd_gap import subtract_pair
from seven_clock import potential as numeric_potential

ROOT=Path(__file__).resolve().parent
NETWORK=((1,2),(4,5),(0,2),(3,5),(0,1),(3,4),(2,5),(0,3),(1,4),(2,4),(1,3),(2,3))
KEYS=('one_step_increase',)+tuple('nine_step_'+''.join(map(str,m)) for m in product((0,1),repeat=4))

def require(ok,label):
    if not ok:raise ArithmeticError(label)


def potential(B,q):
    gs=(B-q[0],q[0]-q[1],q[1]-q[2],q[2])
    return sum(z3.If(4*g>=B,4*g-B,B-4*g) for g in gs)


def prefix(depth,positive_final=True):
    B=z3.Int('B');q=z3.Ints('q0_0 q0_1 q0_2')
    D=[B>=2,B>q[0],q[0]>=q[1],q[1]>=q[2],q[2]>=1];states=[q]
    for t in range(depth):
        u,v,w=states[-1];ds=[u,v,w-1,B-w-1,B-v-1,B-u]
        for k,(i,j) in enumerate(NETWORK):
            a,b=ds[i],ds[j];hi=z3.Int(f'hi_{t}_{k}');lo=z3.Int(f'lo_{t}_{k}')
            D.extend((hi==z3.If(a>=b,a,b),lo==a+b-hi));ds[i],ds[j]=hi,lo
        nq=z3.Ints(f'q{t+1}_0 q{t+1}_1 q{t+1}_2')
        D.extend((nq[0]==B-1-ds[5],nq[1]==ds[0]-ds[4],nq[2]==ds[1]-ds[3]))
        D.extend((nq[0]<B,nq[0]>=nq[1],nq[1]>=nq[2],nq[2]>=0))
        if positive_final or t<depth-1:D.append(nq[2]>=1)
        states.append(nq)
    return B,states,D


def obligation(key):
    if key=='one_step_increase':
        B,qs,D=prefix(1,False);D.append(potential(B,qs[1])>potential(B,qs[0])+8)
    else:
        mask=tuple(map(int,key.removeprefix('nine_step_')))
        B,qs,D=prefix(9);ps=[potential(B,q) for q in qs]
        q=qs[0];g=(B-q[0],q[0]-q[1],q[1]-q[2],q[2])
        D.extend(4*t>=B if bit else 4*t<B for bit,t in zip(mask,g))
        D.extend((ps[0]>200,ps[9]>=ps[0]))
        # Consequences of the independently proved +8 bound. The lower
        # hints also use the negated target ps[9]>=ps[0].
        for i in range(10):
            D.extend((ps[i]<=ps[0]+8*i,ps[i]>=ps[0]-8*(9-i)))
            if i:D.append(ps[i]<=ps[i-1]+8)
    s=z3.SolverFor('QF_LIA');s.add(*D)
    return B,qs,s.to_smt2()


def cvc_status(smt,timeout=180000):
    s=cvc5.Solver();s.setLogic('QF_LIA');s.setOption('tlimit-per',str(timeout))
    pa=cvc5.InputParser(s);pa.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6,smt.replace('(check-sat)',''),'seven-clock');sm=pa.getSymbolManager()
    while True:
        c=pa.nextCommand()
        if c.isNull():break
        c.invoke(s,sm)
    start=time.monotonic();status=str(s.checkSat());return status,time.monotonic()-start


def yices_status(smt,timeout=180):
    executable=shutil.which('yices-smt2')
    require(executable is not None,'Yices 2 command yices-smt2 is required')
    start=time.monotonic()
    try:
        r=subprocess.run([executable,f'--timeout={timeout}'],input='(set-logic QF_LIA)\n'+smt,text=True,capture_output=True,timeout=timeout+5)
    except subprocess.TimeoutExpired:return 'unknown (process timeout)',time.monotonic()-start
    status=r.stdout.strip()
    require(r.returncode==0 and status in ('sat','unsat','unknown'),'Unexpected Yices result: '+r.stdout+r.stderr)
    return status,time.monotonic()-start


def solve(item):
    key,smt=item;s=z3.SolverFor('QF_LIA');s.set(timeout=180000);s.add(z3.parse_smt2_string(smt.replace('(check-sat)','')))
    t=time.monotonic();zs=str(s.check());zt=time.monotonic()-t
    cs,ct=yices_status(smt)
    return {'key':key,'z3':zs,'yices':cs,'z3_seconds':zt,'yices_seconds':ct,'sha256':hashlib.sha256(smt.encode()).hexdigest()}


def accept(records):
    require(tuple(r['key'] for r in records)==KEYS,'Complete ordered obligation inventory')
    require(all(r['z3']==r['yices']=='unsat' for r in records),'Every obligation must be UNSAT in both solvers')


def controls():
    records=[]
    for key,depth,threshold,bound in [('false_increase_four',1,None,4),('false_six_step_descent',6,1000,None)]:
        B,qs,D=prefix(depth);p0=potential(B,qs[0]);pt=potential(B,qs[-1]);s=z3.SolverFor('QF_LIA');s.set(timeout=45000);s.add(*D)
        if threshold is None:s.add(pt>p0+bound)
        else:s.add(p0>threshold,pt>=p0)
        smt=s.to_smt2();zs=str(s.check());cs,_=yices_status(smt,45)
        require(zs==cs=='sat','False strengthening requires SAT in both engines')
        m=s.model();base=m.eval(B).as_long();path=[tuple(m.eval(x).as_long() for x in q) for q in qs]
        require(all(subtract_pair(base,a)==b for a,b in zip(path,path[1:])),'Independent false-claim replay')
        ps=[numeric_potential(base,q) for q in path]
        require(ps[-1]>ps[0]+bound if threshold is None else ps[0]>threshold and ps[-1]>=ps[0],'Numeric violation')
        records.append({'key':key,'z3':zs,'yices':cs,'base':base,'path':path,'potential':ps,'smtlib':smt})
    return records


def main():
    if not __debug__:raise RuntimeError('Optimized Python is not accepted.')
    statements=[{'key':key,'negated_lemma_smtlib':obligation(key)[2]} for key in KEYS]
    (ROOT/'seven_clock_obligations.json').write_text(json.dumps(statements,indent=2)+'\n')
    # Accept the hint lemma before any dependent query is evaluated.
    first=solve((KEYS[0],statements[0]['negated_lemma_smtlib']));require(first['z3']==first['yices']=='unsat','One-step dependency failed');records=[first];print(first,flush=True)
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures=[ex.submit(solve,(r['key'],r['negated_lemma_smtlib'])) for r in statements[1:]]
        for f in as_completed(futures):
            r=f.result();records.append(r);print(r,flush=True)
            (ROOT/'seven_clock_audit_checkpoint.json').write_text(json.dumps(records,indent=2)+'\n')
    records.sort(key=lambda r:KEYS.index(r['key']));accept(records)
    negative=controls();rejected=0
    for kind in ('missing','sat','unknown'):
        bad=[dict(r) for r in records]
        if kind=='missing':bad.pop()
        else:bad[0]['yices']=kind
        try:accept(bad)
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Invalid acceptance record admitted')
    data={'records':records,'negative_controls':negative,'corruptions_rejected':rejected,'actual_solver_calls':38,'z3_version':z3.get_version_string(),'yices_version':subprocess.run(['yices-smt2','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],'scope':'Actual positive seven-slot prefixes: one-step potential increase at most eight, and strict nine-step decrease above 200. Global first-capture and weighted moment deductions are proved in the paper.'}
    (ROOT/'seven_clock_audit_results.json').write_text(json.dumps(data,indent=2)+'\n')
    (ROOT/'seven_clock_audit_checkpoint.json').unlink()
    print('PASS: seventeen universal integer obligations in both solvers; two independently replayed false strengthenings; three rejected record corruptions.',flush=True)
if __name__=='__main__':main()
