"""Independent checks for the global seven-slot actual clock theorem."""
from itertools import permutations,product
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory
import json,shutil,subprocess,sys
import sympy as sp
import z3
from seven_clock_audit import NETWORK,prefix,accept,require
from seven_clock import potential,positive_bound,global_bound
from seven_zero_transience_audit import accept as accept_zero
from odd_gap import subtract_pair

ROOT=Path(__file__).resolve().parent

def raw_number(B,q):
    ds=[*q,0,0,0,0]
    hi=sum(v*B**j for j,v in enumerate(reversed(ds)));lo=sum(v*B**j for j,v in enumerate(ds));x=hi-lo;out=[]
    for _ in range(7):out.append(x%B);x//=B
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(3))


def six_step_family():
    t=sp.Symbol('t',integer=True);B=8*t+4
    qs=[(5*t+3,3*t+2,2*t+1),(6*t+3,3*t+1,2*t+1),(6*t+3,4*t+2,3*t+1),(6*t+2,3*t+3,t+1),(7*t+3,5*t,3*t-1),(7*t+2,4*t+5,2*t+1),(7*t+1,5*t+2,2*t+4)]
    for q,p in zip(qs,qs[1:]):
        u,v,w=q;ds=[u,v,w-1,B-1,B-w-1,B-v-1,B-u];order=sorted(range(7),key=lambda i:int(ds[i].subs(t,100)),reverse=True);sorted_ds=[ds[i] for i in order]
        for a,b in zip(sorted_ds,sorted_ds[1:]):
            z=sp.expand(a-b);require(z.coeff(t)>=0 and z.subs(t,10)>=0,'All-t family sorting guard')
        require(all(sp.expand(sorted_ds[i]-sorted_ds[-1-i]-p[i])==0 for i in range(3)),'All-t family update')
    for q in (qs[0],qs[-1]):
        gs=[B-q[0],q[0]-q[1],q[1]-q[2],q[2]];terms=[sp.expand(4*g-B) for g in gs];signs=[1 if z.subs(t,10)>=0 else -1 for z in terms]
        require(all((s*z).coeff(t)>=0 and (s*z).subs(t,10)>=0 for s,z in zip(signs,terms)),'Family absolute-value guard')
        require(sp.expand(sum(s*z for s,z in zip(signs,terms))-8*t)==0,'Unbounded family plateau')
    return {'base':'8t+4','source':['5t+3','3t+2','2t+1'],'minimum_t':10,'initial_and_sixth_potential':'8t'}


def main():
    if not __debug__:raise RuntimeError('Optimized Python is not accepted.')
    audit=json.loads((ROOT/'seven_clock_audit_results.json').read_text());accept(audit['records'])
    zero=json.loads((ROOT/'seven_zero_transience_results.json').read_text());accept_zero(zero['records'])
    network=0
    for source in list(permutations(range(6)))+list(product((0,1),repeat=6)):
        d=list(source)
        for i,j in NETWORK:d[i],d[j]=max(d[i],d[j]),min(d[i],d[j])
        require(d==sorted(source,reverse=True),'Sorting circuit');network+=1
    rng=Random(541116);numeric=0;encodings=0
    for B in (18,19,400,401,10000,10001,10**12,10**12+1):
        for _ in range(100):
            q=tuple(sorted([rng.randrange(1,B) for _ in range(3)],reverse=True))
            for _ in range(9):
                nxt=subtract_pair(B,q);require(nxt==raw_number(B,q),'Whole-integer primitive');numeric+=1
                if q[-1]>0:require(potential(B,nxt)<=potential(B,q)+8,'One-step potential')
                q=nxt
        if B in (18,19,400,401):
            for _ in range(12):
                q=tuple(sorted([rng.randrange(1,B) for _ in range(3)],reverse=True));base,states,D=prefix(1,False);s=z3.Solver();s.add(*D,base==B,*(a==b for a,b in zip(states[0],q)));require(s.check()==z3.sat,'Fixed-source encoding feasible');m=s.model();image=tuple(m.eval(a).as_long() for a in states[1]);require(image==raw_number(B,q),'Symbolic encoding against independent integer');encodings+=1
    family=six_step_family();finite=[]
    compiler=shutil.which('c++');require(compiler is not None,'C++17 compiler required')
    with TemporaryDirectory(prefix='seven-clock-') as tmp:
        binary=str(Path(tmp)/'holdout');subprocess.run([compiler,'-O2','-std=c++17',str(ROOT/'seven_clock_holdout.cpp'),'-o',binary],check=True)
        for B,prefixes in ((129,True),(130,True),(400,False),(401,False)):
            args=[binary,str(B)]+(['prefixes'] if prefixes else []);r=subprocess.run(args,check=True,capture_output=True,text=True);record=json.loads(r.stdout);finite.append(record);print('Whole graph',record,flush=True)
    optimized=subprocess.run([sys.executable,'-O',str(Path(__file__).resolve())],capture_output=True,text=True);require(optimized.returncode!=0 and 'Optimized Python' in optimized.stderr,'Disabled assertion control')
    result={'network_controls':network,'whole_integer_updates':numeric,'fixed_symbolic_sources':encodings,'unbounded_six_step_counterexample':family,'finite_graphs':finite,'graph_states':sum(r['states'] for r in finite),'column_updates':sum(r['raw_column_updates'] for r in finite),'optimized_python_rejected':True,'scope':'Independent primitive, circuit, affine-family and complete finite graph checks. Unbounded descent and zero clearance rely on the exact solver audits and analytical proof; finite graphs do not prove those statements.'}
    (ROOT/'seven_clock_check_results.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: global clock primitives, unbounded six-step counterexample, complete finite graphs and optimization rejection.',flush=True)
if __name__=='__main__':main()
