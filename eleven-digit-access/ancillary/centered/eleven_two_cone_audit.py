"""Unbounded independent audit of the two-fold return and native-gap clock.

The last step is explicitly computer-assisted: fifty transverse integer
states, each with two unrestricted longitudinal coordinates in a quadrant.
No bounded longitudinal box or inherited full centered-clock tree is used.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import z3

ROOT=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def mx(a,b):return z3.If(a>=b,a,b)
def mn(a,b):return z3.If(a<=b,a,b)


def symbolic_step(v):
    a,w,p,s,r=v
    C=mx(0,-a+2*s+1);A=mx(-1,-r)
    J=mx(0,-r-1)-mx(0,r-1);K=mx(0,r+1)
    H=mx(0,w-J-p-1);b=s-p-C;D=mx(A-1,b)
    E=mx(A,b-1)-mx(-b-1,-A);F=mn(A,b-1)+J+1
    return ((a-6*D,w+J-p-4*D,E-2*D,F-3*D+H,-p-K+D-H),
            r+J+K+H)


def guarded(v):
    a,w,p,s,r=v
    return z3.And(w-a+p<=-1,a-2*w-2*mx(0,r+1)>=-1)


def reduced(v):
    a,w,p,s,r=v;b=s-p
    return z3.And(a>=-27,w<=2,p>=-2,p<=2,b>=-4,b<=-1,r>=-b-4,r<=0)


def offsets(v):
    a,w,p,s,r=v
    return (-a,-2*a+6*p,-3*a+6*s,6*w-4*a,a+6*r)


def raw_step(x):
    pairs=[-6]
    for i in range(5):
        aa=x[i]-6*(i==4);bb=-x[4-i]-6*(i<4)
        pairs.extend((mx(aa,bb),mn(aa,bb)))
    return tuple(pairs[i]-pairs[-1-i] for i in range(5))


def query(key):
    real=key in ('primitive','preparation')
    S=z3.SolverFor('QF_LRA' if real else 'QF_LIA')
    S.set(timeout=30000)
    v=(z3.Reals if real else z3.Ints)('a0 w0 p0 s0 r0')
    source=v
    if key=='primitive':
        target,current=symbolic_step(v)
        raw=raw_step(raw_step(offsets(v)))
        xx=offsets(v);tt=offsets(target)
        S.add(guarded(v),z3.Or(*[a!=b for a,b in zip(raw,tt)],
                              tt[3]-tt[4]-xx[3]+xx[4]!=6*current))
    elif key=='preparation':
        target,_=symbolic_step(v)
        S.add(target[2]>2)
    else:
        h=3 if key in ('failure_reduction','false_two_observations') else 18
        if key=='false_two_observations':h=2
        if key.startswith('quadrant_'):
            p,b,r=map(int,key.split('_')[1:])
            S.add(v[0]>=-27,v[1]<=2,v[2]==p,v[3]==p+b,v[4]==r)
            threshold=2
        else:
            S.add(v[2]<=2);threshold=1
        work=0;states=[v]
        for j in range(h):
            S.add(guarded(v))
            target,current=symbolic_step(v)
            nxt=z3.Ints(f'a{j+1} w{j+1} p{j+1} s{j+1} r{j+1}')
            S.add(*[a==b for a,b in zip(nxt,target)])
            v=nxt;states.append(v);work+=current
            S.add(work<threshold)
        if key in ('failure_reduction','false_two_observations'):
            S.add(z3.Not(reduced(states[1])))
    return S,source


def main():
    require(__debug__,'Optimized Python is not an accepted audit mode')
    exe=shutil.which('yices-smt2');require(exe is not None,'Independent Yices required')
    keys=['primitive','preparation','failure_reduction','false_two_observations']
    keys += [f'quadrant_{p}_{b}_{r}' for p in range(-2,3) for b in range(-4,0)
             for r in range(-b-4,1)]
    require(len(keys)==54 and len(set(keys))==54,'Complete quadrant inventory')
    queries={};rows=[]
    for key in keys:
        S,source=query(key)
        smt='(set-logic '+('QF_LRA' if key in ('primitive','preparation') else 'QF_LIA')+')\n'+S.to_smt2()
        queries[key]=smt;status=str(S.check())
        run=subprocess.run([exe,'--timeout=30'],input=smt,text=True,capture_output=True,timeout=35)
        expected='sat' if key=='false_two_observations' else 'unsat'
        require(run.returncode==0 and status==run.stdout.strip()==expected,
                'Both engines must confirm '+key+': '+status+' / '+run.stdout+' / '+run.stderr)
        row=dict(key=key,z3=status,yices=run.stdout.strip(),
                 sha256=hashlib.sha256(smt.encode()).hexdigest())
        if status=='sat':
            from eleven_two_cone import guards,step,reduced as numerical_reduced,offsets as numeric_offsets
            from odd_base_covariance import centered_step
            M=S.model();v=tuple(int(str(M.eval(x))) for x in source)
            row['source']=v;total=0;first=None;x=numeric_offsets(v)
            for _ in range(2):
                require(guards(v),'SAT witness guards')
                v,aux=step(v);total+=aux['current'];x=centered_step(centered_step(x))
                require(x==numeric_offsets(v) and total<=0,'Independent primitive SAT replay')
                if first is None:first=v
            require(not numerical_reduced(first),'Two observations do not force the reduction')
            row['independent_witness_replay']=True
        rows.append(row)
    def accept(candidate):
        require(len(candidate)==len(keys) and {r['key'] for r in candidate}==set(keys),
                'Complete unique record inventory required')
        for row in candidate:
            expected='sat' if row['key']=='false_two_observations' else 'unsat'
            require(row['z3']==row['yices']==expected,'Unknown or wrong solver status')
            require(row['sha256']==hashlib.sha256(queries[row['key']].encode()).hexdigest(),
                    'Wrong query hash')
    accept(rows);rejected=[]
    for kind in ('missing','unknown','wrong_status','wrong_hash','duplicate'):
        altered=json.loads(json.dumps(rows))
        if kind=='missing':altered.pop()
        elif kind=='unknown':altered[0]['yices']='unknown'
        elif kind=='wrong_status':altered[0]['z3']='sat'
        elif kind=='wrong_hash':altered[0]['sha256']='incorrect'
        else:altered[1]=altered[0]
        try:accept(altered)
        except ArithmeticError:rejected.append(kind)
        else:raise ArithmeticError('Corrupt audit accepted')
    report=dict(status='passed',results=rows,solver_calls=2*len(rows),
                rejected_controls=rejected,
                scope='Unbounded primitive two-fold identity, preparation, failure reduction and all fifty longitudinal quadrants. This is a guarded return clock, not global centered composition.')
    (ROOT/'eleven_two_cone_obligations.json').write_text(json.dumps(queries,indent=2)+'\n')
    (ROOT/'eleven_two_cone_audit_results.json').write_text(json.dumps(report,indent=2)+'\n')
    summary=dict(status='passed',queries=len(rows),solver_calls=2*len(rows),
                 quadrants=50,rejected_controls=rejected)
    text=json.dumps(summary,indent=2)+'\n'
    (ROOT/'eleven_two_cone_audit_verification.txt').write_text(text)
    print(text)


if __name__=='__main__':main()
