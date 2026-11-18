"""Complete zero-face clearance in seven slots, over unbounded bases.

The separate seven_first_access check proves the two-update positive-access bound;
this audit closes zero-origin returns and the explicit rank-one images.
"""
from pathlib import Path
import hashlib,json,time
import z3
from seven_clock_audit import NETWORK,cvc_status,require
from odd_gap import subtract_pair

ROOT=Path(__file__).resolve().parent
KEYS=('even_first_positive','even_second_zero','even_no_third_return','odd_first_positive','odd_second_zero','odd_no_third_return','rank_one_formula')

def sorted_step(B,values,tag):
    ds=list(values);D=[]
    for k,(i,j) in enumerate(NETWORK):
        a,b=ds[i],ds[j];hi=z3.Int(f'{tag}_hi{k}');lo=z3.Int(f'{tag}_lo{k}')
        D.extend((hi==z3.If(a>=b,a,b),lo==a+b-hi));ds[i],ds[j]=hi,lo
    q=z3.Ints(f'{tag}_u {tag}_v {tag}_w')
    D.extend((q[0]==B-1-ds[5],q[1]==ds[0]-ds[4],q[2]==ds[1]-ds[3]))
    return q,D

def positive_step(B,q,tag):
    u,v,w=q;return sorted_step(B,[u,v,w-1,B-w-1,B-v-1,B-u],tag)

def statements():
    a,u,v=z3.Ints('a u v')
    for beta,label in ((0,'even'),(1,'odd')):
        B=2*a+beta;D=[a>=9,B>u,u>=v,v>=1]
        q1,d1=sorted_step(B,[u,v-1,B-1,B-1,B-v-1,B-u],label+'1')
        yield label+'_first_positive',D+d1+[q1[2]<=0]
        q2,d2=positive_step(B,q1,label+'2')
        special=z3.And(u==a,v==a) if beta==0 else z3.Or(z3.And(u==a,v==a),z3.And(u==a+1,v==a),z3.And(u==a+1,v==a+1))
        yield label+'_second_zero',D+d1+[q1[2]>=1]+d2+[z3.Or((q2[2]==0)!=special,z3.And(special,z3.Or(q2[0]!=a+1,q2[1]!=1)))]
        q3,d3=positive_step(B,q2,label+'3')
        yield label+'_no_third_return',D+d1+[q1[2]>=1]+d2+[q2[2]>=1]+d3+[q3[2]==0]
    B,u=z3.Ints('B u')
    q,d=sorted_step(B,[u-1,B-1,B-1,B-1,B-1,B-u],'rank1')
    hi=z3.If(B-u>=u-1,B-u,u-1);lo=B-1-hi
    yield 'rank_one_formula',[B>=18,B>u,u>=1]+d+[z3.Or(q[0]!=hi,q[1]!=lo,q[2]!=0)]

def accept(records):
    require(tuple(r['key'] for r in records)==KEYS,'Complete zero-return inventory')
    require(all(r['z3']==r['cvc5']=='unsat' for r in records),'Every zero identity must pass both engines')


def main():
    if not __debug__:raise RuntimeError('Optimized Python is not accepted.')
    records=[];formula=[]
    for key,D in statements():
        s=z3.SolverFor('QF_LIA');s.set(timeout=45000);s.add(*D);smt=s.to_smt2();t=time.monotonic();zs=str(s.check());zt=time.monotonic()-t;cs,ct=cvc_status(smt,45000)
        r={'key':key,'z3':zs,'cvc5':cs,'z3_seconds':zt,'cvc5_seconds':ct,'sha256':hashlib.sha256(smt.encode()).hexdigest()};records.append(r);formula.append({'key':key,'smtlib':smt});print(r,flush=True)
        require(zs==cs=='unsat','Zero-return identity audit failed')
    accept(records)
    rejected=0
    for kind in ('missing','sat','unknown'):
        bad=[dict(r) for r in records]
        if kind=='missing':bad.pop()
        else:bad[0]['cvc5']=kind
        try:accept(bad)
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Invalid zero audit record accepted')
    controls=[]
    for B in (18,19,1000,1001):
        q=(1,0,0);path=[q]
        for _ in range(3):q=subtract_pair(B,q);path.append(q)
        require(path[2][2]==0 and path[3][2]>0,'Rank-one endpoint takes three steps to positivity')
        controls.append({'base':B,'path':path})
    for B in (19,1001):
        q=(B//2+1,0,0);path=[q]
        for _ in range(4):q=subtract_pair(B,q);path.append(q)
        require(path[3][2]==0 and path[4][2]>0,'Four-step zero-origin clearance is sharp in odd bases')
        controls.append({'base':B,'path':path})
    (ROOT/'seven_zero_transience_obligations.json').write_text(json.dumps(formula,indent=2)+'\n')
    (ROOT/'seven_zero_transience_results.json').write_text(json.dumps({'records':records,'actual_solver_calls':14,'corruptions_rejected':rejected,'z3_version':z3.get_version_string(),'clearance_controls':controls,'scope':'Zero-origin first returns and four-step permanent clearance, using inherited positive first-access bounds. Every seven-slot orbit is permanently positive by update six for B>=18.'},indent=2)+'\n')
    print('PASS: seven zero-return identities in both solvers and six clearance controls.',flush=True)
if __name__=='__main__':main()
