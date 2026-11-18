"""Four universal nine-digit leading entry lemmas, checked with Z3 and cvc5."""
import hashlib,json,time
from pathlib import Path
from fractions import Fraction as Q
import z3,cvc5
from nine_entry import CORE_FACETS
from nine_leading_audit import symbolic
from nine_leading_checks import reference

ROOT=Path(__file__).resolve().parent
EXPECTED=('periodic_four_step_classification','exceptional_first_stabilization',
          'exceptional_transverse_persistence','core_neighborhood_affine_progression')


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def core(q):
    return z3.And(*[sum(c*v for c,v in zip(g,(*q,1)))>=0 for g in CORE_FACETS])


def periodic_union(q):
    a,b,c,d=q
    fixed=z3.And(*[v==z3.RealVal('1/3') for v in q])
    two=z3.And(a==z3.RealVal('3/5'),d==z3.RealVal('1/5'),
               z3.Or(z3.And(b==a,c==a),z3.And(b==d,c==d)))
    three=z3.And(a>=z3.RealVal('5/7'),a<=z3.RealVal('3/4'),z3.Or(
        z3.And(b==2*a-1,c==b,d>0,d<=1-a),
        z3.And(c==3*a-2,d==3-4*a,b>=1-a,b<2-2*a),
        z3.And(b==a,c+d==2-2*a,d>=0,d<1-a)))
    return z3.Or(core(q),fixed,two,three)


def doubled_gaps(q):
    a,b,c,d=q
    return (1-a,1-b,a-c,b-d,c+d)


def sorted_gaps(q,prefix,conditions):
    values=list(doubled_gaps(q))
    for i in range(1,5):
        for j in range(i,0,-1):
            x,y=z3.Reals(f'{prefix}_{i}_{j}_hi {prefix}_{i}_{j}_lo')
            conditions.extend((x==z3.If(values[j-1]>=values[j],values[j-1],values[j]),
                               y==values[j-1]+values[j]-x))
            values[j-1],values[j]=x,y
    return values


def classification(depth):
    q,x,D,E=symbolic(depth);conditions=[*D,*E,x[-1][0]<1,z3.Not(periodic_union(x[-1]))]
    conditions.extend(s[0]==q[0] for s in x[1:])
    g0=sorted_gaps(q,'g0',conditions)
    for j in range(1,depth+1):
        gj=sorted_gaps(x[j],f'g{j}',conditions)
        conditions.extend(v==w for v,w in zip(g0,gj))
    return x,conditions


def exceptional_tube(q,eta,stable=False):
    a,b,c,d=q
    return [eta>0,eta<=z3.RealVal('1/10000'),
            a>=z3.RealVal('5/7')-(0 if stable else eta),a<=z3.RealVal('3/4'),
            d>=100*eta,d<=1-a+eta,b-c<=eta,
            c-2*a+1>=-eta,c-2*a+1<=eta]


def transverse(q):
    a,b,c,d=q;s=c-2*a+1
    return b-c+z3.If(s>=0,s,-s)


def near_core_progression(preparation=3):
    depth=preparation+6;q,x,D,E=symbolic(depth);eta=z3.Real('eta')
    near=[eta>0,*[sum(c*v for c,v in zip(g,(*q,1)))>=-eta for g in CORE_FACETS],
          *[v>=100*eta for v in doubled_gaps(q)]]
    last=preparation+6;mid=preparation+3
    false=z3.Or(*[x[last][j]-x[mid][j]!=x[mid][j]-x[preparation][j] for j in range(4)])
    return x,[*D,*E,*near,z3.Not(core(x[-1])),false]


def obligations():
    x,D=classification(4);yield EXPECTED[0],x,D
    q,x,D,E=symbolic(6);eta=z3.Real('eta')
    yield EXPECTED[1],x,[*D,*E,*exceptional_tube(q,eta),x[6][0]!=x[3][0]]
    q,x,D,E=symbolic(6);eta=z3.Real('eta')
    yield EXPECTED[2],x,[*D,*E,*exceptional_tube(q,eta,True),x[6][0]==q[0],
                         transverse(x[6])<transverse(x[3])]
    x,D=near_core_progression();yield EXPECTED[3],x,D


def cvc_result(smt):
    solver=cvc5.Solver();solver.setLogic('QF_LRA');solver.setOption('tlimit-per','90000')
    parser=cvc5.InputParser(solver)
    parser.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6,smt.replace('(check-sat)',''),'nine-global-entry')
    manager=parser.getSymbolManager()
    while True:
        cmd=parser.nextCommand()
        if cmd.isNull():break
        cmd.invoke(solver,manager)
    return str(solver.checkSat())


def accept(records):
    require(tuple(r['key'] for r in records)==EXPECTED,'Complete global-entry inventory')
    for r in records:require(r['z3']==r['cvc5']=='unsat','Both engines must return UNSAT')


def solve(D):
    solver=z3.SolverFor('QF_LRA');solver.set(timeout=90000);solver.add(*D)
    smt=solver.to_smt2();t=time.time();z=str(solver.check());zt=time.time()-t
    t=time.time();c=cvc_result(smt);ct=time.time()-t
    return solver,smt,z,c,zt,ct


def main():
    records=[];statements=[];controls=[]
    for key,x,D in obligations():
        solver,smt,z,c,zt,ct=solve(D)
        record={'key':key,'z3':z,'cvc5':c,'z3_seconds':zt,'cvc5_seconds':ct,
                'sha256':hashlib.sha256(smt.encode()).hexdigest()}
        print(record,flush=True)
        require(z==c=='unsat','A global-entry obligation failed')
        records.append(record);statements.append({'key':key,'negated_lemma_smtlib':smt})
    accept(records)
    for key,(x,D) in [('false_three_step_classification',classification(3)),
                      ('false_unprepared_affine_progression',near_core_progression(0))]:
        solver,smt,z,c,zt,ct=solve(D);require(z==c=='sat','False strengthening requires SAT in both solvers')
        model=solver.model()
        path=[tuple(Q(model.eval(v).numerator_as_long(),model.eval(v).denominator_as_long()) for v in s) for s in x]
        for a,b in zip(path,path[1:]):require(reference(a)==b,'Independent nine-value control replay')
        controls.append({'key':key,'z3':z,'cvc5':c,'path':[[str(v) for v in s] for s in path],
                         'z3_seconds':zt,'cvc5_seconds':ct})
        statements.append({'key':key,'smtlib':smt});print(key,z,c,flush=True)
    rejected=0
    for kind in ('missing','sat','unknown'):
        bad=[dict(r) for r in records]
        if kind=='missing':bad.pop()
        else:bad[0]['cvc5']=kind
        try:accept(bad)
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Corrupted global-entry record accepted')
    result={'records':records,'controls':controls,'actual_solver_calls':12,
            'corruptions_rejected':rejected,'z3_version':z3.get_version_string(),
            'cvc5_version':cvc5.__version__,
            'scope':'Four universal finite real-affine lemmas. Global periodic classification and almost-everywhere finite core entry additionally use the analytical argument and preceding global majorization theorem.'}
    (ROOT/'nine_global_entry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'nine_global_entry_obligations.json').write_text(json.dumps(statements,indent=2)+'\n')
    print('PASS: four universal lemmas in both solvers, two independently replayed SAT controls, three rejected corruptions.')


if __name__=='__main__':main()
