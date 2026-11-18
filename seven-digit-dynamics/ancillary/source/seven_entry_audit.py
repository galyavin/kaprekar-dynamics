"""Independent Z3/CVC5 checks of finite real-affine leading-entry lemmas.

Each obligation is the negation of one precisely stated universal lemma.
The global almost-everywhere argument is separate, in the paper.
No solver call is hidden in a Python assert; non-UNSAT results fail closed.
"""
from pathlib import Path
import hashlib
import json
import time

import z3
import cvc5


def require(ok,context):
    if not ok:raise RuntimeError(context)


def maximum(x,y):return z3.If(x>=y,x,y)
def minimum(x,y):return z3.If(x<=y,x,y)
def absolute(x):return maximum(x,-x)


def leading(q):
    a,b,c=q
    values=(a,absolute(a+b-1),absolute(b+c-1))
    hi=maximum(values[0],maximum(values[1],values[2]))
    lo=minimum(values[0],minimum(values[1],values[2]))
    return hi,sum(values)-hi-lo,lo


def defect(q):
    a,b,c=q
    return sum(maximum(0,t) for t in
               (1-a-b,b+c-a,2*a-1-b-2*c,2*b+c+a-2))


def two_cycle(q):
    a,b,c=q
    return z3.And(a>=z3.RealVal('1/2'),a<=z3.RealVal('3/5'),
                  c==2*a-1,z3.Or(b==2-3*a,b==a))


def obligations():
    a,b,c=z3.Reals('a b c');q=(a,b,c)
    ordered=[a>=b,b>=c,c>=0,a<1]
    states=[q]
    for _ in range(6):states.append(leading(states[-1]))
    r1=states[1];r2=states[2];r3=states[3];r6=states[6]
    half=z3.RealVal('1/2');three_fifths=z3.RealVal('3/5')
    w0=defect(q);w3=defect(r3);w6=defect(r6)
    yield 'low_preparation', [*ordered,a<half,r1[0]<half,2*r1[0]+r1[1]<1]
    r=a-b;s=a+b+c-1
    expected=(a+maximum(-s,0),a-maximum(s,0),1-a-b)
    yield 'low_recurrence', [*ordered,a<half,2*a+b>=1,
                            z3.Or(*(x!=y for x,y in zip(r1,expected)))]
    yield 'three_step_defect', [*ordered,a>=half,w3>w0]
    yield 'high_top_stabilization', [*ordered,a>three_fifths,r3[0]>r2[0]]
    decrement=minimum(5*a-3,1-a)
    yield 'high_plateau_descent', [*ordered,a>three_fifths,r6[0]==a,
                                  w6>maximum(0,w0-decrement)]
    yield 'six_step_equality', [*ordered,a>half,r6[0]<1,w0>0,w6==w0,
                               z3.Not(two_cycle(q))]
    yield 'near_core_absorption', [*ordered,a>half,
                                  w0<half*minimum(2*a-1,1-a),w6>0]
    near=[absolute(b-(2-3*a))<(2*a-1)/100,
          absolute(c-(2*a-1))<(2*a-1)/100,
          z3.Or(b!=2-3*a,c!=2*a-1)]
    trapped=[]
    for state in states[:5]:
        trapped.extend((state[0]<=three_fifths,defect(state)>=2*state[0]-1))
    yield 'two_cycle_neighborhood', [*ordered,a>half,a<=three_fifths,*near,*trapped]


EXPECTED_KEYS=('low_preparation','low_recurrence','three_step_defect',
               'high_top_stabilization','high_plateau_descent',
               'six_step_equality','near_core_absorption','two_cycle_neighborhood')


def cvc_result(smt):
    solver=cvc5.Solver()
    solver.setLogic('QF_LRA')
    solver.setOption('tlimit-per','45000')
    parser=cvc5.InputParser(solver)
    parser.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6,
                          smt.replace('(check-sat)',''), 'entry-obligation')
    manager=parser.getSymbolManager()
    while True:
        command=parser.nextCommand()
        if command.isNull():break
        command.invoke(solver,manager)
    return str(solver.checkSat())


def accept(records):
    require(tuple(r['key'] for r in records)==EXPECTED_KEYS,'Complete obligation inventory required.')
    for record in records:
        require(record['z3']=='unsat' and record['cvc5']=='unsat',
                ('Non-UNSAT obligation',record['key'],record['z3'],record['cvc5']))


def main():
    root=Path(__file__).resolve().parent
    records=[];statements=[]
    for key,constraints in obligations():
        solver=z3.SolverFor('QF_LRA');solver.set(timeout=45000)
        solver.add(*constraints)
        smt=solver.to_smt2()
        start=time.monotonic();z3_status=str(solver.check());z3_seconds=time.monotonic()-start
        start=time.monotonic();cvc_status=cvc_result(smt);cvc_seconds=time.monotonic()-start
        record={'key':key,'z3':z3_status,'cvc5':cvc_status,
                'z3_seconds':z3_seconds,'cvc5_seconds':cvc_seconds,
                'sha256':hashlib.sha256(smt.encode()).hexdigest()}
        records.append(record);statements.append({'key':key,'negated_lemma_smtlib':smt})
        print(key,z3_status,cvc_status,flush=True)
        require(z3_status=='unsat' and cvc_status=='unsat',('Audit failed',record))
    accept(records)
    # Real SAT control in each engine, then acceptance controls for SAT,
    # UNKNOWN and a missing obligation. These never overwrite proof records.
    sat=z3.SolverFor('QF_LRA');sat.add(z3.BoolVal(True))
    require(str(sat.check())=='sat' and cvc_result(sat.to_smt2())=='sat','SAT controls')
    rejected=0
    for mutation in ('sat','unknown','missing'):
        bad=[dict(r) for r in records]
        if mutation=='missing':bad.pop()
        else:bad[0]['cvc5']=mutation
        try:accept(bad)
        except RuntimeError:rejected+=1
        else:raise RuntimeError('Invalid audit record accepted.')
    result={'logic':'QF_LRA','z3_version':z3.get_version_string(),
            'cvc5_version':cvc5.__version__,'records':records,
            'actual_solver_calls':2*len(records)+2,'rejected_record_controls':rejected,
            'scope':'Finite real-affine leading-entry lemmas. The measure-zero and compact-exhaustion arguments require the accompanying mathematical proof.'}
    (root/'seven_entry_obligations.json').write_text(json.dumps(statements,indent=2)+'\n')
    (root/'seven_entry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    summary=f'{len(records)} universal real-affine obligations UNSAT in both Z3 and CVC5; {result["actual_solver_calls"]} actual solver calls; {rejected} rejected acceptance controls.\n'
    (root/'seven_entry_verification.txt').write_text(summary)
    print(summary,end='')


if __name__=='__main__':main()
