"""Universal integer-prefix exclusions, checked in both Z3 and CVC5."""
import hashlib,json,time
from pathlib import Path
import z3,cvc5
from odd_zero_access import symbolic_positive_prefix
from odd_gap import subtract_pair

ROOT=Path(__file__).resolve().parent
OBLIGATIONS=(('seven_odd_no_third_first_zero',3,2,8),
             ('nine_odd_no_fourth_first_zero',4,3,16))


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def cvc_result(smt):
    solver=cvc5.Solver();solver.setLogic('QF_LIA');solver.setOption('tlimit-per','45000')
    parser=cvc5.InputParser(solver)
    parser.setStringInput(cvc5.InputLanguage.SMT_LIB_2_6,smt.replace('(check-sat)',''),'zero-access')
    manager=parser.getSymbolManager()
    while True:
        command=parser.nextCommand()
        if command.isNull():break
        command.invoke(solver,manager)
    return str(solver.checkSat())


def accept(records):
    require(tuple(r['key'] for r in records)==tuple(r[0] for r in OBLIGATIONS),
            'Complete zero-access obligation inventory')
    for r in records:require(r['z3']==r['cvc5']=='unsat','Both independent solvers must return UNSAT')


def solve(k,depth,minimum):
    a,B,states,D=symbolic_positive_prefix(k,depth,minimum)
    solver=z3.SolverFor('QF_LIA');solver.set(timeout=45000);solver.add(*D)
    smt=solver.to_smt2();start=time.time();z=str(solver.check());c=cvc_result(smt)
    return solver,B,states,smt,z,c,time.time()-start


def main():
    records=[];statements=[]
    for key,k,depth,minimum in OBLIGATIONS:
        solver,B,states,smt,z,c,seconds=solve(k,depth,minimum)
        print(key,z,c,seconds,flush=True)
        require(z==c=='unsat','Unbounded odd-base gate exclusion failed')
        records.append({'key':key,'k':k,'depth':depth,'minimum_half':minimum,
                        'z3':z,'cvc5':c,'seconds':seconds,
                        'sha256':hashlib.sha256(smt.encode()).hexdigest()})
        statements.append({'key':key,'negated_lemma_smtlib':smt})
    accept(records);controls=[]
    for k,depth,minimum in [(3,1,8),(4,2,16)]:
        solver,B,states,smt,z,c,seconds=solve(k,depth,minimum)
        require(z==c=='sat','False stronger cutoff must admit a counterexample in both solvers')
        model=solver.model();base=model.eval(B).as_long()
        path=[tuple(model.eval(v).as_long() for v in q) for q in states]
        for x,y in zip(path,path[1:]):require(subtract_pair(base,x)==y,'Independent SAT path replay')
        final=subtract_pair(base,path[-1]);path.append(final)
        require(all(q[-1]>0 for q in path[:-1]) and final[-1]==0,'Actual first zero at the claimed sharp depth')
        controls.append({'k':k,'base':base,'z3':z,'cvc5':c,'path':path,'seconds':seconds})
        statements.append({'key':f'false_stronger_cutoff_{2*k+1}', 'smtlib':smt})
        print('false_stronger_cutoff',2*k+1,z,c,flush=True)
    rejected=0
    for kind in ('missing','sat','unknown'):
        bad=[dict(r) for r in records]
        if kind=='missing':bad.pop()
        else:bad[0]['cvc5']=kind
        try:accept(bad)
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Corrupted acceptance record accepted')
    result={'records':records,'sharpness_controls':controls,'corruptions_rejected':rejected,
            'actual_solver_calls':8,'z3_version':z3.get_version_string(),'cvc5_version':cvc5.__version__,
            'scope':'First zero entry from positive initial states: at most two steps in odd seven-slot bases B>=17, three in odd nine-slot bases B>=33. No larger odd dimension cap is asserted.'}
    (ROOT/'odd_zero_access_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'odd_zero_access_obligations.json').write_text(json.dumps(statements,indent=2)+'\n')
    print('PASS: two universal integer exclusions, two independently replayed SAT controls, three rejected corruptions.')


if __name__=='__main__':main()
