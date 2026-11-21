"""Solver-free B-facet proof replay and independent integer subtraction.

Optional --generate uses Z3 to find nonnegative rational coefficients. All
acceptance is exact arithmetic, with a complete expected guard inventory.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

from nine_entry import evaluate,compose,bulk_exit,inverse_leading
from nine_b_facet import (B0,NON_B,PHASES,DELTA_B,THREE_ORDERS,phase,
    phase_next,three_formula,six_orders,guard_domain,absorb,facet_exit,leading_exit)
from nine_phase_checks import replay,raw_step

ROOT=Path(__file__).resolve().parent
IDENTITY=tuple(tuple(int(i==j) for j in range(6)) for i in range(4))


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def from_orders(orders):
    from nine_phase_checks import step
    matrix=IDENTITY
    for order in orders:
        matrix,_=step(matrix,(0,0,0,0,1,0),(0,0,0,0,0,1),order)
    return {'orders':orders,'matrix':matrix}


def witness(row,domain):
    import z3
    variables=z3.Reals(' '.join('c'+str(i) for i in range(len(domain))))
    solver=z3.Solver();solver.set(timeout=10000)
    solver.add(*[v>=0 for v in variables])
    solver.add(*[sum(v*int(r[j]) for v,r in zip(variables,domain))==int(row[j])
                 for j in range(len(row))])
    require(solver.check()==z3.sat,'Witness search failed')
    model=solver.model();cs=[]
    for v in variables:
        x=model.eval(v,model_completion=True)
        cs.append(Q(x.numerator_as_long(),x.denominator_as_long()))
    return {'row':row,'coefficients':list(map(str,cs))}


def rational_replay(records,rows,domain):
    require([tuple(r['row']) for r in records]==list(map(tuple,rows)),
            'Complete rational obligation inventory')
    for record,row in zip(records,rows):
        cs=list(map(Q,record['coefficients']))
        require(len(cs)==len(domain) and min(cs)>=0,'Nonnegative coefficients')
        require(tuple(sum(c*r[j] for c,r in zip(cs,domain)) for j in range(len(row)))
                ==tuple(row),'Exact rational identity')
    return len(records)


def domains(orders):
    record=from_orders(orders);matrix,guards,_=replay(record)
    homogeneous=sorted(set(g[:5] for g in guards))
    equalities=(B0,tuple(-x for x in B0),(0,0,0,0,1))
    return record,matrix,guards,homogeneous,(*NON_B,*equalities),(*homogeneous,*equalities)


def generate():
    records=[]
    for key,orders in THREE_ORDERS.items():
        record=from_orders(orders);_,guards,_=replay(record)
        records.append({'key':key,'proofs':[witness(g,guard_domain(key)) for g in guards]})
    six=[]
    for e in PHASES:
        record,M,guards,hom,forward,backward=domains(six_orders(e))
        six.append({'phase':e,'forward':[witness(g,forward) for g in hom],
                    'backward':[witness(g,backward) for g in NON_B]})
    data={'records':records,'six':six}
    (ROOT/'nine_b_facet_certificates.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


def verify(data):
    require([r['key'] for r in data['records']]==list(THREE_ORDERS),'Three-word inventory')
    require([r['phase'] for r in data['six']]==list(PHASES),'Six-word inventory')
    count=0
    expected_high=((1,0,0,0,0,-3),(0,1,0,0,0,-2),
                   (0,0,1,0,0,-1),(0,0,0,1,0,-1))
    expected_low=((1,0,0,0,0,-3),(0,0,-1,0,1,-3),
                  (0,-1,0,0,1,-2),(0,0,0,1,0,-1))
    for record,expected in zip(data['records'],(expected_high,expected_low)):
        key=record['key'];M,guards,_=replay(from_orders(THREE_ORDERS[key]))
        require(M==expected,'Three-step affine map')
        count+=rational_replay(record['proofs'],guards,guard_domain(key))
    for record in data['six']:
        e=record['phase'];_,M,guards,hom,forward,backward=domains(six_orders(e))
        count+=rational_replay(record['forward'],hom,forward)
        count+=rational_replay(record['backward'],NON_B,backward)
        for i,row in enumerate(M):
            diff=tuple(row[j]-IDENTITY[i][j]-(DELTA_B[i] if j==5 else 0) for j in range(6))
            multiplier=diff[1]
            require(diff==tuple(multiplier*x for x in (*B0,-e)),
                    'Six-step displacement on the exact phase plane')
        require(tuple(evaluate(row[:4],DELTA_B) for row in M)==DELTA_B,
                'Repeated affine progression')
        require(phase_next(phase_next(e))==e,'Closed two-return phase')
    return count


def phase_checks():
    max_steps=max_error=max_facet_loss=0;worst=[]
    # This is an exhaustive proof of the finite integer phase domain. The
    # displacement formula is invariant under translations preserving e.
    for e in range(-512,513):
        B=1000000;q=(900000,620000,380000+e,200000);start=q;n=0
        while phase(B,q) not in PHASES:
            q=three_formula(B,q);n+=3
            require(abs(phase(B,q))<=512,'Invariant bounded phase interval')
            max_error=max(max_error,*(abs(a-b) for a,b in zip(q,start)))
            displacement=tuple(a-b for a,b in zip(q,start))
            max_facet_loss=max(max_facet_loss,*(-evaluate(r[:4],displacement) for r in NON_B))
            require(n<=513,'Finite absorption')
        max_steps=max(max_steps,n)
        if n==513:worst.append(e)
    max_norm=max(sum(abs(x) for x in r[:4]) for r in NON_B)
    require(max_steps<=513 and max_error<=513 and max_facet_loss<=4000,
            'Buffer remains at least 10000')
    return {'phases_checked':1025,'max_steps':max_steps,'max_coordinate_error':max_error,
            'facet_max_l1_norm':max_norm,'max_facet_loss':max_facet_loss,
            'worst_initial_phases':worst}


def determinant(matrix):
    if len(matrix)==1:return matrix[0][0]
    return sum((-1)**j*x*determinant([r[:j]+r[j+1:] for r in matrix[1:]])
               for j,x in enumerate(matrix[0]))


def clock_check():
    alpha=tuple(map(Q,('52/45','-23/180','67/18','-59/15')))
    decreases=((3,2,1,1),(6,12,6,4),(10,20,30,20),(2,4,6,-2))
    require([evaluate(alpha,d) for d in decreases]==[3,12,42,32],'Old clock identities')
    det=determinant(decreases);require(det!=0,'Unique old linear functional')
    residual=6+evaluate(alpha,DELTA_B)
    require(residual==Q(647,60),'New closed-phase residual')
    return {'old_decrease_determinant':det,'alpha':list(map(str,alpha)),
            'six_step_residual':str(residual)}


def repetition_error():
    errors=[]
    for e in PHASES:
        _,_,guards,*_=domains(six_orders(e))
        # Project onto B0 by replacing w with w-e. On that plane the proved
        # homogeneous word domain is C(K) intersect B0. Affine offsets and
        # integer rounding account for the entire repetition-count error.
        errors.append(1+max(Q(abs(g[2]*e+g[-1]),-evaluate(g[:4],DELTA_B))
                      for g in guards if evaluate(g[:4],DELTA_B)<0))
    require(max(errors)==5,'Uniform facet repetition error')
    return {'by_phase':dict(zip(map(str,PHASES),map(str,errors))),
            'repetitions':5,'elapsed_steps':30,'endpoint_linf':34}


def actual_checks():
    records=[];total=0
    for B in range(200000,200005):
        for e in PHASES:
            q=(9*B//10,31*B//50,B-31*B//50+e,B//5);start=q
            output=facet_exit(B,q)
            for _ in range(output['elapsed']):q=raw_step(B,q)
            require(q==output['target'],'Independent full facet path')
            total+=output['elapsed']
            records.append({'base':B,'phase':e,'start':start,**output})
    # Extremal and seam absorption phases, in five consecutive base regimes.
    for B in range(1000000,1000005):
        for e in (-512,-511,-5,-4,-3,-2,-1,0,1,510,511,512):
            q=(9*B//10,31*B//50,B-31*B//50+e,B//5);output=absorb(B,q)
            for _ in range(output['elapsed']):q=raw_step(B,q)
            require(q==output['target'],'Independent phase absorption')
            total+=output['elapsed']
    return records,total


def transfer_checks():
    data=json.loads((ROOT/'nine_entry_certificates.json').read_text())
    pieces=data['exit_partitions'];records=[]
    for piece in pieces:
        if not piece['face'].startswith('B'):continue
        rotation=int(piece['face'][1]);r=tuple(map(Q,piece['radius_center']))
        for signs in product((-1,1),repeat=4):
            x=inverse_leading(r,signs)
            for B in range(10000000,10000005):
                source=tuple(int(B*t) for t in x);q=raw_step(B,source)
                segment=bulk_exit(B,q);q=segment['target']
                for _ in range(rotation):q=raw_step(B,q)
                require(abs(phase(B,q))<=212,'Bounded incoming B offset')
                absorption=absorb(B,q);q=absorption['target'];facet=facet_exit(B,q)
                require(facet['repetitions']>0,'Positive facet access')
                records.append({'base':B,'face':piece['face'],'signs':signs,
                                'absorption_steps':absorption['elapsed'],
                                'phase':phase(B,q),'facet_steps':facet['elapsed']})
    require(len(records)==240,'All B faces, inverse branches and five bases')
    return records


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true')
    args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_b_facet_certificates.json').read_text())
    count=verify(data);controls=0
    for kind in ('missing','coefficient','phase'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'][0]['proofs'].pop()
        elif kind=='coefficient':bad['records'][0]['proofs'][0]['coefficients'][0]='-1'
        else:bad['six'][0]['phase']=0
        try:verify(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt proof accepted')
    phases=phase_checks();clock=clock_check();error=repetition_error()
    actual,total=actual_checks();transfer=transfer_checks()
    result={'status':'passed','rational_identities':count,'corruptions_rejected':controls,
            'absorption':phases,'clock':clock,'repetition_error':error,'independent_steps':total,
            'full_facet_segments':actual,'transfer_checks':transfer,
            'scope':'Timed B-facet paths and weighted access; no terminal-arrival or global law claim.'}
    (ROOT/'nine_b_facet_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {count} rational guard identities; {controls} corruptions rejected.\n'
            f'All 1025 incoming integer phases absorbed in at most {phases["max_steps"]} steps.\n'
            f'20 full facet paths and 60 absorption paths: {total} independent subtraction steps.\n'
            '240 incoming transfers across three faces, sixteen inverse branches and five bases.\n'
            'Unique old linear clock has new closed-phase residual 647/60.\n'
            'No first-capture, global-entry, full arrival-law or seven-digit tail conclusion.\n')
    (ROOT/'nine_b_facet_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
