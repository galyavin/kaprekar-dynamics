"""Exact six-cell input weights and independent full numerical arrivals.

Run in .venv-research for exact SymPy geometry. --generate uses Z3 only to
produce rational polyhedral reductions; default acceptance replays those
identities and performs two independent exact integration constructions.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

from nine_b_routes import route_specs,limit_summary,exact_arrival,MASS
from nine_entry import evaluate,compose,inverse_leading
from nine_entry_checks import geometry,cubature,folded_density,simplex_rule,integrate
from nine_b_facet_checks import require,witness,rational_replay

ROOT=Path(__file__).resolve().parent
CONSTANT=(0,0,0,0,1)


def reduce_domain(domain):
    rows=list(domain);records=[]
    for row in list(rows):
        others=[g for g in rows if g!=row];source=(*others,CONSTANT)
        try:proof=witness(row,source)
        except ArithmeticError:continue # Keeping a row requires no implication claim.
        rational_replay([proof],[row],source)
        records.append({'domain':source,'proof':proof});rows.remove(row)
    return tuple(rows),records


def generate():
    records=[]
    for spec in route_specs():
        domain,proofs=reduce_domain(spec['domain']);vertices,simplices=geometry(domain)
        records.append({**spec,'reduced_domain':domain,'reductions':proofs,
                        'vertices':vertices,'simplices':simplices})
    data={'routes':records}
    (ROOT/'nine_b_routes_certificates.json').write_text(json.dumps(data,indent=2,default=str)+'\n')
    return data


def rational_matrix(value):return tuple(tuple(map(Q,row)) for row in value)


def verify_geometry(data):
    specs=route_specs();nodes=simplex_rule();records=[];identities=0;total=Q(0)
    require([(r['first'],r['second']) for r in data['routes']]
            ==[(s['first'],s['second']) for s in specs],'Complete six-route inventory')
    for record,spec in zip(data['routes'],specs):
        require(tuple(map(tuple,record['domain']))==spec['domain'],'Full incoming inequalities')
        for key in ('bb_entry_matrix','bbb_entry_matrix'):
            require(rational_matrix(record[key])==spec[key],'Leading route matrix')
        rows=list(spec['domain'])
        for reduction in record['reductions']:
            row=tuple(reduction['proof']['row']);require(row in rows,'Reduced row is present')
            domain=tuple(g for g in rows if g!=row)+(CONSTANT,)
            require(tuple(map(tuple,reduction['domain']))==domain,'Sequential reduction domain')
            identities+=rational_replay([reduction['proof']],[row],domain);rows.remove(row)
        require(tuple(map(tuple,record['reduced_domain']))==tuple(rows),'Reduced facet inventory')
        vertices,simplices=geometry(tuple(rows))
        require(rational_matrix(record['vertices'])==tuple(vertices),'Exact vertex enumeration')
        require(tuple(map(tuple,record['simplices']))==tuple(simplices),'Exact triangulation')
        center=tuple(sum(v[i] for v in vertices)/len(vertices) for i in range(4))
        require(min(evaluate(g,(*center,1)) for g in spec['domain'])>0,'Full-dimensional strict route')
        target=spec['bbb_entry_matrix'];base=CONSTANT
        require(target[0]==tuple(3*v-b for v,b in zip(target[1],base)), 'Triple-B first coordinate')
        require(target[2]==tuple(b-v for b,v in zip(base,target[1])), 'Triple-B third coordinate')
        require(target[3]==tuple(2*v-b for v,b in zip(target[1],base)), 'Triple-B fourth coordinate')
        for v in vertices:
            t=evaluate(target[1],(*v,1));require(Q(3,5)<=t<=Q(2,3),'Terminal interval on every vertex')
        t=evaluate(target[1],(*center,1));require(Q(3,5)<t<Q(2,3),'Strict terminal-strip access')
        quadrature=cubature(vertices,simplices,nodes,folded_density)
        expanded=sum(integrate(vertices,simplices,lambda r:inverse_leading(r,signs))
                     for signs in product((-1,1),repeat=4))
        require(quadrature==expanded and quadrature>0,'Independent exact original-input integral')
        total+=quadrature
        records.append({'first':spec['first'],'second':spec['second'],
                        'radius_center':list(map(str,center)), 'mass':str(quadrature),
                        'mass_decimal':float(quadrature),'vertices':len(vertices),'simplices':len(simplices)})
    require(total==MASS,'Exact total original-input mass')
    return {'rational_reductions':identities,'routes':records,'mass':str(total),'mass_decimal':float(total)}


def raw_number_capture(B,q):
    """Full nine-slot number orbit; no borrowed digits or pair transition used."""
    value=0
    for digit in (*q,0,0,0,0,0):value=B*value+digit
    numbers={};pairs={};pair_capture=None
    while value not in numbers:
        index=len(numbers);numbers[value]=index
        n=value;digits=[]
        for _ in range(9):n,digit=divmod(n,B);digits.append(digit)
        require(n==0,'Fixed nine-slot numerical state')
        digits.sort();p=tuple(digits[-1-i]-digits[i] for i in range(4))
        if p in pairs and pair_capture is None:pair_capture=(pairs[p],index-pairs[p],p)
        elif p not in pairs:pairs[p]=index
        low=high=0
        for digit in digits:low=B*low+digit
        for digit in reversed(digits):high=B*high+digit
        value=high-low
    require(pair_capture is not None,'Pair repeat precedes numeric repeat')
    return {'actual_time':numbers[value],'period':len(numbers)-numbers[value],
            'pair_time':pair_capture[0],'pair_period':pair_capture[1],
            'first_cycle_state':pair_capture[2],'steps':len(numbers)}


def actual_checks(routes):
    checks=[];full=[];steps=0
    for index,record in enumerate(routes):
        r=tuple(map(Q,record['radius_center']))
        for signs in product((-1,1),repeat=4):
            x=inverse_leading(r,signs)
            for B in range(1000000,1000005):
                q=tuple(int(B*v) for v in x);result=exact_arrival(B,q)
                require(result['route']==(record['first'],record['second']),'Incoming route label')
                target=limit_summary(tuple(Q(v,B) for v in q))
                require(abs(Q(result['actual_time'])-B*target['time'])<=Q(168,5),'Uniform accepted clock error')
                require(max(abs(Q(v)-B*w) for v,w in zip(result['first_cycle_state'],target['terminal']))
                        <=Q(18,5),'Uniform first-capture position error')
                checks.append({'base':B,'route':result['route'],'signs':signs,
                               'source':q,'actual_time':result['actual_time'],'period':result['period'],
                               'first_cycle_state':result['first_cycle_state']})
        for j,B in enumerate((1000000+index%5,2000000+(index+2)%5)):
            signs=((1,-1,1,-1),(-1,1,-1,1))[j];x=inverse_leading(r,signs);q=tuple(int(B*v) for v in x)
            expected=exact_arrival(B,q);actual=raw_number_capture(B,q)
            for key in ('actual_time','pair_time','period','first_cycle_state'):
                require(actual[key]==expected[key],('Independent full number capture',key))
            require(actual['pair_period']==actual['period'],'Pair/numeric periods agree')
            steps+=actual['steps'];full.append({'base':B,'route':expected['route'],'source':q,**actual})
    require(len(checks)==480 and len(full)==12,'Complete finite verification inventory')
    return {'compressed_captures':checks,'independent_full_numbers':full,'independent_steps':steps}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_b_routes_certificates.json').read_text())
    geometry_result=verify_geometry(data);controls=0
    # Reject missing routes and broken reductions before any numerical clocks.
    for kind in ('missing','coefficient'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['routes'].pop()
        else:bad['routes'][0]['reductions'][0]['proof']['coefficients'][0]='-1'
        try:verify_geometry(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt route geometry accepted')
    numeric=actual_checks(geometry_result['routes'])
    result={'status':'passed',**geometry_result,**numeric,'corruptions_rejected':controls,
            'scope':'All-base conditional joint law on six first-image B-to-BBB route cells; not global exhaustion.'}
    (ROOT/'nine_b_routes_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {geometry_result["rational_reductions"]} rational polytope reductions; '
            f'{controls} corruptions rejected.\n'
            'Six route cells integrate identically by barycentric expansion and exact cubature.\n'
            f'Original-input mass {MASS} = {float(MASS):.15g}.\n'
            f'480 captures cover every inverse branch and five bases; 12 full numeric orbits\n'
            f'verify first pair and number capture over {numeric["independent_steps"]} subtraction steps.\n')
    (ROOT/'nine_b_routes_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
