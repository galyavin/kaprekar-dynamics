"""Whole-core route inventory and independent complete first arrivals.

--generate uses rational linear optimization only to choose test points.
Coverage comes from the proved exit geometry and the acyclic route tree, not
from solver acceptance, a random trajectory sample or finite base counts.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from nine_entry import CORE_FACETS,EXIT_FACETS,DELTA,evaluate,inverse_leading
from nine_b_routes import CM,L,flow,multiply,choose
from nine_b_facet import FACETS,DELTA_B
from nine_a_facet import DELTA_A
from nine_bb import DELTA_BB
from nine_core import exact_arrival,FACE_NAMES,ALLOWED
from nine_b_facet_checks import require
from nine_aaa_checks import raw_capture
from nine_aab_checks import actual_number_checks
from nine_releases_checks import exit_geometry

ROOT=Path(__file__).resolve().parent
ROWS={name:row for name,row,speed in EXIT_FACETS}
DELTAS={'A':DELTA_A,'B':DELTA_B,'BB':DELTA_BB}


def rotate(matrix,n):
    for _ in range(n):matrix=multiply(L,matrix)
    return matrix


def route_domains():
    records=[]
    def walk(domain,matrix,mode,route):
        if mode in ('AA','AB','BBB'):
            records.append({'route':(*route,mode),'domain':tuple(sorted(set(domain)))});return
        delta=DELTAS[mode];candidates=[(r,-evaluate(r[:4],delta)) for r in FACETS if evaluate(r[:4],delta)<0]
        for selected in sorted(ALLOWED[mode]):
            row=ROWS[selected];speed=-evaluate(row[:4],delta);new=list(domain)
            choose(new,matrix,row,speed,candidates);target=multiply(flow(row,speed,delta),matrix)
            if mode=='A':
                next_mode='AB' if selected=='B1' else 'AA';turn=1 if selected=='A0' else 0
            elif mode=='B':
                next_mode='A' if selected=='A0' else 'AB' if selected=='A1' else 'BB'
                turn=1 if selected=='A0' else 2 if selected in ('A1','B2') else 0
            else:next_mode='BBB' if selected=='B2' else 'AB';turn=0 if selected=='B2' else 2
            walk(new,rotate(target,turn),next_mode,(*route,mode+':'+selected))
    for first,row,speed in EXIT_FACETS:
        domain=list(CORE_FACETS);choose(domain,CM,row,speed,[(r,s) for name,r,s in EXIT_FACETS])
        matrix=multiply(flow(row,speed,DELTA),CM);mode=first[0]
        turns=(int(first[1])+1)%3 if mode=='A' else int(first[1])
        walk(domain,rotate(matrix,turns),mode,(first,))
    require(len(records)==33 and len({r['route'] for r in records})==33,'Complete finite route tree')
    return records


def generate():
    import z3
    centers=[];diagnostics=[]
    for spec in route_domains():
        r=z3.Reals('a b c d');margin=z3.Real('margin');opt=z3.Optimize();opt.set(timeout=10000)
        opt.add(margin>0,margin<=1,*[sum(g[i]*r[i] for i in range(4))+g[4]>=margin for g in spec['domain']])
        opt.maximize(margin);status=opt.check();diagnostics.append({'route':spec['route'],'status':str(status)})
        require(status==z3.sat,'A route test center was not found; investigate rather than assume emptiness.')
        model=opt.model();values=[]
        for x in r:
            v=model.eval(x);values.append(Q(v.numerator_as_long(),v.denominator_as_long()))
        require(min(evaluate(g,(*values,1)) for g in spec['domain'])>0,'Exact strict test-center verification')
        centers.append({'route':spec['route'],'radius':list(map(str,values))})
    data={'centers':centers,'discovery_only':diagnostics}
    (ROOT/'nine_core_test_points.json').write_text(json.dumps(data,indent=2)+'\n');return data


def verify_centers(data):
    specs=route_domains();require([tuple(r['route']) for r in data['centers']]==[s['route'] for s in specs],
                               'All 33 route test centers present')
    for record,spec in zip(data['centers'],specs):
        r=tuple(map(Q,record['radius']))
        require(min(evaluate(g,(*r,1)) for g in spec['domain'])>0,'Rational center lies in its strict route domain')
    # Re-enumerate exact facet intersections, including lower-dimensional
    # competitors; a missing competing facet cannot be hidden by samples.
    return exit_geometry()


def corruption_checks(data):
    rejected=0
    for kind in ('missing','point','route'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['centers'].pop()
        elif kind=='point':bad['centers'][0]['radius'][0]='-1'
        else:bad['centers'][0]['route']=['unproved']
        try:verify_centers(bad)
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Corrupt complete-route inventory accepted')
    return rejected


def actual_checks(data):
    compressed=[];complete=[];steps=0;failures=[];all_routes=set();samples={}
    for index,record in enumerate(data['centers']):
        r=tuple(map(Q,record['radius']));expected=tuple(record['route'])
        for signs in product((-1,1),repeat=4):
            x=inverse_leading(r,signs)
            for base_residue in range(5):
                # A finite accepted-base search is test setup, not a proof of
                # eventual acceptance or a claim about rejected sources.
                for multiplier in (1,2,5,10):
                    B=1000000*multiplier+base_residue;q=tuple(int(B*v) for v in x)
                    try:answer=exact_arrival(B,q)
                    except ValueError as error:
                        failures.append({'route':expected,'base':B,'reason':str(error)});continue
                    break
                else:raise ArithmeticError(('No accepted finite test source',expected,signs,base_residue))
                require(answer['route']==expected,'Complete intended leading route')
                all_routes.add(expected);compressed.append({'route':expected,'base':B,'signs':signs,
                     'pair_time':answer['pair_time'],'period':answer['period'],'first_cycle_state':answer['first_cycle_state']})
                if signs==(1,-1,1,-1) and base_residue==index%5:samples[expected]=(B,q,answer)
    require(len(all_routes)==len(samples)==33 and len(compressed)==2640,'All routes, inverse branches and base residues checked')
    for route,(B,q,answer) in samples.items():
        h,period,captured,n=raw_capture(B,q)
        require((h,period,captured)==(answer['pair_time'],answer['period'],answer['first_cycle_state']),
                'Independent complete whole-core first capture')
        steps+=n;complete.append({'route':route,'base':B,'source':q,'pair_time':h,'period':period,'first_cycle_state':captured})
    # Reuse the independent full-number reference on each of the three
    # terminal channels; its grouping key is just a test-channel identifier.
    numbers=[]
    for row in complete:numbers.append({**row,'return_length':row['route'][-1]})
    # The earlier helper requires exactly two groups; choose AA and AB, whose
    # phase-dependent terminal capture needs the most composition checks.
    numbers=[r for r in numbers if r['return_length'] in ('AA','AB')]
    numerical=actual_number_checks(numbers)
    return {'compressed_captures':compressed,'independent_complete_paths':complete,
            'independent_steps':steps,'test_setup_rejections':failures,**numerical}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_core_test_points.json').read_text())
    geometry=verify_centers(data);controls=corruption_checks(data);result=actual_checks(data)
    result.update({'status':'passed','route_count':33,'exit_geometry':geometry,'corruptions_rejected':controls,
            'covered_original_mass':'1973075435497331/24981193359360000',
            'scope':'Complete first capture on the generic one-step predecessor of the strict core, eventually on compact route interiors; this is not a phase-weighted limit law or global entry.'})
    (ROOT/'nine_core_results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    report=(f'PASS: all 33 leading routes, 16 inverse branches and five base residues: 2640 compressed first captures.\n'
            f'33 complete independent pair paths over {result["independent_steps"]} subtraction steps.\n'
            f'Two full original-number paths over {result["independent_number_steps"]} steps verify the actual clock.\n'
            f'{controls} corrupted route inventories rejected.\n'
            'Generic all-base capture now covers the existing 0.0789824331894 predecessor population.\n')
    (ROOT/'nine_core_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
