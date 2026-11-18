"""Solver-free A-B corner proof replay and independent numeric segments."""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from nine_ab import (A1,A2,B1,NON_AB,DELTA_AB,PHASES,PHASE_CYCLES,KEYS,
    decode,encode,phase_step,phase_entry_bound,branch,phase_domain,three_orders,
    return_orders,absorb,corner_exit)
from nine_entry import compose,evaluate
from nine_b_facet_checks import require,from_orders,witness,rational_replay
from nine_phase_checks import replay,raw_step
from nine_b_routes import normalized

ROOT=Path(__file__).resolve().parent
CHART=((1,0,0,0,0,0),(2,0,-1,0,-1,0),(0,1,0,0,0,0),
       (Q(3,2),0,Q(-1,2),Q(-1,2),-1,0))
EXTENDED=(*CHART,(0,0,0,0,1,0),(0,0,0,0,0,1))
MACRO_CHANGE=tuple((row[0],row[1],row[4]) for row in EXTENDED[:5])
MACRO_DOMAIN=tuple(normalized(compose(row,MACRO_CHANGE)) for row in NON_AB)+((0,0,1),)


def three_data(key):
    M,guards,offsets=replay(from_orders(three_orders(key)))
    macro=[];phase=[]
    for guard in guards:
        row=compose(guard,EXTENDED)
        if row[0] or row[1] or row[4]:macro.append(normalized((row[0],row[1],row[4])))
        else:phase.append(normalized((row[2],row[3],row[5])))
    return M,guards,macro,phase


def return_data(p):
    M,guards,offsets=replay(from_orders(return_orders(p)))
    homogeneous=sorted(set(normalized(compose(g[:5],MACRO_CHANGE)) for g in guards))
    return M,guards,homogeneous


def exit_rows():
    # Every other facet is strictly positive at the proposed A1 exit.
    rows=[]
    for g in NON_AB:
        if g==A1:continue
        slope=evaluate(g[:4],DELTA_AB)
        rows.append(normalized(compose(tuple(24*a+slope*b for a,b in zip(g,A1)),MACRO_CHANGE)))
    return rows


def generate():
    records=[]
    for key in KEYS:
        _,_,macro,phase=three_data(key)
        records.append({'key':key,'macro':[witness(g,MACRO_DOMAIN) for g in macro],
                        'phase':[witness(g,(*phase_domain(key),(0,0,1))) for g in phase]})
    recurrent=[]
    for p in PHASES:
        _,_,hom=return_data(p)
        recurrent.append({'phase':p,'forward':[witness(g,MACRO_DOMAIN) for g in hom],
                          'backward':[witness(g,(*hom,(0,0,1))) for g in MACRO_DOMAIN[:-1]]})
    data={'records':records,'recurrent':recurrent,
          'exit':[witness(g,MACRO_DOMAIN) for g in exit_rows()]}
    (ROOT/'nine_ab_certificates.json').write_text(json.dumps(data,indent=2)+'\n');return data


def verify(data):
    require([tuple(r['key']) for r in data['records']]==list(KEYS),'Six integer-cell inventory')
    require([tuple(r['phase']) for r in data['recurrent']]==list(PHASES),'Recurrent AB inventory')
    count=0;errors=[]
    for record in data['records']:
        key=tuple(record['key']);M,guards,macro,phase=three_data(key)
        count+=rational_replay(record['macro'],macro,MACRO_DOMAIN)
        count+=rational_replay(record['phase'],phase,(*phase_domain(key),(0,0,1)))
        chart=tuple(compose(row,EXTENDED) for row in M)
        outputs=(tuple(2*chart[0][j]-chart[1][j]-(j==4) for j in range(6)),
                 tuple(chart[0][j]+chart[1][j]-2*chart[3][j]-(j==4) for j in range(6)),
                 tuple(chart[0][j]-(j==0) for j in range(6)),
                 tuple(chart[2][j]-(j==1) for j in range(6)))
        a,b=key
        p={'low':(-1,0,2),'middle':(-2,0,2),'high':(1,0,-4)}[a]
        r=(-1 if a=='low' else 0,1 if b=='high' else -1,-3 if b=='high' else -7)
        du=(-1,0,-1) if a!='high' else (0,0,-3)
        dw=(1 if a=='low' else 0,0,-1)
        require(outputs==tuple((0,0,x,y,0,z) for x,y,z in (p,r,du,dw)),'Exact three-step formulas')
    for record in data['recurrent']:
        p=tuple(record['phase']);M,guards,hom=return_data(p)
        count+=rational_replay(record['forward'],hom,MACRO_DOMAIN)
        count+=rational_replay(record['backward'],MACRO_DOMAIN[:-1],(*hom,(0,0,1)))
        chart=tuple(compose(row,EXTENDED) for row in M)
        for i,row in enumerate(chart):
            diff=tuple(row[j]-CHART[i][j]-(DELTA_AB[i] if j==5 else 0) for j in range(6))
            require(not(diff[0] or diff[1] or diff[4]) and diff[5]==-diff[2]*p[0]-diff[3]*p[1],
                    'Twenty-four-step displacement on the exact phase plane')
        require(tuple(evaluate(row[:4],DELTA_AB) for row in M)==DELTA_AB,'Affine progression')
        errors.append(1+max(Q(abs(Q(g[-1])-g[1]*p[0]-Q(g[3]*(p[0]+p[1]),2)),
                              -evaluate(g[:4],DELTA_AB))
                      for g in guards if evaluate(g[:4],DELTA_AB)<0))
    require(all(any(row) for row in exit_rows()),'No forced simultaneous exit')
    count+=rational_replay(data['exit'],exit_rows(),MACRO_DOMAIN)
    return count,errors


def phase_checks():
    finite=set(product((-2,-1,0,2,3,4),range(-5,1)));maximum=0
    require(all(phase_step(p)[0] in finite for p in finite),'Invariant absorbing phase set')
    for initial in finite:
        p=initial;n=0
        while p not in PHASES:p,_,_=phase_step(p);n+=1;require(n<=36,'Finite absorption')
        maximum=max(maximum,n)
    for cycle in PHASE_CYCLES:
        du=dw=0
        for i,p in enumerate(cycle):
            new,u,w=phase_step(p);require(new==cycle[(i+1)%len(cycle)],'Recurrent phase cycle')
            du+=u;dw+=w
        require(du==dw==-3*len(cycle)//2,'Common average drift')
    grid=0
    for initial in product(range(-64,65),repeat=2):
        p=initial;n=0
        while p not in PHASES:
            p,_,_=phase_step(p);n+=1
            require(n<=phase_entry_bound(initial),'Coarse unbounded phase bound')
        grid=max(grid,n)
    return {'absorbing_states':len(finite),'maximum_finite_transient':maximum,
            'grid_states':129**2,'maximum_grid_returns':grid}


def actual_checks():
    paths=[];prefixes=[];steps=0
    for B in range(200000,200005):
        for p in PHASES:
            U=21*B//25;U+=(sum(p)-U)%2;q=encode(B,U,49*B//100,p);start=q
            out=corner_exit(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent full AB segment')
            steps+=out['elapsed'];paths.append({'base':B,'source':start,**out})
    for B in (10000000,10000001):
        for p in ((-64,64),(64,-64),(-64,-64),(64,64),(0,0),(1,-1)):
            U=21*B//25;U+=(sum(p)-U)%2;q=encode(B,U,49*B//100,p);start=q;out=absorb(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent AB phase absorption')
            steps+=out['elapsed'];prefixes.append({'base':B,'source':start,**out})
    return paths,prefixes,steps


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_ab_certificates.json').read_text())
    count,errors=verify(data);controls=0
    for kind in ('missing','coefficient','phase'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'].pop()
        elif kind=='coefficient':bad['records'][0]['macro'][0]['coefficients'][0]='-1'
        else:bad['recurrent'][0]['phase']=[0,0]
        try:verify(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt AB certificate accepted')
    phases=phase_checks();paths,prefixes,steps=actual_checks()
    result={'status':'passed','rational_identities':count,'corruptions_rejected':controls,
            'repetition_error':str(max(errors)),'elapsed_error':str(24*max(errors)),
            'phase_checks':phases,'independent_paths':paths,'independent_prefixes':prefixes,
            'independent_steps':steps,'scope':'AB corner timed transport to A1=A2=B1; no terminal-arrival claim.'}
    (ROOT/'nine_ab_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {count} rational identities; {controls} corrupted proofs rejected.\n'
            'Every integer phase absorbs; three cycles share the twenty-four-step drift (-12,-24,-12,-18).\n'
            f'{len(paths)} complete corner segments and {len(prefixes)} absorption paths: {steps} independent subtraction steps.\n'
            f'Leading elapsed error at most {24*max(errors)} steps; unique relative-interior exit is A1.\n')
    (ROOT/'nine_ab_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
