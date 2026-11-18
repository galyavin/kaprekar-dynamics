"""Solver-free BB corner proof replay and independent actual paths.

Optional --generate uses Z3 to find rational identities. Four explicit integer
cases partition the entire phase plane; no discovery-solver coverage result is
trusted by the checker or needed by the accompanying absorption argument.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from nine_bb import (B0,B1,NON_BB,DELTA_BB,THREE_ORDERS,PHASE_CYCLES,PHASES,
    decode,encode,phase_step,phase_domain,branch,twelve_orders,absorb,corner_exit,phase_entry_bound)
from nine_entry import compose,evaluate
from nine_b_facet_checks import require,from_orders,witness,rational_replay
from nine_phase_checks import replay,raw_step

ROOT=Path(__file__).resolve().parent
# Coordinates (U,Z,p,r,B,1).
CHART=((1,0,0,0,0,0),(-1,2,0,1,1,0),(1,-2,1,-1,0,0),(0,1,0,0,0,0))
EXTENDED=(*CHART,(0,0,0,0,1,0),(0,0,0,0,0,1))
MACRO_CHANGE=tuple((row[0],row[1],row[4]) for row in EXTENDED[:5])
MACRO_DOMAIN=tuple(compose(row,MACRO_CHANGE) for row in NON_BB)+((0,0,1),)


def three_data(key):
    M,guards,offsets=replay(from_orders(THREE_ORDERS[key]))
    macro=[];phase=[]
    for guard in guards:
        row=compose(guard,EXTENDED)
        if row[0] or row[1] or row[4]:macro.append((row[0],row[1],row[4]))
        else:phase.append((row[2],row[3],row[5]))
    return M,guards,offsets,macro,phase


def twelve_data(p):
    M,guards,offsets=replay(from_orders(twelve_orders(p)))
    homogeneous=sorted(set(compose(g[:5],MACRO_CHANGE) for g in guards))
    return M,guards,offsets,homogeneous


def generate():
    records=[]
    for key in range(4):
        _,_,_,macro,phase=three_data(key)
        records.append({'key':key,'macro':[witness(g,MACRO_DOMAIN) for g in macro],
                        'phase':[witness(g,(*phase_domain(key),(0,0,1))) for g in phase]})
    recurrent=[]
    for p in PHASES:
        _,_,_,hom=twelve_data(p)
        recurrent.append({'phase':p,'forward':[witness(g,MACRO_DOMAIN) for g in hom],
                          'backward':[witness(g,(*hom,(0,0,1))) for g in MACRO_DOMAIN[:-1]]})
    data={'records':records,'recurrent':recurrent}
    (ROOT/'nine_bb_certificates.json').write_text(json.dumps(data,indent=2)+'\n');return data


def verify(data):
    require([r['key'] for r in data['records']]==list(range(4)),'Four branch inventory')
    require([tuple(r['phase']) for r in data['recurrent']]==list(PHASES),'Recurrent phase inventory')
    targets=(((1,0,-3),(0,1,-3),(0,0,-1)),
             ((-1,0,-5),(-1,1,-4),(0,0,-1)),
             ((-1,0,-5),(1,-1,-6),(-1,1,0)),
             ((1,0,-3),(0,-1,-7),(0,1,1)))
    count=0;errors=[]
    for record,target in zip(data['records'],targets):
        key=record['key'];M,guards,offsets,macro,phase=three_data(key)
        count+=rational_replay(record['macro'],macro,MACRO_DOMAIN)
        count+=rational_replay(record['phase'],phase,(*phase_domain(key),(0,0,1)))
        chart=tuple(compose(r,EXTENDED) for r in M)
        pnew=tuple(chart[1][j]+chart[2][j]-(1 if j==4 else 0) for j in range(6))
        rnew=tuple(chart[0][j]+chart[1][j]-2*chart[3][j]-(1 if j==4 else 0) for j in range(6))
        dz=tuple(chart[3][j]-(1 if j==1 else 0) for j in range(6))
        require((pnew,rnew,dz)==tuple((0,0,a,b,0,c) for a,b,c in target),
                'Exact three-step phase formula')
        for j,offset in enumerate([*offsets,M]):
            require(offset[0]==(1,0,0,0,0,-j),'One unit of first-coordinate decrease per step')
    for record in data['recurrent']:
        p=tuple(record['phase']);M,guards,offsets,hom=twelve_data(p)
        count+=rational_replay(record['forward'],hom,MACRO_DOMAIN)
        count+=rational_replay(record['backward'],MACRO_DOMAIN[:-1],(*hom,(0,0,1)))
        chart=tuple(compose(r,EXTENDED) for r in M)
        for i,row in enumerate(chart):
            diff=tuple(row[j]-CHART[i][j]-(DELTA_BB[i] if j==5 else 0) for j in range(6))
            require(not(diff[0] or diff[1] or diff[4]) and diff[5]==-diff[2]*p[0]-diff[3]*p[1],
                    'Twelve-step displacement on the exact two-phase plane')
        require(tuple(evaluate(row[:4],DELTA_BB) for row in M)==DELTA_BB,'Affine repeated progression')
        for j,offset in enumerate([*offsets,M]):
            require(offset[0]==(1,0,0,0,0,-j),'Exact longitudinal corner clock')
        # q=q0+(0,r,p-r,0), with q0 on the exact two-facet intersection.
        errors.append(1+max(Q(abs(g[-1]+g[1]*p[1]+g[2]*(p[0]-p[1])),
                              -evaluate(g[:4],DELTA_BB))
                      for g in guards if evaluate(g[:4],DELTA_BB)<0))
    return count,errors


def phase_checks():
    # This finite check corroborates the elementary unbounded fold argument;
    # its grid is not used as a substitute for that argument.
    largest=0
    for initial in product(range(-64,65),repeat=2):
        p=initial;n=0
        while p not in PHASES:
            p,_=phase_step(p);n+=1
            require(n<1000,'Phase corroboration cap')
        require(n<=phase_entry_bound(initial),'Coarse unbounded phase bound')
        largest=max(largest,n)
    for cycle in PHASE_CYCLES:
        total=0
        for i,p in enumerate(cycle):
            new,dz=phase_step(p);require(new==cycle[(i+1)%4],'Exact phase cycle');total+=dz
        require(total==-7,'Common twelve-step Z displacement')
    # Finite corroboration of the two section folds derived algebraically in
    # the proof note; these samples do not replace that unbounded argument.
    for section in (-1,-2):
        probes=(-100,-6,-5,-4,-3,-2,-1,0,100)
        for r in probes:
            p,_=phase_step((section,r));p,_=phase_step(p)
            target=abs(r+2)-5 if section==-1 else abs(abs(r+3)-1)-5
            require(p==(section,target),'Exact two-return section fold')
    return {'initial_pairs_checked':129**2,'maximum_returns_observed':largest,'cycles':PHASE_CYCLES}


def actual_checks():
    paths=[];steps=0
    for B in range(200000,200005):
        for p in PHASES:
            q=encode(B,22*B//25,B//4,p);start=q;out=corner_exit(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent complete BB path')
            steps+=out['elapsed'];paths.append({'base':B,'source':start,**out})
    for B in (10000000,10000001):
        for p in ((-64,64),(64,-64),(-64,-64),(64,64),(0,0),(1,-1)):
            q=encode(B,22*B//25,B//4,p);out=absorb(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent BB phase absorption');steps+=out['elapsed']
    return paths,steps


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_bb_certificates.json').read_text())
    count,errors=verify(data);controls=0
    for kind in ('missing','coefficient','phase'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'].pop()
        elif kind=='coefficient':bad['records'][0]['macro'][0]['coefficients'][0]='-1'
        else:bad['recurrent'][0]['phase']=[0,0]
        try:verify(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt BB certificate accepted')
    phases=phase_checks();paths,steps=actual_checks()
    result={'status':'passed','rational_identities':count,'corruptions_rejected':controls,
            'repetition_error':str(max(errors)),'elapsed_error':str(12*max(errors)),
            'phase_checks':phases,'independent_paths':paths,'independent_steps':steps,
            'scope':'BB corner timed transport; general terminal arrival requires the subsequent route.'}
    (ROOT/'nine_bb_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {count} rational guard identities; {controls} corrupted proofs rejected.\n'
            'Three recurrent phase cycles share the twelve-step displacement (-12,-2,2,-7).\n'
            f'60 full corner segments and 12 absorption paths: {steps} independent subtraction steps.\n'
            f'Leading corner elapsed error at most {12*max(errors)} steps.\n'
            'Global phase absorption uses the accompanying elementary fold argument.\n')
    (ROOT/'nine_bb_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
