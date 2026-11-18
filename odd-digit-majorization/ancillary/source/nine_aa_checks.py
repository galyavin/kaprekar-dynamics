"""Unbounded AA phase contraction, guard replay and independent segments."""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from nine_aa import (A0,A1,A2,NON_AA,DELTAS,PHASES,PHASE_CYCLES,KEYS,decode,encode,
    branch,phase_domain,phase_matrix,phase_step,phase_entry_bound,three_orders,
    return_spec,absorb,corner_exit)
from nine_entry import compose,evaluate
from nine_b_facet_checks import require,from_orders,witness,rational_replay
from nine_phase_checks import replay,raw_step
from nine_b_routes import normalized

ROOT=Path(__file__).resolve().parent
CHART=((1,0,0,0,0,0),(2,0,-1,0,-1,0),(3,0,-1,-1,-2,0),(0,1,0,0,0,0))
EXTENDED=(*CHART,(0,0,0,0,1,0),(0,0,0,0,0,1))
MACRO_CHANGE=tuple((row[0],row[1],row[4]) for row in EXTENDED[:5])
MACRO_DOMAIN=tuple(normalized(compose(row,MACRO_CHANGE)) for row in NON_AA)+((0,0,1),)


def three_data(key):
    M,guards,_=replay(from_orders(three_orders(key)));macro=[];phase=[]
    for guard in guards:
        row=compose(guard,EXTENDED)
        if row[0] or row[1] or row[4]:macro.append(normalized((row[0],row[1],row[4])))
        else:phase.append(normalized((row[2],row[3],row[5])))
    return M,macro,phase


def return_data(p):
    orders,delta=return_spec(p);M,guards,_=replay(from_orders(orders))
    homogeneous=sorted(set(normalized(compose(g[:5],MACRO_CHANGE)) for g in guards))
    return M,guards,homogeneous,delta,len(orders)


def lower_rows(key):
    return tuple((a,b,c+2) for a,b,c in phase_matrix(key)[:2])


def contraction_domain(first,second):
    M=(*phase_matrix(first)[:2],(0,0,1))
    rows=(*phase_domain(first),*(compose(row,M) for row in phase_domain(second)))
    # Coordinates (p,r,M,1), with p,r>=-2 and p,r<=M, M>=5.
    domain=tuple((a,b,0,c) for a,b,c in rows)+((1,0,0,2),(0,1,0,2),
                (-1,0,1,0),(0,-1,1,0),(0,0,1,-5),(0,0,0,1))
    outputs=[compose(row,M) for row in phase_matrix(second)[:2]]
    targets=tuple((-a,-b,1,-1-c) for a,b,c in outputs)
    return domain,targets


def exit_rows(delta):
    speed=-evaluate(A0[:4],delta);rows=[]
    for g in NON_AA:
        if g==A0:continue
        slope=evaluate(g[:4],delta)
        rows.append(normalized(compose(tuple(speed*a+slope*b for a,b in zip(g,A0)),MACRO_CHANGE)))
    return rows


def generate():
    records=[]
    for key in KEYS:
        _,macro,phase=three_data(key)
        records.append({'key':key,'macro':[witness(g,MACRO_DOMAIN) for g in macro],
          'phase':[witness(g,(*phase_domain(key),(0,0,1))) for g in phase],
          'lower':[witness(g,(*phase_domain(key),(0,0,1))) for g in lower_rows(key)]})
    recurrent=[]
    for p in PHASES:
        _,_,hom,_,_=return_data(p)
        recurrent.append({'phase':p,'forward':[witness(g,MACRO_DOMAIN) for g in hom],
                          'backward':[witness(g,(*hom,(0,0,1))) for g in MACRO_DOMAIN[:-1]]})
    contraction=[]
    for first,second in product(KEYS,repeat=2):
        domain,targets=contraction_domain(first,second)
        contraction.append({'first':first,'second':second,'proofs':[witness(g,domain) for g in targets]})
    data={'records':records,'recurrent':recurrent,'contraction':contraction,
          'exit':[[witness(g,MACRO_DOMAIN) for g in exit_rows(delta)] for delta in DELTAS]}
    (ROOT/'nine_aa_certificates.json').write_text(json.dumps(data,indent=2)+'\n');return data


def verify(data):
    require([tuple(r['key']) for r in data['records']]==list(KEYS),'Nine integer-cell inventory')
    require([tuple(r['phase']) for r in data['recurrent']]==list(PHASES),'Recurrent AA inventory')
    require([(tuple(r['first']),tuple(r['second'])) for r in data['contraction']]
            ==list(product(KEYS,repeat=2)),'Complete two-step contraction inventory')
    require(len(data['exit'])==2,'Both outgoing phase classes')
    count=0;errors={9:[],42:[]}
    for record in data['records']:
        key=tuple(record['key']);M,macro,phase=three_data(key)
        count+=rational_replay(record['macro'],macro,MACRO_DOMAIN)
        count+=rational_replay(record['phase'],phase,(*phase_domain(key),(0,0,1)))
        count+=rational_replay(record['lower'],lower_rows(key),(*phase_domain(key),(0,0,1)))
        chart=tuple(compose(row,EXTENDED) for row in M)
        outputs=(tuple(2*chart[0][j]-chart[1][j]-(j==4) for j in range(6)),
                 tuple(chart[0][j]+chart[1][j]-chart[2][j]-(j==4) for j in range(6)),
                 tuple(chart[0][j]-(j==0) for j in range(6)),
                 tuple(chart[3][j]-(j==1) for j in range(6)))
        require(outputs==tuple((0,0,a,b,0,c) for a,b,c in phase_matrix(key)),'Exact three-step formulas')
    for record in data['contraction']:
        domain,targets=contraction_domain(tuple(record['first']),tuple(record['second']))
        count+=rational_replay(record['proofs'],targets,domain)
    for record in data['recurrent']:
        p=tuple(record['phase']);M,guards,hom,delta,length=return_data(p)
        count+=rational_replay(record['forward'],hom,MACRO_DOMAIN)
        count+=rational_replay(record['backward'],MACRO_DOMAIN[:-1],(*hom,(0,0,1)))
        chart=tuple(compose(row,EXTENDED) for row in M)
        for i,row in enumerate(chart):
            diff=tuple(row[j]-CHART[i][j]-(delta[i] if j==5 else 0) for j in range(6))
            require(not(diff[0] or diff[1] or diff[4]) and diff[5]==-diff[2]*p[0]-diff[3]*p[1],
                    'Return displacement on the exact phase plane')
        require(tuple(evaluate(row[:4],delta) for row in M)==delta,'Affine progression')
        error=1+max(Q(abs(g[-1]-g[1]*p[0]-g[2]*(p[0]+p[1])),
                       -evaluate(g[:4],delta)) for g in guards if evaluate(g[:4],delta)<0)
        errors[length].append(error)
    for record,delta in zip(data['exit'],DELTAS):
        rows=exit_rows(delta);require(all(any(row) for row in rows),'No forced simultaneous exit')
        count+=rational_replay(record,rows,MACRO_DOMAIN)
    return count,errors


def phase_checks():
    finite=set(product(range(-2,5),repeat=2));maximum=0
    require(all(phase_step(p)[:2] in finite for p in finite),'Invariant absorbing phase set')
    for initial in finite:
        p=initial;n=0
        while p not in PHASES:p=phase_step(p)[:2];n+=1;require(n<=9,'Finite absorption')
        maximum=max(maximum,n)
    for cycle,delta in zip(PHASE_CYCLES,DELTAS):
        du=dz=0
        for i,p in enumerate(cycle):
            pn,rn,u,z=phase_step(p);require((pn,rn)==cycle[(i+1)%len(cycle)],'Recurrent phase cycle')
            du+=u;dz+=z
        require((du,2*du,3*du,dz)==delta,'Distinct cycle velocities')
    grid=0
    for initial in product(range(-64,65),repeat=2):
        p=initial;n=0
        while p not in PHASES:
            p=phase_step(p)[:2];n+=1
            require(n<=phase_entry_bound(initial),'Unbounded contraction bound')
        grid=max(grid,n)
    return {'absorbing_states':len(finite),'maximum_finite_transient':maximum,
            'grid_states':129**2,'maximum_grid_returns':grid}


def actual_checks():
    paths=[];prefixes=[];steps=0
    for B in range(200000,200005):
        for p in PHASES:
            q=encode(B,41*B//50,21*B//100,p);start=q;out=corner_exit(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent full AA segment')
            steps+=out['elapsed'];paths.append({'base':B,'source':start,**out})
    for B in (10000000,10000001):
        for p in ((-64,64),(64,-64),(-64,-64),(64,64),(0,1),(1,-1)):
            q=encode(B,41*B//50,21*B//100,p);start=q;out=absorb(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent AA phase absorption')
            steps+=out['elapsed'];prefixes.append({'base':B,'source':start,**out})
    return paths,prefixes,steps


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_aa_certificates.json').read_text())
    count,errors=verify(data);controls=0
    for kind in ('missing','coefficient','phase','contraction'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'].pop()
        elif kind=='coefficient':bad['records'][0]['macro'][0]['coefficients'][0]='-1'
        elif kind=='phase':bad['recurrent'][0]['phase']=[0,1]
        else:bad['contraction'].pop()
        try:verify(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt AA certificate accepted')
    phases=phase_checks();paths,prefixes,steps=actual_checks()
    result={'status':'passed','rational_identities':count,'corruptions_rejected':controls,
            'repetition_errors':{k:str(max(v)) for k,v in errors.items()},
            'elapsed_errors':{k:str(k*max(v)) for k,v in errors.items()},
            'phase_checks':phases,'independent_paths':paths,'independent_prefixes':prefixes,
            'independent_steps':steps,'scope':'Two-speed AA corner transport to A0=A1=A2; no complete capture or phase-weight claim.'}
    (ROOT/'nine_aa_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {count} rational identities; {controls} corrupted proofs rejected.\n'
            'Every integer phase reaches one of two cycles, with distinct nine- and forty-two-step velocities.\n'
            f'{len(paths)} complete corner segments and {len(prefixes)} absorption paths: {steps} independent subtraction steps.\n'
            f'Leading elapsed errors: {result["elapsed_errors"]}; unique relative-interior exit is A0.\n')
    (ROOT/'nine_aa_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
