"""Exact A-facet guard replay, finite phase absorption and independent paths.

Default verification uses the standard library. Optional --generate uses Z3
only to find coefficients that are subsequently replayed as rational identities.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

from nine_entry import evaluate,compose,bulk_exit,inverse_leading
from nine_a_facet import (A2,NON_A,PHASES,DELTA_A,THREE_ORDERS,phase,phase_next,
    three_formula,twelve_orders,guard_domain,absorb,facet_exit,leading_exit)
from nine_b_facet_checks import IDENTITY,require,from_orders,witness,rational_replay
from nine_phase_checks import replay,raw_step

ROOT=Path(__file__).resolve().parent


def domains(e):
    record=from_orders(twelve_orders(e));M,guards,_=replay(record)
    homogeneous=sorted(set(g[:5] for g in guards))
    equalities=(A2,tuple(-x for x in A2),(0,0,0,0,1))
    return M,guards,homogeneous,(*NON_A,*equalities),(*homogeneous,*equalities)


def generate():
    records=[]
    for key,orders in THREE_ORDERS.items():
        _,guards,_=replay(from_orders(orders))
        records.append({'key':key,'proofs':[witness(g,guard_domain(key)) for g in guards]})
    recurrent=[]
    for e in PHASES:
        _,_,hom,forward,backward=domains(e)
        recurrent.append({'phase':e,'forward':[witness(g,forward) for g in hom],
                          'backward':[witness(g,backward) for g in NON_A]})
    data={'records':records,'recurrent':recurrent}
    (ROOT/'nine_a_facet_certificates.json').write_text(json.dumps(data,indent=2)+'\n')
    return data


def verify(data):
    require([r['key'] for r in data['records']]==list(THREE_ORDERS),'Three-word inventory')
    require([r['phase'] for r in data['recurrent']]==list(PHASES),'Recurrent-word inventory')
    expected=(
        ((-1,1,0,0,1,-1),(0,1,0,0,0,-4),(2,-1,1,0,-1,-1),(0,0,0,1,0,-1)),
        ((-1,1,0,0,1,-1),(2,0,0,0,-1,-4),(0,0,1,0,0,-1),(0,0,0,1,0,-1)),
        ((1,0,0,0,0,-3),(0,1,0,0,0,-2),(0,0,1,0,0,-1),(0,0,0,1,0,-1)))
    count=0;errors=[]
    for record,target in zip(data['records'],expected):
        key=record['key'];M,guards,_=replay(from_orders(THREE_ORDERS[key]))
        require(M==target,'Three-step exact map')
        count+=rational_replay(record['proofs'],guards,guard_domain(key))
    for record in data['recurrent']:
        e=record['phase'];M,guards,hom,forward,backward=domains(e)
        count+=rational_replay(record['forward'],hom,forward)
        count+=rational_replay(record['backward'],NON_A,backward)
        for i,row in enumerate(M):
            diff=tuple(row[j]-IDENTITY[i][j]-(DELTA_A[i] if j==5 else 0) for j in range(6))
            multiplier=-diff[1]
            require(diff==tuple(multiplier*x for x in (*A2,-e)),
                    'Twelve-step displacement on the exact phase plane')
        require(tuple(evaluate(row[:4],DELTA_A) for row in M)==DELTA_A,'Affine iteration')
        f=e
        for _ in range(4):f=phase_next(f)
        require(f==e,'Closed recurrent phase')
        # q0=(u,v+e,w,z) projects onto A2. Constants are c-e*g_v.
        errors.append(1+max(Q(abs(g[-1]-e*g[1]),-evaluate(g[:4],DELTA_A))
                      for g in guards if evaluate(g[:4],DELTA_A)<0))
    return count,errors


def absorption_checks():
    max_steps=max_error=max_loss=0
    for e in range(-512,513):
        B=10000000;start=(8200000,6400000-e,4300000,2100000);q=start;elapsed=0
        while phase(B,q) not in PHASES:
            q=three_formula(B,q);elapsed+=3
            require(-512<=phase(B,q)<=514,'Bounded intermediate phase')
            d=tuple(a-b for a,b in zip(q,start))
            max_error=max([max_error,*map(abs,d)])
            max_loss=max([max_loss,*[-evaluate(r[:4],d) for r in NON_A]])
            require(elapsed<=390,'Finite absorption bound')
        require(absorb(B,start)['target']==q,'Runtime absorption')
        max_steps=max(max_steps,elapsed)
    require(max_loss<=4000,'14000 margin retains the 10000 guard buffer')
    return {'initial_phases':1025,'max_steps':max_steps,
            'max_coordinate_error':max_error,'max_facet_loss':max_loss}


def actual_checks():
    records=[];steps=0
    for B in range(200000,200005):
        for e in PHASES:
            u=41*B//50;q=(u,2*u-B-e,43*B//100,21*B//100);start=q
            result=facet_exit(B,q)
            for _ in range(result['elapsed']):q=raw_step(B,q)
            require(q==result['target'],'Independent complete A-facet path')
            steps+=result['elapsed'];records.append({'base':B,'start':start,**result})
    for B in range(10000000,10000005):
        for e in (-512,-511,-3,-2,-1,0,1,2,3,4,511,512):
            u=41*B//50;q=(u,2*u-B-e,43*B//100,21*B//100);result=absorb(B,q)
            for _ in range(result['elapsed']):q=raw_step(B,q)
            require(q==result['target'],'Independent A-facet absorption')
            steps+=result['elapsed']
    return records,steps


def transfer_checks():
    data=json.loads((ROOT/'nine_entry_certificates.json').read_text());records=[]
    # A2 L=A0, A2 L^2=A1, A2 L^0=A2.
    for piece in data['exit_partitions']:
        if not piece['face'].startswith('A'):continue
        rotation=(int(piece['face'][1])+1)%3;r=tuple(map(Q,piece['radius_center']))
        for signs in product((-1,1),repeat=4):
            x=inverse_leading(r,signs)
            for B in range(100000000,100000005):
                source=tuple(int(B*t) for t in x);q=raw_step(B,source)
                bulk=bulk_exit(B,q);q=bulk['target']
                for _ in range(rotation):q=raw_step(B,q)
                require(abs(phase(B,q))<=318,'Bounded incoming A offset')
                absorption=absorb(B,q);q=absorption['target'];facet=facet_exit(B,q)
                require(facet['repetitions']>0,'Positive A-facet access')
                records.append({'base':B,'face':piece['face'],'signs':signs,
                                'absorption_steps':absorption['elapsed'],
                                'phase':phase(B,q),'facet_steps':facet['elapsed']})
    require(len(records)==240,'All A faces, inverse branches and five bases')
    return records


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_a_facet_certificates.json').read_text())
    count,errors=verify(data);controls=0
    for kind in ('missing','coefficient','phase'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['records'][0]['proofs'].pop()
        elif kind=='coefficient':bad['records'][0]['proofs'][0]['coefficients'][0]='-1'
        else:bad['recurrent'][0]['phase']=1
        try:verify(bad)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupt A-facet proof accepted')
    absorption=absorption_checks();actual,steps=actual_checks();transfer=transfer_checks()
    error=max(errors)
    result={'status':'passed','rational_identities':count,'corruptions_rejected':controls,
            'absorption':absorption,'repetition_errors':dict(zip(map(str,PHASES),map(str,errors))),
            'elapsed_error':str(12*error),'endpoint_linf_error':str(4+12*error),
            'independent_steps':steps,'full_facet_segments':actual,'transfer_checks':transfer,
            'scope':'First A-facet timed paths; no terminal-arrival or global law claim.'}
    (ROOT/'nine_a_facet_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {count} rational identities; {controls} corruptions rejected.\n'
            f'All 1025 incoming integer phases absorbed in at most {absorption["max_steps"]} steps.\n'
            f'30 full facet paths and 60 absorption paths: {steps} independent subtraction steps.\n'
            f'Leading facet elapsed error at most {12*error} steps.\n'
            '240 transfers across three A faces, sixteen inverse branches and five bases.\n')
    (ROOT/'nine_a_facet_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
