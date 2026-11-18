"""Independent finite arithmetic and actual-arrival checks of the core kernel.

The law rests on exact lattice indices, guard proofs and weighted Riemann sums.
Seeded finite-base frequencies below are diagnostics, not convergence proofs.
"""
import argparse,copy,json
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import nine_kernel as kernel
import nine_lattice as lattice
from nine_entry import C,inverse_C,inverse_leading,evaluate
from nine_core import exact_arrival
from nine_core_checks import route_domains
from nine_b_facet_checks import require,from_orders
from nine_phase_checks import replay
from odd_affine import word

ROOT=Path(__file__).resolve().parent


def integer_checks():
    rng=Random(97631);checks=0
    for _ in range(200):
        m=rng.randrange(1,5);n=rng.randrange(1,6)
        A=tuple(tuple(rng.randrange(-10,11) for j in range(n)) for i in range(m))
        U=[list(row) for row in lattice.identity(n)]
        for _ in range(10):
            if n==1:U[0][0]*=-1;continue
            i,j=rng.sample(range(n),2);a=rng.randrange(-3,4)
            for row in U:row[j]+=a*row[i]
        H=lattice.column_basis(A);changed=lattice.multiply(A,tuple(map(tuple,U)))
        require(lattice.column_basis(changed)==H and lattice.column_basis(H)==H,
                'Integer image basis is invariant under independent unimodular changes')
        checks+=1
    # Same real AAB line, different integer image: the release retains index 18.
    a=(0,-5,0,-5);c=(0,0,0,0);h=kernel.aa.A1;delta=kernel.ab.DELTA_AB
    ordinary=((2,0),(4,0),(0,1),(3,0));released=((4,),(8,),(7,),(6,))
    total=0
    for label,M,wanted in (('ordinary',ordinary,((2,),(4,),(6,),(3,))),
                           ('released',released,((36,),(72,),(108,),(54,)))):
        state=lattice.canonical(a,c,M);out=lattice.return_split(state,h,24,11,delta,5,2)
        require(all(s[2]==wanted for e,w,s in out),'Exact outgoing primitive/index-18 lattice')
        lookup={e:(w,s) for e,w,s in out};counts=Counter();b=3;B=5*b+2
        # A complete residue cube, computed by direct integer floor division.
        for z in product(range(24),repeat=len(state[2][0])):
            q=tuple(v*b+w+t for v,w,t in zip(state[0],state[1],lattice.matvec(state[2],z)))
            N=(evaluate(h,(*q,B))-11)//24+1
            target=tuple(v+N*d for v,d in zip(q,delta));e=evaluate(h,(*target,B));counts[e]+=1
            weight,s=lookup[e]
            residual=tuple(v-aa*b-cc for v,aa,cc in zip(target,s[0],s[1]))
            require(lattice.reduce_vector(residual,s[2])==(0,0,0,0),'Direct stopped point belongs to its exact image coset')
            total+=1
        require(all(Q(counts[e],sum(counts.values()))==w for e,(w,s) in lookup.items()),
                'Independent complete residue count agrees with g/s')
    # Refinement must be explicit, never silently replace B by its mean phase.
    caught=0
    bad=((1,0,0,0),(0,0,0,0),((2,),(0,),(0,),(0,)))
    try:lattice.return_split(bad,(1,0,0,0,0),2,1,(-2,0,0,0),1,0)
    except ArithmeticError as error:
        require('refinement' in str(error),'Specific progression-refinement rejection');caught+=1
    else:raise ArithmeticError('Invalid fixed base progression was accepted')
    return {'integer_basis_checks':checks,'direct_residue_points':total,'refinement_control':caught}


def first_image_checks():
    center=kernel.core.C(tuple(map(Q,json.loads((ROOT/'nine_core_test_points.json').read_text())['centers'][0]['radius'])))
    radius=inverse_C(center);vertices=kernel.faces(());count=0;words=set()
    for signs in product((-1,1),repeat=4):
        x=inverse_leading(radius,signs);B=10**8;q=tuple(int(B*v) for v in x)
        orders=word(B,q,1)['orders'];M,guards,_=replay(from_orders(orders));words.add(orders)
        lattice.unimodular_check(tuple(row[:4] for row in M))
        # On the entire strict branch, sorting the signed radii has a fixed
        # permutation; its continuous closure includes every core vertex.
        for v in vertices:
            source=inverse_leading(inverse_C(v),signs)
            require(lattice.matvec(tuple(row[:5] for row in M),(*source,1))==v,
                    'Exact leading first image on each complete inverse branch')
        for g in guards:
            values=[evaluate(g[:5],(*inverse_leading(inverse_C(v),signs),1)) for v in vertices]
            require(min(values)>=0 and (any(values) or g[5]>=0),'Original strict-branch sorting guards')
            count+=1
    require(len(words)==16,'All sixteen distinct inverse-branch words')
    return {'unimodular_first_image_branches':16,'first_image_guards':count}


def verify_table(data,recompile=True):
    specs=route_domains();centers=json.loads((ROOT/'nine_core_test_points.json').read_text())['centers']
    require(data['modulus']==5 and len(data['kernels'])==165,'Complete 33 by five kernel inventory')
    for i,(spec,record) in enumerate(zip(specs,centers)):
        require(tuple(record['route'])==spec['route'],'Expected core route order')
        radius=tuple(map(Q,record['radius']));require(min(evaluate(g,(*radius,1)) for g in spec['domain'])>0,'Exact strict route center')
        for r in range(5):
            row=data['kernels'][5*i+r]
            require(tuple(row['route'])==spec['route'] and row['base_residue']==r,'Every route and base progression retained')
            atoms={(tuple(map(Q,a['clock'])),tuple(tuple(map(Q,z)) for z in a['terminal'])):Q(a['weight']) for a in row['atoms']}
            require(len(atoms)==len(row['atoms']) and sum(atoms.values())==1 and min(atoms.values())>0,
                    'Distinct atoms and exact positive mass one')
            if recompile:require(atoms==kernel.compile_route(spec['route'],C(radius),5,r),'Saved law equals full exact guard/lattice reconstruction')
    require(Counter(len(row['atoms']) for row in data['kernels'])=={5:60,4:75,1:30},'Five, four and one observable atoms')
    return len(specs)*5


def corruptions(data):
    rejected=0
    for kind in ('missing','weight','base','clock'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['kernels'].pop()
        elif kind=='weight':bad['kernels'][0]['atoms'][0]['weight']='1/7'
        elif kind=='base':bad['kernels'][0]['base_residue']=2
        else:bad['kernels'][0]['atoms'][0]['clock'][0]='-999'
        try:verify_table(bad,recompile=kind=='clock')
        except ArithmeticError:rejected+=1
        else:raise ArithmeticError('Corrupt kernel table accepted')
    return rejected


def separation_checks(data):
    distinct=Counter()
    for i in range(33):
        rows=data['kernels'][5*i:5*i+5];kind=rows[0]['route'][-1]
        n=len({json.dumps(row['atoms'],sort_keys=True) for row in rows})
        require(n==(5 if kind=='AB' else 1),'Exact base-class separation or independence on each route')
        distinct[kind]+=1
    # Reconstruct the displayed example independently from its three drifts.
    T=kernel.I;clock=kernel.ZERO
    for name,delta,length,turns in (('A0',kernel.core.DELTA,3,1),
                                    ('B1',kernel.af.DELTA_A,12,0),
                                    ('A1',kernel.ab.DELTA_AB,24,0)):
        h=kernel.ROWS[name];speed=-evaluate(h[:4],delta)
        clock=tuple(c+Q(length,speed)*v for c,v in zip(clock,kernel.core.compose(h,T)))
        T=kernel.multiply(kernel.flow(h,speed,delta),T)
        for _ in range(turns):T=kernel.multiply(kernel.L,T)
    require(clock==(Q(7,4),-1,0,Q(-1,4),Q(-3,4)) and T[0]==(Q(-1,4),Q(1,2),0,Q(-1,4),Q(3,4)),
            'Displayed separating route clock and longitudinal coordinate')
    expected=(('11/48','5/48','137/288','55/288'),('29/128','5/48','47/96','23/128'),
              ('29/144','5/48','607/1152','193/1152'),('107/576','5/48','69/128','197/1152'),
              ('233/1152','5/48','73/144','215/1152'))
    distance=tuple(v-(Q(4,5) if i==4 else 0) for i,v in enumerate(T[0]))
    for beta,row in enumerate(data['kernels'][10:15]):
        lookup={tuple(map(Q,a['clock'])):a for a in row['atoms']}
        for k,rate in enumerate((Q(7,2),Q(9,2),Q(50,9),Q(37,6))):
            a=lookup[tuple(c+rate*v for c,v in zip(clock,distance))]
            require(Q(a['weight'])==Q(expected[beta][k]),'Exact displayed five-class probability table')
            z=tuple(tuple(map(Q,g)) for g in a['terminal'])
            if k==0:
                u=tuple(v/6+(Q(2,3) if i==4 else 0) for i,v in enumerate(T[0]))
                target=(u,tuple(2*v-(1 if i==4 else 0) for i,v in enumerate(u)),
                        tuple(3*v-(2 if i==4 else 0) for i,v in enumerate(u)),
                        tuple(-v+(1 if i==4 else 0) for i,v in enumerate(u)))
            else:target=kernel.ENDPOINT
            require(z==target,'Displayed separating terminal component')
    return {'base_independent_AA_routes':distinct['AA'],'base_independent_BBB_routes':distinct['BBB'],
            'five_distinct_AB_routes':distinct['AB'],'separating_probability_table':expected}


def actual_checks(data,samples):
    rng=Random(465120);centers=json.loads((ROOT/'nine_core_test_points.json').read_text())['centers'];out=[]
    max_error=Q(0);total=0;frequency_error=Q(0)
    for i,rec in enumerate(centers):
        x=C(tuple(map(Q,rec['radius'])))
        for r in range(5):
            row=data['kernels'][5*i+r];atoms=[(tuple(map(Q,a['clock'])),tuple(tuple(map(Q,z)) for z in a['terminal'])) for a in row['atoms']]
            B=10**9+r;counts=[0]*len(atoms)
            for _ in range(samples):
                # A large integer box around a strict rational interior point;
                # no kernel probabilities enter source selection.
                q=tuple(int(B*v)+rng.randrange(-100000,100001) for v in x)
                answer=exact_arrival(B,q,first_image=False)
                require(tuple(answer['route'])==tuple(rec['route']),'Independent actual source retains its intended route')
                xn=tuple(Q(v,B) for v in q);observed=(Q(answer['actual_time'],B),*(Q(v,B) for v in answer['first_cycle_state']))
                distances=[]
                for clock,z in atoms:
                    predicted=(evaluate(clock,(*xn,1)),*(evaluate(g,(*xn,1)) for g in z))
                    distances.append(max(abs(a-b) for a,b in zip(observed,predicted)))
                matches=[j for j,d in enumerate(distances) if d<Q(1,100000)]
                require(len(matches)==1,('Unique independent actual joint atom',i,r,distances))
                j=matches[0];counts[j]+=1;max_error=max(max_error,distances[j]);total+=1
            err=max(abs(Q(n,samples)-Q(a['weight'])) for n,a in zip(counts,row['atoms']))
            frequency_error=max(frequency_error,err)
            out.append({'route_index':i,'base_residue':r,'counts':counts,'maximum_frequency_error':str(err)})
        print('Actual route',i,'passed',flush=True)
    # Diagnostic tolerance only. Every mathematical weight is derived by the
    # exact index proof; no inference of convergence uses this finite sample.
    require(frequency_error<Q(1,4),'Gross finite-frequency diagnostic disagreement')
    return {'actual_captures':total,'maximum_joint_error':str(max_error),
            'maximum_frequency_error':str(frequency_error),'frequency_diagnostics':out}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--samples',type=int,default=96);args=parser.parse_args()
    require(args.samples>=64,'Use at least 64 points per route/base diagnostic')
    data=json.loads((ROOT/'nine_kernel_results.json').read_text());result={**integer_checks(),**first_image_checks()}
    result['reconstructed_kernels']=verify_table(data);result['corruptions_rejected']=corruptions(data)
    result.update(separation_checks(data))
    result.update(actual_checks(data,args.samples));result['status']='passed';result['audit']=dict(kernel.STATS)
    (ROOT/'nine_kernel_check_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f"PASS: 165 exact weighted kernels on all 33 core routes and five base classes.\n"
            f"{result['integer_basis_checks']} integer basis changes; {result['direct_residue_points']} complete residue-grid points.\n"
            f"16 unimodular first-image branches; four corrupt kernel tables rejected.\n"
            f"{result['actual_captures']} independent actual compressed first captures agree with a unique joint atom.\n"
            f"Maximum joint error {result['maximum_joint_error']}; frequencies are diagnostics only.\n")
    (ROOT/'nine_kernel_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
