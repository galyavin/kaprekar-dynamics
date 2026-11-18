"""Finite actual-prefix controls after global entry."""
from fractions import Fraction as Q
from itertools import product, combinations
import sympy as sy
from nine_entry import inverse_leading,evaluate
from nine_aaa_checks import raw_capture
from nine_phase_checks import raw_step
from odd_affine import word


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def actual_checks(points,kernel,arrival,factors=(1,2,5,10)):
    rows=[];rejected=[];maximum_error=Q(0);prefix_words=set();samples=[]
    for i,r in enumerate(points):
        for signs in product((-1,1),repeat=4):
            x=inverse_leading(r,signs)
            for beta in range(5):
                for factor in factors:
                    B=1000000*factor+beta;q=tuple(int(B*v) for v in x)
                    try:
                        expected=kernel(tuple(Q(v,B) for v in q),beta)
                        answer=arrival(B,q)
                        if answer['route']!=expected['route']:
                            raise ValueError(f'Finite-base route differs: actual {answer["route"]}, leading {expected["route"]}')
                        prefix=answer['prefix_steps'];certificate=word(B,q,prefix);M=certificate['matrix']
                        require(tuple(evaluate(g,(*q,B,1)) for g in M)==answer['core_source'],'Complete affine prefix matches actual subtraction')
                        # At a finite base, a thin borrowed sorting cell can
                        # still differ from the source's strict leading word.
                        # Unimodularity is asserted only after that hypothesis
                        # is checked; such setups require a larger base.
                        if any(evaluate(g[:5],(*[Q(v,B) for v in q],1))<=0
                               for g in certificate['guards']):
                            raise ValueError('Finite-base prefix violates a strict leading guard')
                        require(abs(sy.Matrix([g[:4] for g in M]).det())==1,'Complete incoming lattice stays full')
                        if tuple(evaluate(g[:5],(*[Q(v,B) for v in q],1)) for g in M)!=expected['core_source']:
                            raise ValueError('Finite-base prefix does not follow the strict leading word')
                        terminal = answer['first_cyclic_number_differences']
                        require(terminal == raw_step(B, answer['first_cycle_state']),
                                'First cyclic number difference phase from raw subtraction')
                        observed=(Q(answer['actual_time'],B),*(Q(v,B) for v in terminal))
                        atoms=sorted(set((atom['time'],*atom['terminal']) for atom in expected['atoms']))
                        separation=min((max(abs(a-b) for a,b in zip(u,v)) for u,v in combinations(atoms,2)),default=Q(1))
                        radius=min(Q(1,500),separation/3)
                        errors=[max(abs(a-b) for a,b in zip(observed,atom)) for atom in atoms]
                        matches=[v for v in errors if v<radius]
                        if len(matches)!=1:
                            raise ValueError(f'Finite-base joint error {min(errors)} exceeds separated-atom radius {radius}')
                    except ValueError as error:
                        rejected.append({'center':i,'base':B,'signs':signs,'reason':str(error)})
                        continue
                    break
                else:raise ArithmeticError(('No accepted finite source',i,signs,beta,rejected[-len(factors):]))
                prefix_words.add(certificate['orders'])
                maximum_error=max(maximum_error,matches[0])
                rows.append({'center':i,'base':B,'signs':signs,'prefix_steps':prefix,'route':answer['route'],
                             'pair_time':answer['pair_time'],'period':answer['period'],'first_cycle_state':answer['first_cycle_state'],
                             'atom_comparison_radius':str(radius)})
                if i<2 and signs==(1,-1,1,-1) and beta==i:samples.append((B,q,answer))
        print('Actual center',i,'all branches and base classes passed',flush=True)
    complete=[];raw_steps=0
    for B,q,answer in samples:
        h,period,first,n=raw_capture(B,q)
        require((h,period,first)==(answer['pair_time'],answer['period'],answer['first_cycle_state']),
                'Independent complete first-capture path agrees')
        raw_steps+=n;complete.append({'base':B,'source':q,'pair_time':h,'period':period,'first_cycle_state':first})
        print('Independent full first capture',B,h,flush=True)
    return {'compressed_captures':rows,'finite_setup_rejections':rejected,'distinct_actual_prefix_words':len(prefix_words),
            'maximum_joint_error':str(maximum_error),'independent_complete_captures':complete,'independent_pair_steps':raw_steps}
