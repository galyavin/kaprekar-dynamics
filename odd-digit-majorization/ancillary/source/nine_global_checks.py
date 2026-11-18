"""Independent leading, local-lemma and actual joint-transfer controls."""
import hashlib,json
from fractions import Fraction as Q
from itertools import combinations_with_replacement,product
from pathlib import Path
from random import Random
import sympy as sy
from nine_entry import C,H,CORE_FACETS,evaluate,in_core_radius,inverse_leading
from nine_leading_checks import reference,known
from nine_shear_checks import actual_checks
from nine_global import leading_target,limit_kernel,exact_arrival

ROOT=Path(__file__).resolve().parent


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def dg(q):
    a,b,c,d=q
    return (1-a,1-b,a-c,b-d,c+d)


def nu(q):
    a,b,c,d=q
    return b-c+abs(c-2*a+1)


def violation(q):
    return max(Q(0),*[-evaluate(g,(*q,1)) for g in CORE_FACETS])


def ordered(q):
    return 1>q[0]>=q[1]>=q[2]>=q[3]>=0


def allowed(q):
    return in_core_radius(q) or known(q,2) or known(q,3)


def orbit(q,steps):
    out=[q]
    for _ in range(steps):out.append(reference(out[-1]))
    return out


def leading_checks():
    states=steps=zero_blocks=cycles=0;periods=set();exceptions=[]
    for asc in combinations_with_replacement(range(18),4):
        q=tuple(Q(v,18) for v in reversed(asc));path=orbit(q,4)
        require(all(H(a)==b for a,b in zip(path,path[1:])),
                'Independent nine-value map equals radius conjugacy')
        require(all(min(dg(b))>=min(dg(a)) for a,b in zip(path,path[1:])),
                'Leading minimum-gap barrier')
        if all(p[0]==q[0] and sorted(dg(p))==sorted(dg(q)) for p in path):
            require(allowed(path[-1]),'Four-step zero-loss classification');zero_blocks+=1
        states+=1;steps+=4
        if min(dg(q))<=0:continue
        seen={};r=q
        for t in range(1000):
            if in_core_radius(r):break
            if r in seen:
                period=t-seen[r];periods.add(period);cycles+=1
                require(period in (1,2,3) and allowed(r),'Noncore rational cycle classification')
                require(any(s[1]==s[2]==2*s[0]-1 for s in seen),
                        'Positive-gap nonentry has a finite exceptional-plane antecedent')
                if len(exceptions)<6:exceptions.append(list(map(str,q)))
                break
            seen[r]=t;r=reference(r);steps+=1
        else:raise ArithmeticError('Finite grid control exceeded its declared resource cap')
    return {'grid_states':states,'independent_grid_steps':steps,'zero_loss_blocks':zero_blocks,
            'positive_gap_noncore_grid_orbits':cycles,'noncore_grid_periods':sorted(periods),
            'exceptional_grid_sources':exceptions}


def local_checks():
    stabilizations=persistences=steps=0
    for eta in (Q(1,10000),Q(1,100000),Q(1,1000000)):
        for a in (Q(5,7)-eta,Q(5,7),Q(5,7)+eta,Q(11,15),Q(3,4)-eta,Q(3,4)):
            for d,p,s in product((Q(1,10),1-a-eta,1-a,1-a+eta),
                                 (Q(0),eta/2,eta),(-eta,-eta/2,Q(0),eta/2,eta)):
                q=(a,2*a-1+s+p,2*a-1+s,d)
                if not ordered(q) or d<100*eta:continue
                path=orbit(q,6);steps+=6;stabilizations+=1
                require(path[6][0]==path[3][0],'Exceptional first-coordinate stabilization')
                if a>=Q(5,7) and path[6][0]==a:
                    require(nu(path[6])>=nu(path[3]),'Prepared transverse persistence');persistences+=1
    vertices=[tuple(map(Q,p)) for p in json.loads((ROOT/'nine_entry_certificates.json').read_text())['vertices']]
    centers=[tuple((u+v)/2 for u,v in zip(a,b)) for a,b in combinations_with_replacement(vertices,2)]
    rng=Random(4800917);sources=[]
    for eps in (Q(1,100000),Q(1,1000000)):
        for center in centers:
            for _ in range(12):
                q=tuple(v+eps*rng.randrange(-3,4) for v in center);sources.append((q,20*eps))
        for a,d,p,s in product((Q(5,7),Q(11,15),Q(3,4)),
                              (20*eps,100*eps),(eps,2*eps),(-eps,eps/2)):
            sources.append(((a,2*a-1+s+p,2*a-1+s,d),500*eps))
    near=outside=absorbed=0
    for q,eta in sources:
        if not ordered(q) or min(dg(q))<100*eta or violation(q)>eta:continue
        path=orbit(q,9);steps+=9;near+=1
        if in_core_radius(path[9]):absorbed+=1;continue
        require(tuple(a-b for a,b in zip(path[9],path[6]))==
                tuple(a-b for a,b in zip(path[6],path[3])),'Prepared affine progression near the core')
        outside+=1
    require(stabilizations>500 and persistences>100 and near>300 and outside>10,
            'Local controls include nontrivial nonabsorbed cases')
    return {'exceptional_stabilization_sources':stabilizations,'transverse_persistence_sources':persistences,
            'near_core_sources':near,'near_core_absorbed_sources':absorbed,
            'near_core_nonabsorbed_progressions':outside,'independent_local_steps':steps}


def centers():
    rng=Random(4816709);out=[]
    # The first two and the final center lie below the old high-radius extensions.
    windows=[(Q(2,5),Q(1,2)),(Q(61,100),Q(69,100)),(Q(9,10),Q(97,100))]
    for lower,upper in windows:
        for _ in range(200):
            a=lower+(upper-lower)*Q(rng.randrange(1,100000),100001)
            tail=sorted((a*Q(rng.randrange(1,100000),100001) for _ in range(3)),reverse=True)
            r=(a,*tail)
            if in_core_radius(r):continue
            try:limit_kernel(inverse_leading(r,(1,-1,1,-1)),0)
            except ValueError:continue
            out.append(r);break
        else:raise ArithmeticError('No generic actual-test center')
    a=Q(11,15);s=-Q(1,7000);p=Q(1,3000)
    r=(a,2*a-1+s+p,2*a-1+s,Q(1001,10000))
    limit_kernel(inverse_leading(r,(1,-1,1,-1)),0);out.append(r)
    return out


def irrational_checks():
    def raw(r):
        half=sy.Rational(1,2)
        q=(half*(1+r[0]),*[half*(a+b) for a,b in zip(r,r[1:])])
        digits=sorted((sy.Integer(1),*q,*[1-v for v in q]),reverse=True)
        pair=[sy.expand(digits[i]-digits[-1-i]) for i in range(4)]
        out=[];previous=sy.Integer(1)
        for v in pair:
            previous=sy.expand(2*v-previous);out.append(previous)
        return tuple(out)
    def in_closed_core(r):
        return all(bool(sy.expand(sum(c*v for c,v in zip(g,(*r,1))))>=0) for g in CORE_FACETS)
    u=sy.sqrt(2);a=sy.Rational(11,15)
    sources=[(sy.Rational(43,100)+u/10000,sy.Rational(31,100),sy.Rational(17,100),sy.Rational(7,100)),
             (a,2*a-1+u/2000,2*a-1-u/2000,sy.Rational(1,10))]
    rows=[];steps=0
    for r in sources:
        source=r
        for t in range(2000):
            if in_closed_core(r):break
            r=raw(r);steps+=1
        else:raise ArithmeticError('Algebraic control reached its resource cap')
        require(raw(raw(raw(r)))==r,'Exact irrational core periodicity')
        rows.append({'source':list(map(str,source)),'entry_steps':t,'target':list(map(str,r))})
    a=sy.Rational(5,7)+u/100;d=sy.Rational(1,10)+u/1000
    r=(a,2*a-1,2*a-1,d)
    require(not in_closed_core(r) and raw(raw(raw(r)))==r,'Exact irrational exceptional cycle')
    # Without integer slopes the final scalar argument is false: t -> t/2
    # converges to zero without finite arrival from every positive t.
    return {'irrational_entry_controls':rows,'independent_irrational_entry_steps':steps,
            'irrational_exceptional_cycle':list(map(str,r)),
            'integer_slope_scope_control':'The scalar map t/2 has asymptotic but no finite zero entry.'}


def prefix_scope_control(points,rejections):
    from odd_affine import word
    for row in rejections:
        if row['reason']!='Finite-base prefix violates a strict leading guard':continue
        B=row['base'];x=inverse_leading(points[row['center']],row['signs'])
        q=tuple(int(B*v) for v in x);_,n=leading_target(tuple(Q(v,B) for v in q))
        cert=word(B,q,n);det=int(sy.Matrix([r[:4] for r in cert['matrix']]).det())
        if abs(det)==1:continue
        margin=min(evaluate(g[:5],(*[Q(v,B) for v in q],1)) for g in cert['guards'])
        require(margin<=0,'Non-unimodular actual prefix fails the strict leading hypothesis')
        return {'base':B,'source':q,'steps':n,'actual_branch_determinant':det,
                'minimum_leading_guard':str(margin),
                'scope':'A finite borrowed prefix need not be unimodular before its strict leading guards hold.'}
    raise ArithmeticError('Missing retained counterexample to unguarded actual-prefix unimodularity')


def main():
    audit=json.loads((ROOT/'nine_global_entry_results.json').read_text())
    statements=json.loads((ROOT/'nine_global_entry_obligations.json').read_text())
    require(len(audit['records'])==4 and audit['actual_solver_calls']==12,'Complete accepted audit')
    for r,s in zip(audit['records'],statements):
        require(r['key']==s['key'] and r['z3']==r['cvc5']=='unsat','Both engines accepted each lemma')
        require(r['sha256']==hashlib.sha256(s['negated_lemma_smtlib'].encode()).hexdigest(),'Exact statement hash')
    result=leading_checks();print('Leading controls',result,flush=True)
    result.update(local_checks());print('Local controls passed',flush=True)
    result.update(irrational_checks());print('Exact irrational controls passed',flush=True)
    points=centers();result['radius_centers']=[[str(v) for v in r] for r in points]
    result['leading_prefix_lengths']=[leading_target(inverse_leading(r,(1,-1,1,-1)))[1] for r in points]
    factors=(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,
             100000,200000,500000,1000000)
    result.update(actual_checks(points,kernel=limit_kernel,arrival=exact_arrival,factors=factors))
    result['strict_prefix_hypothesis_control']=prefix_scope_control(points,result['finite_setup_rejections'])
    h=Q(1,8);s=1-2*h
    complement=s**9+18*h*s**8+72*h*h*s**7
    u,v=sy.symbols('u v');density=3024*u*v*(1-u-v)**5
    integral=sy.integrate(sy.integrate(density,(v,sy.Rational(1,8),1-u)),
                          (u,sy.Rational(1,8),sy.Rational(7,8)))
    require(Q(integral)==complement,'Independent two-gap integral of the original high-radius mass')
    require(sy.integrate(sy.integrate(density,(v,0,1-u)),(u,0,1))==1,'Original marginal density normalization')
    result['high_initial_radius_original_mass']=str(1-complement)
    negative=0
    for r in [(Q(1,3),)*4,(Q(3,5),)*3+(Q(1,5),),
              (Q(11,15),Q(7,15),Q(7,15),Q(1,10))]:
        require(not in_core_radius(r) and reference(reference(reference(reference(reference(reference(r))))))==r,
                'Every-state entry is false');negative+=1
    x=inverse_leading(points[-1],(1,-1,1,-1))
    try:leading_target(x,max_steps=2)
    except ValueError as e:require('cap' in str(e),'Cap is a resource failure');negative+=1
    else:raise ArithmeticError('Short evaluation cap was ignored')
    result.update(negative_controls=negative,status='passed',
                  scope='Independent controls for universal leading entry lemmas and actual global joint transfer. Almost-everywhere coverage and weights follow from the analytical exhaustion proof, not these samples.')
    (ROOT/'nine_global_check_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: independent leading and local controls, complete joint-prefix checks, actual captures and negative controls.')


if __name__=='__main__':main()
