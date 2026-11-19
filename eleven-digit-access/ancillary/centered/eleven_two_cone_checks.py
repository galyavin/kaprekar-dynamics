"""Independent primitive, physical, failure-width and actual-path controls."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from random import Random
import json
import eleven_two_cone as cone
from odd_base_covariance import centered_step, centered, in_region
from odd_gap import subtract_pair

ROOT=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def whole_number(B,q):
    digits=(*q,0,0,0,0,0,0);high=low=0
    for a,b in zip(digits,reversed(digits)):
        high=high*B+a;low=low*B+b
    value=high-low;out=[]
    for _ in range(11):
        value,r=divmod(value,B);out.append(r)
    require(value==0,'Whole-number width')
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(5))


def check_step(v,counts):
    x=cone.offsets(v);target,aux=cone.step(v)
    y=centered_step(x);z=centered_step(y)
    require(x[1]+x[3]<=-6 and y[2]>=-3,'Raw deep guards')
    require(z==cone.offsets(target),'Independent raw two-step return')
    require(z[3]-z[4]-x[3]+x[4]==6*aux['current'],'Native fourth-gap current')
    require(z[2]-z[3]-x[2]+x[3]==aux['gap3_change'],'Native third-gap balance')
    require(target[2]<=2,'Primitive preparation')
    if all(Q(q).denominator==1 for q in v):
        require(all(Q(q).denominator==1 for q in target),'Physical integer-coordinate closure')
        require(all((z[i]-(i+1)*z[0])%6==0 for i in range(5)), 'Physical lattice preserved')
    counts['primitive_updates']+=2
    return target,aux


def main():
    require(__debug__,'Optimized Python is not an accepted verification mode')
    counts=Counter();rng=Random(171119)
    for k in range(3,15):
        m=k+1
        for _ in range(20):
            beta=rng.randrange(m)
            x=tuple((i+1)*beta+m*rng.randrange(-30,31) for i in range(k))
            for i in range(1,k-1):
                j=k-1-i
                if i>j:continue
                y=list(x);y[i],y[j]=-x[j]-m,-x[i]-m
                require(centered_step(y)==centered_step(x),'All-dimension interior symmetry')
                require(all((y[q]-(q+1)*beta)%m==0 for q in range(k)),
                        'All-dimension physical base coset')
                counts['all_dimension_reflections']+=1
    for _ in range(16000):
        v=tuple(rng.randrange(-60,61) for _ in range(5));x=cone.offsets(v)
        require(tuple(cone.coordinates(x))==v,'Physical coordinate bijection')
        for indices in ((2,), (1,3)):
            reflected=list(x)
            if len(indices)==1:reflected[2]=-x[2]-6
            else:reflected[1],reflected[3]=-x[3]-6,-x[1]-6
            require(centered_step(reflected)==centered_step(x),'Exact primitive reflection symmetry')
            require(all((reflected[i]-(i+1)*reflected[0])%6==0 for i in range(5)),
                    'Each reflection preserves the physical base coset')
            counts['reflection_controls']+=1
        if cone.guards(v):
            check_step(v,counts);counts['primitive_sources']+=1
    # Exact rational raw controls are distinct from the integer clock claim.
    for _ in range(1000):
        v=tuple(Q(rng.randrange(-60,61),rng.randrange(1,8)) for _ in range(5))
        if cone.guards(v):
            check_step(v,counts);counts['real_primitive_sources']+=1

    for _ in range(20000):
        v=(rng.randrange(-40,81),rng.randrange(-80,21),rng.randrange(-8,3),
           rng.randrange(-20,15),rng.randrange(-8,5))
        work=0;first=None
        for j in range(3):
            if not cone.guards(v):break
            v,aux=check_step(v,counts);work+=aux['current']
            if first is None:first=v
            if work>0:break
        else:
            require(cone.reduced(first),'Three failures force the proved unbounded quadrant')
            counts['failed_prefix_reductions']+=1

    phases=cone.reduced_phases()
    expected={(p,p+b,r) for p in range(-2,3) for b in range(-4,0)
              for r in range(-b-4,1)}
    require(len(phases)==50 and set(phases)==expected,'Complete transverse inventory')
    hist=Counter()
    for p,s,r in phases:
        for a in (-27,0,27,100,10**20):
            for w in (2,0,-27,-100,-10**20):
                v=(a,w,p,s,r);start=cone.fourth_gap(v);work=0
                for h in range(1,19):
                    if not cone.guards(v):
                        counts['quadrant_guard_exits']+=1;break
                    v,aux=check_step(v,counts);work+=aux['current']
                    require(cone.fourth_gap(v)-start==6*work,'Cumulative native-gap current')
                    if work>=2:
                        hist[h]+=1;counts['quadrant_positive_hits']+=1;break
                else:raise ArithmeticError('Accepted quadrant conclusion failed')
                counts['quadrant_sources']+=1
    require(counts['quadrant_sources']==1250,'Complete numerical quadrant control inventory')

    # Actual boundary paths in every base class. The whole-number and column
    # algorithms do not import the new return formula.
    for beta in range(6):
        for label,aa,ww,p,s,r in [('central',0,-1000,0,4,0),
                                  ('inner',6000,2,0,0,-1),
                                  ('interior',6000,-1000,-4,-10,5)]:
            a=aa+(-beta-aa)%6;v=(a,ww,p,s,r);states=[v];raw=[cone.offsets(v)]
            source_gap=cone.fourth_gap(v);hit=None;Cseen=Hseen=False
            for j in range(19):
                require(cone.guards(v),'Entire declared actual boundary path is guarded')
                target,aux=cone.step(v);Cseen|=aux['C']>0;Hseen|=aux['H']>0
                raw.append(centered_step(raw[-1]));raw.append(centered_step(raw[-1]))
                require(raw[-1]==cone.offsets(target),'Prepared whole actual path')
                v=target;states.append(v)
                if hit is None and cone.fourth_gap(v)>source_gap:hit=j+1
            require(hit is not None and hit<=19,'Actual boundary lower hit')
            if label=='central':require(Cseen,'Actual central seam exercised')
            if label=='inner':require(Hseen,'Actual inner seam exercised')
            M=max(abs(q) for x in raw for q in x);B=2*M+6;B+=(beta-B)%6
            q=tuple((B*(5-i)+raw[0][i])//6 for i in range(5))
            for j in range(38):
                require(centered(B,q)==raw[j] and in_region(B,q),'Full actual paired guards')
                q1=subtract_pair(B,q);q2=whole_number(B,q)
                require(q1==q2 and centered(B,q1)==raw[j+1],
                        'Independent whole-number, column and centered agreement')
                q=q1;counts['actual_column_updates']+=1;counts['actual_number_updates']+=1
            counts['actual_paths']+=1

    # False extensions: missing either reflection and dropping physical arithmetic.
    for v,fold in [((0,-1000,0,0,0),'C'),((6000,2,0,0,-1),'H')]:
        target,aux=cone.step(v)
        if fold=='C':
            from eleven_two_cycle_clock import phase_step
            a,w,p,s,r=v;u,(da,dd)=phase_step((p,s,r))
            wrong=(a+6*da,w+dd+4*da,*u)
        else:
            a,w,p,s,r=target;wrong=(a,w,p,s-aux['H'],r+aux['H'])
        raw=centered_step(centered_step(cone.offsets(v)))
        require(raw!=cone.offsets(wrong),'Omitted reflection is a false recurrence')
        counts['negative_controls']+=1
    v=(600,-1000,Q(1,2),Q(-7,6),Q(-5,6));x=cone.offsets(v);start=cone.fourth_gap(v)
    require(all(Q(a).denominator==1 for a in x) and
            any((x[i]-(i+1)*x[0])%6 for i in range(5)),
            'Adverse source is enlarged-integer but not physical')
    for j in range(20):
        require(cone.guards(v),'Adverse finite-window deep guards')
        v,aux=check_step(v,counts)
        require(aux['C']==aux['H']==0 and aux['current']==Q(-2,3),
                'Real arithmetic relaxation has adverse native work')
        require(cone.fourth_gap(v)==start-4*(j+1),'Independent adverse gap balance')
        counts['adverse_returns']+=1
    counts['negative_controls']+=1
    report=dict(status='passed',counts=dict(counts),quadrant_hit_histogram=dict(sorted(hist.items())),
                scope='Independent physical, rational, reflection, reduction, quadrant and actual-boundary controls; no global composition claim.')
    text=json.dumps(report,indent=2)+'\n'
    (ROOT/'eleven_two_cone_results.json').write_text(text)
    (ROOT/'eleven_two_cone_verification.txt').write_text(text)
    print(text)


if __name__=='__main__':main()
