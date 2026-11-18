"""Independent finite populations and exact-prefix checks for the global law.

Finite distributions and quadrature are diagnostics, not the entry proof.
The unbounded affine lemmas are checked separately by seven_entry_audit.py.
"""
from array import array
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from random import Random
import json
import math

from odd_gap import ordered_multiplicity
from seven_facet_checks import raw_step,raw_arrival
from seven_global import limit_kernel,capture_arrival


def require(ok,context):
    if not ok:raise ArithmeticError(context)


def index(q):
    u,v,w=q
    return u*(u+1)*(u+2)//6+v*(v+1)//2+w-1


def finite_population(B):
    size=B*(B+1)*(B+2)//6-1
    times=array('i',[-1])*size
    for u in range(1,B):
        for v in range(u+1):
            for w in range(v+1):
                q=(u,v,w)
                if times[index(q)]>=0:continue
                path=[];seen={};p=q
                while times[index(p)]<0 and p not in seen:
                    seen[p]=len(path);path.append(p);p=raw_step(B,p)
                    require(p[0]>0,('nonconstant image',B,q))
                if p in seen:
                    cut=seen[p]
                    for state in path[cut:]:times[index(state)]=0
                    path=path[:cut]
                h=times[index(p)]
                for state in reversed(path):h+=1;times[index(state)]=h
    histogram=Counter();mass=0;cycles=0
    for u in range(1,B):
        for v in range(u+1):
            for w in range(v+1):
                if u>v>w>0:
                    weight=(B-u)*(5040*(u-v)*(v-w)*w-840*(u-v)-840*w)
                else:weight=ordered_multiplicity(B,(u,v,w))
                h=times[index((u,v,w))]
                require(h>=0,('unresolved state',B,u,v,w))
                histogram[h+1]+=weight;mass+=weight
                if h==0:
                    # Exactly one numerical cycle input has this pair state.
                    histogram[1]-=1;histogram[0]+=1;cycles+=1
    require(mass==B**7-B==sum(histogram.values()),('population normalization',B,mass))
    mean=sum(t*w for t,w in histogram.items())/(B*mass)
    transforms={str(s):sum(w*math.exp(-s*t/B) for t,w in histogram.items())/mass
                for s in (1,4,10)}
    return {'base':B,'pair_states':size,'cycle_states':cycles,'input_mass':mass,
            'mean_T_over_B':mean,'laplace':transforms,
            'maximum_T_over_B':max(histogram)/B}


def main():
    root=Path(__file__).resolve().parent
    rng=Random(7162026);checks=[];rejections=0;steps=0
    for B in (1000000,1000001,1000006,1000007):
        accepted=0
        for _ in range(100):
            ds=sorted(rng.randrange(B) for _ in range(7))
            q=(ds[6]-ds[0],ds[5]-ds[1],ds[4]-ds[2])
            try:answer=capture_arrival(B,q)
            except ValueError:rejections+=1;continue
            h,period=raw_arrival(B,q)
            require((answer['pair_time'],answer['period'])==(h,period),('global prefix',B,q))
            steps+=h+period
            checks.append({'base':B,'q':q,'leading_prefix':answer['leading_prefix'],
                           'pair_time':h,'period':period,'side':answer['side']})
            accepted+=1
            if accepted==3:break
        require(accepted==3,('insufficient buffered samples',B))
    # Fixed rational Monte Carlo approximation of the continuous input model.
    # Its numerical agreement is not a proof or a certified integration error.
    sums={p:{'mean':0.0,**{str(s):0.0 for s in (1,4,10)}} for p in (0,1)}
    sample_count=4000;denominator=1000000007;seams=0;prefixes=Counter()
    for _ in range(sample_count):
        ds=sorted(rng.randrange(denominator) for _ in range(7))
        x=tuple(Q(ds[-1-i]-ds[i],denominator) for i in range(3))
        try:k=limit_kernel(x,0)
        except ValueError:seams+=1;continue
        f,t=float(k['fast']),float(k['slow']);prefixes[k['leading_prefix']]+=1
        require(0<f<t<0.25,('limiting support',x,k))
        for parity in (0,1):
            p=float(k['fast_weight']) if parity==0 else 3/8
            sums[parity]['mean']+=p*f+(1-p)*t
            for s in (1,4,10):sums[parity][str(s)]+=p*math.exp(-s*f)+(1-p)*math.exp(-s*t)
    quadrature={str(p):{k:v/(sample_count-seams) for k,v in data.items()} for p,data in sums.items()}
    finite=[]
    for B in (96,97,192,193):
        record=finite_population(B);finite.append(record)
        print('finite population',B,record['pair_states'],record['mean_T_over_B'],flush=True)
    result={'independent_prefix_clocks':checks,'buffered_sample_rejections':rejections,
            'independent_digit_steps':steps,'quadrature_samples':sample_count,
            'quadrature_seams':seams,'quadrature':quadrature,'leading_prefix_counts':dict(prefixes),
            'finite_populations':finite,
            'scope':'Exact finite input distributions and independent first-capture checks. Monte Carlo values have no certified integration error; global convergence rests on the separate entry and weighted-exhaustion proof.'}
    (root/'seven_global_results.json').write_text(json.dumps(result,indent=2)+'\n')
    summary=(f'{len(checks)} independent leading-prefix clocks, {steps} digit steps; '
             f'{sum(r["pair_states"] for r in finite)} complete pair states in four bases; '
             f'{sample_count-seams} rational quadrature samples.\n')
    (root/'seven_global_verification.txt').write_text(summary)
    print(summary,end='')


if __name__=='__main__':main()
