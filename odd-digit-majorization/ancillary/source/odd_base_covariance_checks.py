"""Independent checks for base covariance; universal proof is in the note."""
import json
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
from random import Random

from odd_base_covariance import (anchor,in_region,paired_step,centered,centered_step,
                                 paired_orders,order_is_paired)

ROOT=Path(__file__).resolve().parent


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def raw(B,q):
    # Full integer subtraction, including all leading-zero slots.
    k=len(q);digits=sorted((*q,*(0 for _ in range(k+1))))
    small=large=0
    for d in digits:small=B*small+d
    for d in reversed(digits):large=B*large+d
    v=large-small;out=[]
    for _ in range(2*k+1):
        v,d=divmod(v,B);out.append(d)
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(k))


def symbolic_digits(k):
    # Rows in (q_1,...,q_k,B,1), reconstructed independently.
    unit=lambda j:tuple(int(i==j) for i in range(k+2))
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    q=[unit(i) for i in range(k)];B=unit(k);one=unit(k+1)
    return (*q[:-1],sub(q[-1],one),sub(B,one),
            *[sub(sub(B,x),one) for x in reversed(q[1:])],sub(B,q[0]))


def verify_orders(k,orders,v):
    require(set(map(tuple,orders))==set(paired_orders(k)),'Complete paired order inventory')
    ds=symbolic_digits(k);direction=(*v,k+1,0);guards=0
    ev=lambda a:sum(x*y for x,y in zip(a,direction))
    for order in orders:
        for i,j in zip(order,order[1:]):
            require(ev(ds[i])>=ev(ds[j]),'Nonnegative sorting-guard base slope')
            guards+=1
        require(tuple(ev(ds[order[i]])-ev(ds[order[-1-i]]) for i in range(k))==tuple(v),
                'Common fixed ray in every affine branch')
    return guards


def formula_checks():
    counts={};noninvariance=[]
    for k in range(2,7):
        B=14 if k<=4 else 10;n=0;witness=None;v=anchor(k)
        for asc in combinations_with_replacement(range(1,B),k):
            q=tuple(reversed(asc))
            if not in_region(B,q):continue
            r=raw(B,q);translated=tuple(x+a for x,a in zip(q,v))
            require(paired_step(B,q)==r,'Paired formula vs full integer subtraction')
            require(centered_step(centered(B,q))==centered(B,r),'Autonomous centered formula')
            require(in_region(B+k+1,translated),'Translation preserves paired guards')
            require(raw(B+k+1,translated)==tuple(x+a for x,a in zip(r,v)),
                    'Exact one-step base translation')
            n+=1
            if not in_region(B,r) and witness is None:witness={'k':k,'base':B,'q':q,'image':r}
        require(witness is not None,'Do not mistake paired region for an invariant region')
        counts[k]=n;noninvariance.append(witness)
    rng=Random(176)
    cosets=0
    for k in range(1,13):
        m=k+1;v=anchor(k)
        for beta in range(m):
            for _ in range(30):
                x=tuple(m*rng.randrange(-1000,1001)-a*beta for a in v)
                y=centered_step(x)
                require(all((z+a*beta)%m==0 for z,a in zip(y,v)),'Each base coset is preserved')
                cosets+=1
    return {'exhaustive_paired_inputs':counts,'centered_coset_checks':cosets,
            'noninvariance_witnesses':noninvariance}


def extract_orders(o):
    if isinstance(o,list):
        if len(o)==9 and all(isinstance(x,int) for x in o) and sorted(o)==list(range(9)):
            yield tuple(o)
        else:
            for v in o:yield from extract_orders(v)
    elif isinstance(o,dict):
        for v in o.values():yield from extract_orders(v)


def nine_inventory():
    import nine_entry as bulk
    import nine_a_facet as af
    import nine_b_facet as bf
    import nine_ab as ab
    import nine_aa as aa
    import nine_bb as bb
    import nine_aaa as aaa
    orders={bulk.ORDER};sources={}
    for words in af.THREE_ORDERS.values():orders.update(words)
    for words in bf.THREE_ORDERS.values():orders.update(words)
    for words in bb.THREE_ORDERS:orders.update(words)
    for key in ab.KEYS:orders.update(ab.three_orders(key))
    for key in aa.KEYS:orders.update(aa.three_orders(key))
    for phase in ab.PHASES:orders.update(ab.return_orders(phase))
    for phase in aa.PHASES:orders.update(aa.return_spec(phase)[0])
    orders.update(aaa.ORDERS.values())
    for name in ('nine_ab','nine_aa','nine_aab','nine_aaa','nine_a_facet',
                 'nine_b_facet','nine_bb','nine_bbb','nine_releases'):
        found=set(extract_orders(json.loads((ROOT/(name+'_certificates.json')).read_text())))
        sources[name]=len(found);orders.update(found)
    require(all(order_is_paired(o) for o in orders),'Every inherited core sorting word is paired')
    # Independent nine-digit symbolic guard constructor.
    from nine_phase_checks import replay
    from nine_b_facet_checks import from_orders
    v=(*anchor(4),5)
    checks=0
    for order in orders:
        M,guards,_=replay(from_orders((order,)))
        require(tuple(sum(a*b for a,b in zip(row[:5],v)) for row in M)==anchor(4),
                'Inherited nine-slot affine base covariance')
        for g in guards:
            require(sum(a*b for a,b in zip(g[:5],v))>=0,'Inherited guard monotonicity')
            checks+=1
    return {'certificates':sources,'distinct_orders':len(orders),'orders':sorted(orders),
            'independent_guard_checks':checks}


def actual_orbit_checks():
    from nine_phase_checks import raw_step
    captures=json.loads((ROOT/'nine_selection_results.json').read_text())['actual']['captures']
    total=0;out=[];v=anchor(4)
    for rec in captures:
        B=rec['base'];q=tuple(rec['source']);shifted=tuple(a+b for a,b in zip(q,v));seen={}
        while q not in seen:
            require(in_region(B,q),'Whole actual source orbit remains paired')
            require(in_region(B+5,shifted),'Whole translated orbit remains paired')
            require(shifted==tuple(a+b for a,b in zip(q,v)),'Pointwise orbit translation')
            seen[q]=len(seen)
            q=raw_step(B,q);shifted=raw_step(B+5,shifted)
        require(shifted==tuple(a+b for a,b in zip(q,v)),'Translated cycle closes')
        h=seen[q];period=len(seen)-h
        require((h,period,q)==(rec['pair_time'],rec['period'],tuple(rec['first_cycle_state'])),
                'First capture and period preserved')
        total+=2*len(seen)
        out.append({'kind':rec['kind'],'base':B,'translated_base':B+5,
                    'pair_time':h,'period':period,'translated_first_cycle_state':shifted})
    return {'raw_pair_steps':total,'translated_captures':out}


def negative_controls():
    bads=[('wrong anchor',list(paired_orders(4)),(5,3,2,1)),
          ('missing sorting branch',list(paired_orders(4))[1:],anchor(4))]
    for label,orders,v in bads:
        try:verify_orders(4,orders,v)
        except ArithmeticError:continue
        raise ArithmeticError('Accepted corrupted '+label)
    return [b[0] for b in bads]


def main():
    result={'symbolic_guard_checks':{k:verify_orders(k,list(paired_orders(k)),anchor(k))
                                     for k in range(1,9)},
            'formulas':formula_checks(),'nine_digit_inventory':nine_inventory(),
            'actual':actual_orbit_checks(),'rejected_controls':negative_controls()}
    (ROOT/'odd_base_covariance_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
