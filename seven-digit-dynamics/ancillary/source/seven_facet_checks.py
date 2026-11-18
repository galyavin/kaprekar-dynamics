"""Solver-free proof replay and independent digit checks for seven transport."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import copy
import json

from odd_gap import ordered_multiplicity
from seven_transport import trace
from seven_facet import (RETURN_DOMAIN,SOURCE_DOMAIN,VERTICES,TETRAHEDRA,
                         source_domain,phase_iterate,edge_entry,exact_arrival,
                         residue_rate,conditional_rate_weights,source_mass,tetra_weight)


def require(condition,context):
    if not condition:
        raise AssertionError(context)


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def verify_certificate(data):
    require([r['case'] for r in data['records']]==list(range(-6,7)),'case inventory')
    count=0
    for record in data['records']:
        case=record['case']
        domain=[*RETURN_DOMAIN,(0,2,0,-1,-case)]
        if case<6:
            domain.append((0,-2,0,1,case))
        require([tuple(r) for r in record['domain']]==domain,('domain',case))
        one,base=(0,0,0,0,1),(0,0,0,1,0)
        matrix=((1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0))
        guards=set()
        require(len(record['orders'])==3,('word length',case))
        for order in record['orders']:
            require(sorted(order)==list(range(7)),('order inventory',case))
            u,v,w=matrix
            guards.update((sub(w,one),sub(u,v),sub(v,w),sub(sub(base,one),u)))
            ds=(u,v,sub(w,one),sub(base,one),sub(sub(base,w),one),sub(sub(base,v),one),sub(base,u))
            for i,j in zip(order,order[1:]):
                guard=sub(ds[i],ds[j])
                if i>j:
                    guard=sub(guard,one)
                guards.add(guard)
            matrix=tuple(sub(ds[order[i]],ds[order[-1-i]]) for i in range(3))
        expected=((1,0,0,0,-3),(0,-1,0,1,-3) if case<=-2 else (0,1,0,0,-2),(0,0,1,0,-1))
        require(matrix==expected==tuple(map(tuple,record['matrix'])),('return matrix',case))
        saved=[tuple(p['guard']) for p in record['certificates']]
        require(len(saved)==len(set(saved)) and set(saved)==guards,('guard inventory',case))
        for proof in record['certificates']:
            coefficients=list(map(Q,proof['coefficients']))
            require(len(coefficients)==len(domain) and all(c>=0 for c in coefficients),('coefficient sign',case))
            require(all(sum(c*r[j] for c,r in zip(coefficients,domain))==proof['guard'][j] for j in range(5)),('coefficient identity',case))
            count+=1
    require(count==data['guard_certificates'], 'total obligation count')
    return count


def raw_step(base,q):
    # Numeric subtraction of a seven-slot representative, not borrowed forms.
    ds=sorted((*q,0,0,0,0))
    high=low=0
    for d in ds:
        low=base*low+d
    for d in reversed(ds):
        high=base*high+d
    value,out=high-low,[]
    for _ in range(7):
        value,d=divmod(value,base);out.append(d)
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(3))


def raw_arrival(base,q):
    seen={}
    while q not in seen:
        seen[q]=len(seen);q=raw_step(base,q)
    return seen[q],len(seen)-seen[q]


def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/'seven_facet_certificates.json').read_text())
    certificates=verify_certificate(data)
    negative_controls=0
    for mutation in ('delete','coefficient'):
        broken=copy.deepcopy(data)
        if mutation=='delete':
            broken['records'][0]['certificates'].pop()
        else:
            c=broken['records'][0]['certificates'][0]['coefficients']
            c[0]=str(Q(c[0])+1)
        try:
            verify_certificate(broken)
        except AssertionError:
            negative_controls+=1
        else:
            raise AssertionError('Broken certificate accepted.')
    phase_checks=0
    for m in range(-6,401):
        p=m
        for n in range(101):
            require(phase_iterate(m,n)==p,('phase iterate',m,n))
            p=abs(p+1)-5;phase_checks+=1
    rate_tables=[]
    for parity in (0,1):
        counts=Counter();table=[]
        for sm in range(8):
            row=[]
            for mm in range(parity,8,2):
                rate=residue_rate(parity,sm,mm)
                counts[rate]+=1;row.append(str(rate))
            table.append(row)
        require({r:Q(c,32) for r,c in counts.items()}==conditional_rate_weights(parity),('phase weights',parity))
        rate_tables.append({'parity':parity,'rates_by_S_mod8_and_m_mod8':table,
                            'weights':{str(r):str(Q(c,32)) for r,c in counts.items()}})
    # Independently check the exact polynomial weight and simplex normalization.
    for g1,g2,g3 in product(range(1,8),repeat=3):
        B=20+g1+g2+g3;q=(g1+g2+g3,g2+g3,g3)
        expected=20*(5040*g1*g2*g3-840*g1-840*g3)
        require(ordered_multiplicity(B,q)==expected,('weight',q))
    unit=((Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(1),Q(1),Q(0)),(Q(1),Q(1),Q(1)))
    require(tetra_weight(unit)==1,'density normalization')
    mass=source_mass()
    require(mass==Q(109859671,18662400000),('source mass',mass))
    # Independently check first capture even inside compressed exploratory words.
    exploration_rng=Random(729311)
    compressed_checks=0
    for B in (17,18,31,32,101,102):
        for _ in range(30):
            q=tuple(sorted(exploration_rng.sample(range(1,B),3),reverse=True))
            compressed=trace(B,q)
            require((compressed['arrival'],len(compressed['cycle']))==raw_arrival(B,q),
                    ('compressed first capture',B,q))
            compressed_checks+=1
    rng=Random(61742026)
    checks=[];direct_steps=0
    for B in (100000,100001,100006,100007,1000000,1000001):
        for index in range(6):
            q=(rng.randrange(83*B//100,84*B//100),
               rng.randrange(51*B//100,103*B//200),
               rng.randrange(23*B//100,24*B//100))
            require(source_domain(B,q),('source membership',B,q))
            result=exact_arrival(B,q)
            actual,period=raw_arrival(B,q)
            require((actual,period)==(result['pair_time'],result['period']),('direct arrival',B,q,result,actual,period))
            direct_steps+=actual+period
            entry=edge_entry(B,q)
            S,m=q[0]+q[2]-B,2*q[1]-B
            require(residue_rate(B%2,S%8,m%8)==entry['rate'],('residue classifier',B,q))
            tau=Q(3,4)*(Q(q[0]+q[2],B)-1)+entry['rate']*Q(q[0]-3*q[2],4*B)
            require(abs(actual+1-B*tau)<=216,('clock bound',B,q,actual,tau))
            shifted=(q[0],q[1]+2,q[2]+2)
            shifted_result=exact_arrival(B,shifted)
            if B%2==0 and S%8 in (1,2,5,6) and B>=200000:
                require(abs(shifted_result['actual_time']-result['actual_time'])>=Q(B,100),('sensitivity',B,q))
            checks.append({'base':B,'q':q,'actual_time':actual+1,'period':period,
                           'edge_rate':str(entry['rate']),'shifted_actual_time':shifted_result['actual_time']})
    result={'rational_guard_certificates':certificates,'negative_controls':negative_controls,
            'phase_iterate_checks':phase_checks,'source_mass':str(mass),'source_mass_float':float(mass),
            'rate_tables':rate_tables,'direct_arrival_checks':checks,'raw_digit_steps':direct_steps,
            'compressed_exploration_checks':compressed_checks,
            'scope':'All-base return guards replayed by exact rational arithmetic; finite traces independently corroborate the complete clocks. The global law also requires leading entry.'}
    (root/'seven_facet_results.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=[f'{certificates} rational guard certificates replayed without a solver; {negative_controls} negative controls rejected.',
           f'{phase_checks} phase iterates; exact source mass {mass}.',
           f'{len(checks)} independent full digit trajectories, {direct_steps} raw steps.',
           f'{compressed_checks} independent compressed exploratory first-capture checks.',
           'Conditional fast/slow weights: even bases 1/2,1/2; odd bases 3/8,5/8.',
           'Facet replay complete.']
    (root/'seven_facet_verification.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__=='__main__':
    main()
