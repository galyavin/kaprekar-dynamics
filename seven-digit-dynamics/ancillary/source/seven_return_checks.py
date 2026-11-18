"""Solver-free guard replay, exact geometry and independent clock checks."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import copy
import json

from seven_return_law import (RETURN_DOMAIN,NEGATIVE_EXIT_DOMAIN,SOURCE_POSITIVE,
    SOURCE_NEGATIVE,POSITIVE_TETRAHEDRA,NEGATIVE_TETRAHEDRA,source_side,
    source_clock,source_arrival,period_region_arrival,leading_bulk,inverse_leading,
    population_masses,bulk_exit,negative_offsets)
from seven_facet import phase_iterate
from seven_boundary import phase_data
from seven_facet_checks import raw_step,raw_arrival
from seven_transport import trace
from odd_affine import word


def require(condition,context):
    if not condition:
        raise AssertionError(context)


def add(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(n,a):return tuple(n*x for x in a)
def sub(a,b):return add(a,mul(-1,b))


def replay(data):
    expected_keys=[f'return:{m}' for m in range(-6,7)]
    expected_keys += [f'exit:{e}:{m}' for e in range(12,16)
                     for m in ((-4,-2) if e%2==0 else (-5,-3,-1))]
    require([r['key'] for r in data['records']]==expected_keys,'case inventory')
    count=0
    for record in data['records']:
        key=record['key'];parts=key.split(':')
        if parts[0]=='return':
            m=int(parts[1]);domain=(*RETURN_DOMAIN,(0,2,0,-1,-m))
            if m<6:domain+=((0,-2,0,1,m),)
            expected=((1,0,0,0,-3),(0,-1,0,1,-3) if m<=-2 else (0,1,0,0,-2),(0,0,1,0,-1))
            length=3
        else:
            e,m=map(int,parts[1:])
            domain=(*NEGATIVE_EXIT_DOMAIN,(2,0,-2,-1,-e),(-2,0,2,1,e),
                    (0,2,0,-1,-m),(0,-2,0,1,m))
            expected=((1,0,0,0,-1),(1,0,-1,0,1),
                      (0,1,-1,0,0) if m>=-1 else (0,-1,-1,1,-1))
            length=1
        require(tuple(map(tuple,record['domain']))==domain,('domain',key))
        B=(0,0,0,1,0);one=(0,0,0,0,1)
        matrix=((1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0));guards=set()
        require(len(record['orders'])==length,('word length',key))
        for order in record['orders']:
            require(sorted(order)==list(range(7)),('sorting permutation',key))
            u,v,w=matrix
            guards.update((sub(w,one),sub(u,v),sub(v,w),sub(sub(B,one),u)))
            digits=(u,v,sub(w,one),sub(B,one),sub(sub(B,w),one),sub(sub(B,v),one),sub(B,u))
            for i,j in zip(order,order[1:]):
                guard=sub(digits[i],digits[j])
                if i>j:guard=sub(guard,one)
                guards.add(guard)
            matrix=tuple(sub(digits[order[i]],digits[order[-1-i]]) for i in range(3))
        require(matrix==expected==tuple(map(tuple,record['matrix'])),('matrix',key))
        if parts[0]=='exit':
            u,v,w=matrix;mm=sub(mul(2,v),B);ss=sub(add(u,w),B)
            guards.update((add(mul(20,one),mm),sub(mul(20,one),mm),
                           add(mul(20,one),ss),sub(mul(20,one),ss),
                           sub(sub(mul(4,u),mul(3,B)),mul(4000,one)),
                           sub(sub(mul(5,B),mul(6,u)),mul(6000,one))))
        saved=[tuple(p['guard']) for p in record['certificates']]
        require(len(saved)==len(set(saved)) and set(saved)==guards,('guard inventory',key))
        for proof in record['certificates']:
            cs=list(map(Q,proof['coefficients']))
            require(len(cs)==len(domain) and all(c>=0 for c in cs),('sign',key))
            require(all(sum(c*r[j] for c,r in zip(cs,domain))==proof['guard'][j]
                        for j in range(5)),('identity',key))
            count+=1
    require(count==data['guard_certificates'],'total count')
    return count


def solve_three(rows):
    a=[[Q(x) for x in (*r[:3],-r[3])] for r in rows]
    for i in range(3):
        pivot=next((j for j in range(i,3) if a[j][i]),None)
        if pivot is None:return None
        a[i],a[pivot]=a[pivot],a[i]
        z=a[i][i];a[i]=[v/z for v in a[i]]
        for j in range(3):
            if j!=i:
                z=a[j][i];a[j]=[x-z*y for x,y in zip(a[j],a[i])]
    return tuple(r[3] for r in a)


def vertices(domain):
    rows=set(tuple(r[i] for i in range(4)) for r in domain)
    rows.update(((1,-1,0,0),(0,1,-1,0),(0,0,1,0),(-1,0,0,1)))
    found=set()
    for triple in combinations(rows,3):
        x=solve_three(triple)
        if x is not None and all(sum(a*b for a,b in zip(r,x))+r[3]>=0 for r in rows):
            found.add(x)
    return found


def center(tetra):
    return tuple(sum(q[i] for q in tetra)/4 for i in range(3))


def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/'seven_return_certificates.json').read_text())
    count=replay(data);negative_controls=0
    for mutation in ('delete','coefficient','domain'):
        broken=copy.deepcopy(data)
        if mutation=='delete':broken['records'][0]['certificates'].pop()
        elif mutation=='domain':broken['records'][0]['domain'][0][4]+=1
        else:
            c=broken['records'][0]['certificates'][0]['coefficients'];c[0]=str(Q(c[0])+1)
        try:replay(broken)
        except AssertionError:negative_controls+=1
        else:raise AssertionError('Corrupt certificate accepted.')
    for domain,tetrahedra in ((SOURCE_POSITIVE,POSITIVE_TETRAHEDRA),(SOURCE_NEGATIVE,NEGATIVE_TETRAHEDRA)):
        require(vertices(domain)==set(q for t in tetrahedra for q in t),'exact vertex inventory')
    masses=population_masses()
    require(masses['period_region']==Q(1339759,43200000),'period region mass')
    require(masses['one_step_predecessor']==Q(937297,5400000),'predecessor mass')
    phases=[]
    for parity in (0,1):
        counts=Counter();table=[]
        for er in range(parity,8,2):
            row=[]
            for mr in range(parity,8,2):
                e0=8000+er;m0=800+mr;N=(e0-12)//4;e=e0-4*N
                m=phase_iterate(m0,N)
                rate=phase_data(parity,negative_offsets(e,m))['time_per_unit_u']
                counts[rate]+=1;row.append(str(rate))
            table.append(row)
        require(counts[Q(3,2)]==(12 if parity==0 else 6),('negative fast weight',parity))
        phases.append({'parity':parity,'rates':table,'fast_weight':str(Q(counts[Q(3,2)],16))})
    bulk_checks=0
    for c,b,a in combinations_with_replacement(range(97),3):
        slacks=(2*a-96,a-b,b-c,c,a+b-96,a-b-c,b+2*c-2*a+96,192-a-2*b-c,96-a)
        if min(slacks)<12:continue
        q=(96+a,a+b,b+c);answer=bulk_exit(192,q);p=q
        for _ in range(answer['elapsed']):p=raw_step(192,p)
        require(p==answer['target'],('bulk endpoint',q))
        try:W=word(192,p,3)
        except ValueError:W=None
        require(W is None or W['orders']!=((3,0,4,1,5,2,6),)*3,('maximality',q))
        bulk_checks+=1
    clocks=[];raw_steps=0;prefix_checks=0
    representatives=[center(t) for t in (*POSITIVE_TETRAHEDRA,*NEGATIVE_TETRAHEDRA)]
    for B in (1000000,1000001,1000006,1000007):
        for x in representatives:
            q=tuple(int(B*v) for v in x)
            require(source_side(B,q) is not None,('source sample',B,q))
            expected=source_arrival(B,q);actual,period=raw_arrival(B,q)
            require((expected['pair_time'],expected['period'])==(actual,period),('raw clock',B,q))
            c=source_clock(B,q)
            require(abs(actual+1-B*c['limiting_time'])<=216,('clock bound',B,q))
            raw_steps+=actual+period
            clocks.append({'base':B,'q':q,'side':c['side'],'pair_time':actual,'period':period,'rate':str(c['rate'])})
            y=x
            for rotation in range(3):
                for signs in product((-1,1),repeat=3):
                    p=tuple(int(B*v) for v in inverse_leading(y,signs))
                    answer=period_region_arrival(B,p,first_image=True)
                    direct=trace(B,p)
                    require((answer['pair_time'],answer['period'])==(direct['arrival'],len(direct['cycle'])),('prefix first capture',B,p))
                    prefix_checks+=1
                y=leading_bulk(y)
    result={'rational_certificates':count,'negative_controls':negative_controls,
            'bulk_first_exit_checks':bulk_checks,'masses':{k:str(v) for k,v in masses.items()},
            'negative_phase_tables':phases,'independent_full_clocks':clocks,
            'independent_digit_steps':raw_steps,'prefix_capture_checks':prefix_checks,
            'scope':'Expanded canonical return and conditional law on the leading period-three region and its one-step predecessors. Global entry uses the separate real-affine lemmas.'}
    (root/'seven_return_results.json').write_text(json.dumps(result,indent=2)+'\n')
    summary=(f'{count} rational guards; {negative_controls} corruption controls; {bulk_checks} maximal bulk exits; '
             f'{len(clocks)} independent full clocks, {raw_steps} digit steps; {prefix_checks} prefix/capture checks.\n')
    (root/'seven_return_verification.txt').write_text(summary)
    print(summary,end='')


if __name__=='__main__':main()
