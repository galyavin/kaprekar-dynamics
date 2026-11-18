"""Rational phase contraction, exact all-base terminal paths and full clocks.

Default acceptance is solver-free. --generate uses Z3 only for rational
coefficient searches and records every sorting word for exact reconstruction.
"""
import argparse,copy,json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

from nine_bbb import (CANONICAL,TERMINALS,CAPTURE_OFFSETS,phase_step,coordinates,inverse_coordinates,
    folded_step,encode,edge_state,phase_entry_bound,exact_arrival)
from nine_entry import ORDER,evaluate,compose
from nine_b_facet import LOW_ORDER
from nine_b_facet_checks import require,from_orders,witness,rational_replay
from nine_phase_checks import step,replay,raw_arrival

ROOT=Path(__file__).resolve().parent
ONE=(0,0,0,0,1);MROW=(0,0,0,1,0)
sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
EDGE=((3,-1,-3),(1,0,0),(-1,1,2),(2,-1,0)) # (V,B,1)
CHART=((3,1,0,0,-1,0),(1,0,0,0,0,0),(-1,0,1,0,1,0),(2,0,0,1,-1,0))


def contraction_domain(small,signs):
    q=((1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0))
    domain=[(1,0,0,0,3),(0,1,0,0,2),(0,0,1,0,1),
            (-1,0,0,1,0),(0,-1,0,1,0),(0,0,-1,1,0),ONE]
    domain+= [(0,0,0,1,-8),(0,0,0,-1,8)] if small else [(0,0,0,1,-9)]
    for high in signs:
        x,y,z=q;domain.append(x if high else neg(x))
        q=(sub(y,ONE) if high else sub(sub(y,x),ONE),sub(z,ONE),
           sub(x,ONE) if high else sub(neg(x),ONE))
    goals=[sub(MROW,r) if small else sub(sub(MROW,ONE),r) for r in q]
    return domain,goals


def terminal_generate(r,e):
    matrix=((4,3*e-r-3),(3,e),(2,r-e+2),(1,2*e-r));seen={};orders=[]
    base=(5,r);one=(0,1)
    while matrix not in seen:
        seen[matrix]=len(orders)
        ds=(*matrix[:-1],sub(matrix[-1],one),sub(base,one),
            *[sub(sub(base,x),one) for x in reversed(matrix[1:])],sub(base,matrix[0]))
        order=tuple(sorted(range(9),key=lambda i:(-evaluate(ds[i],(200000,1)),i)))
        orders.append(order);matrix,_=step(matrix,base,one,order)
        require(len(orders)<10000,'Terminal generation cap')
    return {'base_residue':r,'offset':e,'arrival':seen[matrix],
            'period':len(orders)-seen[matrix],'orders':orders}


def generate():
    from odd_affine import word
    contraction=[]
    for small in (False,True):
        for signs in product((False,True),repeat=6):
            domain,goals=contraction_domain(small,signs)
            contraction.append({'small':small,'signs':signs,
                                'proofs':[witness(g,domain) for g in goals]})
    r=word(1000000,edge_state(1000000,640000),9)
    data={'contraction':contraction,'edge':{'orders':r['orders'],'matrix':r['matrix']},
          'terminals':[terminal_generate(r,e) for r in range(5) for e in (256,257,258)]}
    (ROOT/'nine_bbb_certificates.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    return data


def phase_word_checks():
    """Reconstruct both one-step branches in independent chart coordinates."""
    count=0;extended=(*CHART,(0,0,0,0,1,0),(0,0,0,0,0,1))
    for high,order in ((True,ORDER),(False,LOW_ORDER)):
        M,guards,_=replay(from_orders((order,)))
        target=tuple(compose(row,extended) for row in M)
        # New q from V'=V+a-c+1 and the declared offset recurrence.
        p=((0,-2,0,3,0,-4),
           (0,1,0,-2,0,1) if high else (0,1,-1,-2,0,0),
           (0,-2,-1,2,0,-2))
        V=(1,1,0,-1,0,1);base=(0,0,0,0,1,0)
        expected=(tuple(3*v-b+a for v,b,a in zip(V,base,p[0])),V,
                  tuple(b-v+a for b,v,a in zip(base,V,p[1])),
                  tuple(2*v-b+a for v,b,a in zip(V,base,p[2])))
        require(target==expected,'Actual strip phase formula')
        require(target[0]==tuple(v-(1 if j==5 else 0) for j,v in enumerate(CHART[0])),
                'First coordinate decreases exactly one')
        for g in guards:
            row=compose(g,extended);a,b=row[0],row[4]
            if a or b:
                require(Q(3,5)*a+b>=0 and Q(2,3)*a+b>=0,
                        'Leading guard on the whole open strip')
                require(Q(19,30)*a+b>0,'Strict leading guard in the strip interior')
            else:
                require(row[1]==row[3]==0,'Only the b phase selects a tied order')
                k,c=row[2],row[5]
                require((k>=0 and c-k>=0) if high else (k<=0 and c-2*k>=0),
                        'Phase-only tie guard')
            count+=1
    return count


def conjugacy_checks():
    change=((0,1,0,1),(1,0,-2,3),(2,2,-1,4))
    one=(0,0,0,1)
    inverse=((Q(-4,3),Q(-1,3),Q(2,3),Q(-1,3)),
             (1,0,0,-1),(Q(-2,3),Q(-2,3),Q(1,3),Q(4,3)))
    require(tuple(compose(row,(*change,one)) for row in inverse)
            ==tuple(tuple(int(i==j) for j in range(4)) for i in range(3)),
            'Exact affine inverse of phase coordinates')
    require(tuple(sum(r[j] for r in change) for j in range(4))==(3,3,-3,8),
            'Realizable integer phases have sum two modulo three')
    x,y,z=change
    for high in (False,True):
        p=((-2,0,3,-4),(1,0,-2,1) if high else (1,-1,-2,0),(-2,-1,2,-2))
        target=(sub(y,one) if high else sub(sub(y,x),one),sub(z,one),
                sub(x,one) if high else sub(neg(x),one))
        require(tuple(compose(row,(*p,one)) for row in change)==target,
                'Global piecewise-affine folded conjugacy')


def contraction_checks(records):
    expected=[(small,signs) for small in (False,True) for signs in product((False,True),repeat=6)]
    require([(r['small'],tuple(r['signs'])) for r in records]==expected,'Complete 128 contraction cells')
    count=0
    for record,(small,signs) in zip(records,expected):
        domain,goals=contraction_domain(small,signs)
        count+=rational_replay(record['proofs'],goals,domain)
    return count


def finite_core():
    known={};cycles=[];maximum=0;initials=0;longest=0
    for initial in product(range(-3,9),range(-2,9),range(-1,9)):
        if sum(initial)%3!=2:continue
        initials+=1;q=initial;seen={};path=[]
        while q not in known and q not in seen:
            seen[q]=len(path);path.append(q);q=folded_step(q)
            require(len(path)<1000,'Finite core verification cap')
        if q not in known:
            cycle=path[seen[q]:];cid=len(cycles);cycles.append(cycle)
            for state in cycle:known[state]=(cid,0)
        for state in reversed(path):
            if state not in known:
                cid,n=known[folded_step(state)];known[state]=(cid,n+1)
        longest=max(longest,known[initial][1])
        maximum=max([maximum,*[max(map(abs,s)) for s in path]])
    require((initials,len(known),len(cycles),len(cycles[0]),longest,maximum)==(440,453,1,9,22,10),
            'Finite absorbing-core inventory')
    phases=[inverse_coordinates(q) for q in cycles[0]]
    require(CANONICAL in phases,'Canonical phase lies on the unique cycle')
    require(sum(phase_step(p)[1] for p in phases)==-3,'Nine-step longitudinal decrease')
    # Independent conjugacy and congruence checks on a larger integer cube.
    for p in product(range(-8,9),repeat=3):
        q=coordinates(p);new,_=phase_step(p)
        require(sum(q)%3==2 and inverse_coordinates(q)==p,'Integer chart and congruence')
        require(coordinates(new)==folded_step(q),'Folded phase conjugacy')
        s=q
        for _ in range(3):s=folded_step(s)
        require(s[0]>=-3 and s[1]>=-2 and s[2]>=-1,'Three-step lower preparation')
    return {'initial_states':initials,'visited_states':len(known),'max_transient':longest,
            'max_coordinate':maximum,'cycle':cycles[0],'offset_cycle':phases}


def edge_bound(row):
    a,b,c=row
    slope=b+a*(Q(3,5) if a>=0 else Q(2,3))
    constant=c+256*abs(a)
    require(slope>=0 and 1000000*slope+constant>=0,'Unbounded edge guard')


def verify(data):
    conjugacy_checks()
    result={'rational_contraction_identities':contraction_checks(data['contraction']),
            'strip_guards':phase_word_checks()}
    require(len(data['edge']['orders'])==9,'Nine-step edge word')
    M,guards,offsets=replay(data['edge']);extended=(*EDGE,(0,1,0),(0,0,1))
    require(tuple(compose(row,extended) for row in M)==tuple((a,b,c-3*a) for a,b,c in EDGE),
            'Exact nine-step return')
    for g in guards:edge_bound(compose(g,extended))
    for j,offset in enumerate(offsets):
        first=compose(offset[0],extended)
        require(first==(3,-1,-3-j),'Exact longitudinal clock along the whole edge')
    require([(r['base_residue'],r['offset']) for r in data['terminals']]
            ==[(r,e) for r in range(5) for e in (256,257,258)],'All fifteen terminal families')
    terminal_guards=0;terminal_results=[];errors=[];terminal_error=Q(0)
    for record in data['terminals']:
        r,e=record['base_residue'],record['offset'];base=(5,r);one=(0,1)
        matrix=((4,3*e-r-3),(3,e),(2,r-e+2),(1,2*e-r));seen={};path=[]
        for order in record['orders']:
            require(matrix not in seen,'No premature terminal repeat');seen[matrix]=len(path);path.append(matrix)
            matrix,gs=step(matrix,base,one,order)
            for a,c in gs:
                require(a>=0 and 200000*a+c>=0,'All-base terminal guard');terminal_guards+=1
        require(matrix in seen,'Terminal closure');h=seen[matrix];period=len(path)-h
        require((h,period)==(record['arrival'],record['period'])==TERMINALS[r][e-256],
                'Exact terminal time and period')
        require(all(tuple(row[0] for row in s)==(4,3,2,1) for s in path),
                'Symbolically distinct states are numerically distinct')
        require(all(s[0][1]<=6 for s in path[h:]),'Uniform terminal first-coordinate bound')
        require(tuple(c for k,c in path[h])==CAPTURE_OFFSETS[r][e-256],
                'Exact first cyclic pair in runtime')
        terminal_error=max([terminal_error,*[abs(Q(c)-Q(k*r,5)) for s in path[h:] for k,c in s]])
        errors.append(Q(h-3*e+3)+Q(9*r,5))
        terminal_results.append({'base_residue':r,'offset':e,'pair_time':h,'period':period,
                                 'cycle':path[h:]})
    result.update({'edge_guards':len(guards),'terminal_guards':terminal_guards,
                   'terminals':terminal_results,'pair_error_range':[str(min(errors)),str(max(errors))],
                   'normalized_terminal_error_constant':str(terminal_error)})
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_bbb_certificates.json').read_text())
    result=verify(data);result['finite_phase_core']=finite_core();controls=0
    for kind in ('missing','coefficient','terminal'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['contraction'].pop()
        elif kind=='coefficient':bad['contraction'][0]['proofs'][0]['coefficients'][0]='-1'
        else:bad['terminals'].pop()
        try:verify(bad)
        except (ArithmeticError,AssertionError):controls+=1
        else:raise ArithmeticError('Corrupt triple-B proof accepted')
    paths=[];steps=0
    for B in range(1000000,1000005):
        for p in ((0,0,0),(7,-8,5),(-6,9,-7)):
            q=encode(B,63*B//100,p);answer=exact_arrival(B,q);actual,period=raw_arrival(B,q)
            require((actual,period)==(answer['pair_time'],answer['period']),'Independent complete strip clock')
            steps+=actual+period;paths.append({'base':B,'initial_phase':p,**answer})
    result.update({'status':'passed','corruptions_rejected':controls,'independent_paths':paths,
                   'independent_steps':steps,'scope':'All-base bounded-width triple-B strip; population access remains separate.'})
    (ROOT/'nine_bbb_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {result["rational_contraction_identities"]} exact contraction identities; '
            f'{result["strip_guards"]} one-step strip guards; {result["edge_guards"]} edge guards.\n'
            f'{result["terminal_guards"]} terminal guards cover all fifteen base/entry regimes.\n'
            'Every integer strip phase converges to one nine-cycle; 440 finite-core starts visit 453 states.\n'
            f'15 complete independent clocks, {steps} subtraction steps; {controls} corruptions rejected.\n'
            f'Pair-clock correction range: {result["pair_error_range"]}.\n')
    (ROOT/'nine_bbb_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
