"""Complete AA exit inventory, strip phases and all-base first pair capture."""
import argparse,copy,json
from fractions import Fraction as Q
from collections import Counter
from math import ceil
from pathlib import Path
import nine_aa as aa
from nine_aaa import (THRESHOLDS,INCOMING,CANONICAL,DECREMENTS,ORDERS,EARLY_CAPTURE,encode,
    phase_step,phase_path,branch,edge_orders,exact_pair_arrival,terminal_table,from_recurrent_aa)
from nine_entry import compose,evaluate
from nine_b_facet_checks import require,from_orders
from nine_phase_checks import replay,step,raw_step

ROOT=Path(__file__).resolve().parent
CHART=((1,0,0,0,0,0),(2,1,0,0,-1,0),(3,0,1,0,-2,0),(-1,0,0,1,1,0))
EXTENDED=(*CHART,(0,0,0,0,1,0),(0,0,0,0,0,1))
INCOMING_CHART=((1,0,0,0,0,0),(2,0,-1,0,-1,0),(3,0,-1,-1,-2,0),
                (-1,1,0,0,1,0),(0,0,0,0,1,0),(0,0,0,0,0,1))


def seed(period,r,e):
    a,b,c=CANONICAL[period]
    return (4,e),(3,2*e-r+a),(2,3*e-2*r+b),(1,r-e+c)


def terminal_generate(period,r,e,initial=None):
    matrix=seed(period,r,e) if initial is None else initial
    seen={};path=[];orders=[];base=(5,r);one=(0,1)
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    while matrix not in seen:
        seen[matrix]=len(path);path.append(matrix)
        ds=(*matrix[:-1],sub(matrix[-1],one),sub(base,one),
            *[sub(sub(base,x),one) for x in reversed(matrix[1:])],sub(base,matrix[0]))
        order=tuple(sorted(range(9),key=lambda i:(-evaluate(ds[i],(200000,1)),i)))
        orders.append(order);matrix,_=step(matrix,base,one,order)
        require(len(orders)<10000,'Terminal generation cap')
    h=seen[matrix]
    return {'edge_period':period,'base_residue':r,'offset':e,'arrival':h,'period':len(path)-h,
            'capture_offsets':[c for k,c in path[h]],'cycle_u_max':max(s[0][1] for s in path[h:]),
            'orders':orders}


def generate():
    records=[terminal_generate(period,r,e) for period,d in DECREMENTS.items()
             for r in range(5) for e in range(256,256+d)]
    data={'terminals':records}
    (ROOT/'nine_aaa_certificates.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    table=[{k:v for k,v in r.items() if k!='orders'} for r in records]
    (ROOT/'nine_aaa_terminals.json').write_text(json.dumps(table,indent=2)+'\n')
    terminal_table.cache_clear();return data


def strip_checks():
    count=0
    for key,order in ORDERS.items():
        M,guards,_=replay(from_orders((order,)));actual=tuple(compose(row,EXTENDED) for row in M)
        if key=='low':phase=((0,0,0,1,0,-1),(0,1,0,3,0,1),(0,1,-1,-1,0,0));du=(0,0,0,-1,0,0)
        elif key=='middle':phase=((0,0,0,2,0,0),(0,1,0,2,0,0),(0,1,-1,-1,0,0));du=(0,0,0,-1,0,0)
        else:phase=((0,0,0,-1,0,3),(0,1,0,-1,0,3),(0,1,-1,0,0,-1));du=(0,0,0,0,0,-1)
        V=tuple((j==0)+x for j,x in enumerate(du));B=(0,0,0,0,1,0)
        expected=(V,tuple(2*u-b+a for u,b,a in zip(V,B,phase[0])),
                  tuple(3*u-2*b+a for u,b,a in zip(V,B,phase[1])),
                  tuple(b-u+a for b,u,a in zip(B,V,phase[2])))
        require(actual==expected,'Exact AAA strip phase formula')
        for g in guards:
            row=compose(g,EXTENDED);a,b=row[0],row[4]
            if a or b:
                require(Q(4,5)*a+b>=0 and Q(5,6)*a+b>=0 and Q(49,60)*a+b>0,
                        'Whole open AAA strip leading guard')
            else:
                require(row[1]==row[2]==0,'Only the c phase selects tied sorting orders')
                k,c=row[3],row[5]
                ok=(k<=0 and c-2*k>=0) if key=='low' else (
                    min(c-k,c)>=0 if key=='middle' else k>=0 and k+c>=0)
                require(ok,'Exact phase-only tie guard')
            count+=1
    return count


def incoming_checks():
    records=[];zero_guards=0
    for p,r in aa.PHASES:
        orders,delta=aa.return_spec((p,r));_,guards,_=replay(from_orders(orders));ts=[]
        for g in guards:
            row=compose(g,INCOMING_CHART)
            if row[0] or row[4]:
                a,b=row[0],row[4]
                require(Q(4,5)*a+b>=0 and Q(5,6)*a+b>=0 and Q(49,60)*a+b>0,
                        'Other AA exit guards have strict leading margin')
                continue
            c=row[2]*p+row[3]*r+row[5];e=row[1];zero_guards+=1
            if e>0:ts.append(ceil(-c/e))
            else:require(e==0 and c>=0,'No omitted competing transverse exit guard')
        records.append((p,r,max(ts),-delta[0]-delta[3]))
    require(tuple(records)==THRESHOLDS and len(INCOMING)==len(set(INCOMING))==435,
            'Exact complete AA exit inventory')
    return zero_guards


def phase_checks():
    counts=Counter();maximum=shift_bound=0;cycles=set();visited=set();largest=0
    for initial in INCOMING:
        path=phase_path(initial);states=path['states'];shifts=path['shifts'];n=len(states)
        require(len(set(states))==n,'No premature phase repeat')
        for i,p in enumerate(states):
            q,du=phase_step(p)
            require(shifts[i+1]==shifts[i]+du,'Every phase longitudinal displacement')
            require(q==(states[i+1] if i+1<n else states[n-path['period']]),'Complete finite phase path')
        h=n-path['period'];cycle=states[h:];k=min(range(len(cycle)),key=lambda i:cycle[i]);cycle=cycle[k:]+cycle[:k]
        cycles.add((cycle,shifts[-1]-shifts[h]));counts[path['period']]+=1
        maximum=max(maximum,h);shift_bound=max(shift_bound,*map(abs,shifts));visited.update(states)
        largest=max(largest,*[max(map(abs,p)) for p in states])
    require(sorted((len(c),d) for c,d in cycles)==[(1,0),(2,0),(6,0),(32,-2),(37,-3)],'Five complete reached cycles')
    require(dict(counts)=={37:238,6:42,2:45,32:96,1:14} and maximum==62 and largest==33,
            'Finite source classification inventory')
    return {'initials':435,'visited_states':len(visited),'maximum_transient':maximum,
            'maximum_offset':largest,'maximum_longitudinal_shift':shift_bound,'source_counts_by_period':dict(counts),
            'cycles':[{'phases':c,'delta_U':d} for c,d in sorted(cycles,key=lambda t:len(t[0]))]}


def edge_checks():
    counts=0;bounds={}
    for period,p in CANONICAL.items():
        d=DECREMENTS[period];M,guards,offsets=replay(from_orders(edge_orders(period)))
        chart=((1,0,0),(2,-1,p[0]),(3,-2,p[1]),(-1,1,p[2]),(0,1,0),(0,0,1))
        target=tuple(compose(row,chart) for row in M)
        require(target==tuple((a,b,c-d*a) for a,b,c in chart[:4]),'Exact drifting AAA return')
        for g in guards:
            a,b,c=compose(g,chart);slope=b+a*(Q(4,5) if a>=0 else Q(5,6));constant=c+(252*a if a>=0 else -256*a)
            require(slope>=0 and 1000000*slope+constant>=0,'Unbounded actual AAA edge guard')
            counts+=1
        us=[compose(m[0],chart) for m in offsets]
        require(all(a==1 and b==0 for a,b,c in us),'Bounded longitudinal edge offsets')
        bounds[period]=min(c for a,b,c in us)
    return counts,bounds


def seam_checks():
    """Exclude all earlier AA repeats, not just the last sampled repeat."""
    found={};comparisons=0
    for initial in INCOMING:
        path=phase_path(initial);a,b,c=initial;orders,delta=aa.return_spec((-a,a-b))
        _,_,offsets=replay(from_orders(orders))
        base=((1,0,0),(2,-1,a),(3,-2,b),(-1,1,c))
        start=tuple((r[0],r[1],r[2]-d) for r,d in zip(base,delta));ext=(*start,(0,1,0),(0,0,1))
        if path['stationary']:
            h=path['prefix'];cycle=[encode(0,path['shifts'][i],path['states'][i]) for i in range(h,h+path['period'])]
        for j,off in enumerate(offsets):
            state=[compose(row,ext) for row in off]
            direction=[-evaluate(row[:4],delta) for row in off]
            require([r[:2] for r in state]==[r[:2] for r in base], 'Same leading terminal-line coefficients')
            require(direction[0]>0,'Earlier AA repeats have strictly larger first coordinate')
            if not path['stationary']:continue
            for end in cycle:
                comparisons+=1;n=Q(end[0]-state[0][2],direction[0])
                if n>=0 and n.denominator==1 and all(r[2]+n*d==v for r,d,v in zip(state,direction,end)):
                    require(n==0 and path['prefix']==0,'Only the final AA word can capture early')
                    found.setdefault(initial,(len(orders)-j,tuple(r[2] for r in state)))
    require(found==EARLY_CAPTURE,'Complete exact seam-correction inventory')
    return comparisons


def verify(data):
    expected=[(period,r,e) for period,d in DECREMENTS.items() for r in range(5) for e in range(256,256+d)]
    require([(r['edge_period'],r['base_residue'],r['offset']) for r in data['terminals']]==expected,
            'Complete twenty-five all-base terminal families')
    result={'strip_guards':strip_checks(),'incoming_transverse_guards':incoming_checks(),
            'finite_phase_graph':phase_checks(),'seam_comparisons':seam_checks()}
    edge_count,bounds=edge_checks();count=0;tables=[];time_errors=[];terminal_error=Q(0)
    for record,(P,r,e) in zip(data['terminals'],expected):
        matrix=seed(P,r,e);seen={};path=[];base=(5,r);one=(0,1)
        for order in record['orders']:
            require(matrix not in seen,'First terminal repeat inventory');seen[matrix]=len(path);path.append(matrix)
            matrix,guards=step(matrix,base,one,order)
            for a,c in guards:
                require(a>=0 and 200000*a+c>=0,'All-base affine terminal guard');count+=1
        require(matrix in seen,'Terminal closure');h=seen[matrix];period=len(path)-h
        require(all(tuple(a for a,c in m)==(4,3,2,1) for m in path),'Distinct affine states remain distinct for all bases')
        capture=[c for a,c in path[h]];umax=max(m[0][1] for m in path[h:])
        require((h,period,capture,umax)==(record['arrival'],record['period'],record['capture_offsets'],record['cycle_u_max']),
                'Exact first capture and cycle range')
        require(e+bounds[P]>umax,'Every earlier repeated-edge state is outside the eventual cycle')
        table={k:v for k,v in record.items() if k!='orders'};tables.append(table)
        time_errors.append(Q(h)-Q(P,DECREMENTS[P])*(Q(e)-Q(4*r,5)))
        terminal_error=max(terminal_error,*[abs(Q(c)-Q(a*r,5)) for a,c in path[h]])
    require(tables==json.loads((ROOT/'nine_aaa_terminals.json').read_text()),'Runtime terminal table matches proved paths')
    result.update({'edge_guards':edge_count,'terminal_guards':count,'terminals':tables,
                   'edge_clock_correction_range':list(map(str,(min(time_errors),max(time_errors)))),
                   'first_terminal_error_constant':str(terminal_error)})
    return result


def raw_capture(B,q):
    seen={}
    while q not in seen:seen[q]=len(seen);q=raw_step(B,q)
    return seen[q],len(seen)-seen[q],q,len(seen)


def actual_checks():
    compressed=0;representatives={};paths=[];steps=0
    for phase in INCOMING:
        path=phase_path(phase);representatives.setdefault(path['period'],phase)
        for B in range(1000000,1000005):
            q=encode(B,41*B//50,phase);answer=exact_pair_arrival(B,q)
            require(answer['phase_period']==path['period'],'Preserved incoming phase class');compressed+=1
    for P,p in representatives.items():
        bases=range(1000000,1000005) if P in (32,37) else (1000000,1000001)
        for B in bases:
            q=encode(B,41*B//50,p);answer=exact_pair_arrival(B,q);h,period,captured,n=raw_capture(B,q)
            require((h,period,captured)==(answer['pair_time'],answer['period'],answer['first_cycle_state']),
                    'Independent exact first pair capture')
            steps+=n;paths.append({'base':B,'source':q,**answer})
    return {'compressed_captures':compressed,'independent_paths':paths,'independent_steps':steps}


def actual_seam_checks():
    checks=0;steps=0;compressed=0;drifting={};full=[]
    # Every stationary exit, including all eight early-capture cases, is
    # joined to three complete preceding words and independently traced.
    for phase in INCOMING:
        if not phase_path(phase)['stationary']:continue
        a,b,c=phase;orders,delta=aa.return_spec((-a,a-b));B=1000000
        target=encode(B,41*B//50,phase);source=tuple(x-3*d for x,d in zip(target,delta))
        result=from_recurrent_aa(B,source);h,period,captured,n=raw_capture(B,source)
        require((h,period,captured)==(result['pair_time'],result['period'],result['first_cycle_state']),
                'Exact first capture across the preceding AA seam')
        steps+=n;checks+=1
    for phase in aa.PHASES:
        for B in range(1000000,1000005):
            source=aa.encode(B,41*B//50,21*B//100,phase);answer=from_recurrent_aa(B,source);compressed+=1
            P=answer['tail']['phase_period']
            if P in (32,37):drifting.setdefault((P,B%2),(B,source,answer))
    require(len(drifting)==4,'Both drifting classes in both base parities')
    for B,source,answer in drifting.values():
        h,period,captured,n=raw_capture(B,source)
        require((h,period,captured)==(answer['pair_time'],answer['period'],answer['first_cycle_state']),
                'Independent full AA-to-terminal capture')
        full.append({'base':B,'source':source,'pair_time':h,'period':period,
                     'first_cycle_state':captured,'phase_period':answer['tail']['phase_period']});steps+=n
    return {'joined_stationary_captures':checks,'joined_compressed_captures':compressed,
            'joined_drifting_paths':full,'joined_independent_steps':steps}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_aaa_certificates.json').read_text())
    result=verify(data);controls=0
    for kind in ('missing','order','capture'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['terminals'].pop()
        elif kind=='order':bad['terminals'][0]['orders'].pop()
        else:bad['terminals'][0]['capture_offsets'][0]+=1
        try:verify(bad)
        except (ArithmeticError,AssertionError):controls+=1
        else:raise ArithmeticError('Corrupt AAA proof accepted')
    result.update(actual_checks());result.update(actual_seam_checks());result.update({'status':'passed','corruptions_rejected':controls,
      'scope':'All-base first pair capture on the complete 435-state AA exit inventory and its preceding recurrent AA segment, including eight seam corrections; no original-input phase weights or arbitrary-phase theorem.'})
    (ROOT/'nine_aaa_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {result["strip_guards"]} strip guards, {result["incoming_transverse_guards"]} incoming transverse guards, '
            f'{result["edge_guards"]} edge guards and {result["terminal_guards"]} terminal guards.\n'
            'All 435 AA exit states reach the five proved phase cycles; 25 terminal families cover every base residue.\n'
            f'{result["compressed_captures"]} compressed captures; {len(result["independent_paths"])} independent full pair paths '
            f'over {result["independent_steps"]} subtraction steps; {controls} corruptions rejected.\n'
            f'{result["seam_comparisons"]} symbolic seam comparisons and {result["joined_stationary_captures"]} independent joined captures '
            f'prove all eight earlier-capture corrections; four full joined drifting paths add to {result["joined_independent_steps"]} joined independent steps.\n')
    (ROOT/'nine_aaa_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
