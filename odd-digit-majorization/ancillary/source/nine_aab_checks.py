"""Complete AB exit inventory and exact finite triple-corner transitions."""
import argparse,copy,json
from fractions import Fraction as Q
from math import ceil
from collections import Counter
from pathlib import Path
import nine_ab as ab
import nine_aa as aa
from nine_aab import (encode,decode,transfer,from_recurrent_ab,transitions,
    nine_terminals,complete_from_recurrent_ab)
from nine_entry import compose,evaluate
from nine_b_facet_checks import require,from_orders
from nine_phase_checks import replay,raw_step,step
from odd_affine import word,progression

ROOT=Path(__file__).resolve().parent
# Coordinates (U,e,p,s,B,1), where e=A1, p=A2 and s=B1.
EXIT_CHART=((1,0,0,0,0,0),(2,0,-1,0,-1,0),(3,-1,-1,0,-2,0),
            (Q(3,2),0,Q(-1,2),Q(-1,2),-1,0),
            (0,0,0,0,1,0),(0,0,0,0,0,1))


def nine_seed(s,r):
    e=256+(s-1)%2
    return e,((4,e),(3,2*e-r+1),(2,3*e-2*r+3),(1,(3*e-2*r+1-s)//2))


def nine_terminal_generate():
    from nine_aaa_checks import terminal_generate
    records=[]
    for s in (-2,-1,0,1):
        for r in range(5):
            e,matrix=nine_seed(s,r);record=terminal_generate(9,r,e,initial=matrix);record['phase_s']=s;records.append(record)
    return records


def incoming():
    out=[];thresholds=[];count=0
    for p,s in ab.PHASES:
        _,guards,_=replay(from_orders(ab.return_orders((p,s))));ts=[]
        for g in guards:
            row=compose(g,EXIT_CHART)
            if row[0] or row[4]:
                a,b=row[0],row[4]
                require(Q(4,5)*a+b>=0 and Q(5,6)*a+b>=0 and Q(49,60)*a+b>0,
                        'Unused AB exit guard has strict leading margin')
                continue
            c=row[2]*p+row[3]*s+row[5];e=row[1];count+=1
            if e>0:ts.append(ceil(-c/e))
            else:require(e==0 and c>=0,'Complete transverse AB exit guards')
        threshold=max(ts);thresholds.append((p,s,threshold))
        out.extend((p,e,s) for e in range(threshold-24,threshold))
    require(len(out)==len(set(out))==384,'Complete unique AB exit inventory')
    return out,thresholds,count


def chart(phase):
    p,r,s=phase
    return ((1,0,0),(2,-1,-p),(3,-2,-p-r),
            (Q(3,2),-1,Q(-p-s,2)),(0,1,0),(0,0,1))


def strip_guard(row):
    a,b,c=row
    if a or b:
        require(Q(4,5)*a+b>=0 and Q(5,6)*a+b>=0 and Q(49,60)*a+b>0,
                'Whole open AAB-line leading guard')
    else:require(c>=0,'Exact constant AAB guard')


def generate():
    initials,_,_=incoming();B=1000000000;records=[]
    for initial in initials:
        U=41*B//50;U+=(initial[0]+initial[2]-U)%2;q=encode(B,U,initial);prefix=[]
        for n in range(1000):
            phase=aa.decode(B,q)[2]
            if phase in aa.PHASES:
                orders,delta=aa.return_spec(phase);rec=word(B,q,len(orders))
                if rec['orders']==orders:
                    count,d=progression(B,q,rec)
                    if d==delta and count>10000:
                        records.append({'initial':initial,'orders':prefix,'return_length':len(orders)});break
            rec=word(B,q,3);prefix.extend(rec['orders']);q=tuple(evaluate(row,(*q,B,1)) for row in rec['matrix'])
        else:raise ArithmeticError('Unresolved finite triple-corner discovery')
    data={'transitions':records,'terminals':nine_terminal_generate()}
    (ROOT/'nine_aab_certificates.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    table=[{k:v for k,v in r.items() if k!='orders'} for r in data['terminals']]
    (ROOT/'nine_aab_terminals.json').write_text(json.dumps(table,indent=2)+'\n')
    transitions.cache_clear();nine_terminals.cache_clear();return data


def verify(data):
    initials,thresholds,exit_guards=incoming()
    require([tuple(r['initial']) for r in data['transitions']]==initials,'Complete incoming transition inventory')
    guards_count=0;max_prefix=0;counts=Counter();longitudinal=0;canonical_s=set()
    for record,initial in zip(data['transitions'],initials):
        M,guards,offsets=replay(from_orders(tuple(map(tuple,record['orders']))));source=chart(initial)
        for g in guards:strip_guard(compose(g,source));guards_count+=1
        target=tuple(compose(row,source) for row in M)
        require(target[0][:2]==(1,0) and Q(target[0][2]).denominator==1,'Integer longitudinal prefix shift')
        require(tuple(row[:2] for row in target)==tuple(row[:2] for row in source[:4]),
                'Unchanged leading AAB line and constant outgoing integer phases')
        phase=(2*target[0][2]-target[1][2],target[0][2]+target[1][2]-target[2][2])
        require(phase in aa.PHASES,'Outgoing AA recurrent phase')
        orders,delta=aa.return_spec(phase);R,gs,return_offsets=replay(from_orders(orders))
        require(len(orders)==record['return_length'],'Outgoing return length')
        extended=(*target,(0,1,0),(0,0,1))
        after=tuple(compose(row,extended) for row in R)
        require(after==tuple((a,b,c+d) for (a,b,c),d in zip(target,delta)), 'Exact outgoing displacement')
        for g in gs:
            row=compose(g,extended);strip_guard(row);guards_count+=1
            if row[:2]==(0,0):require(evaluate(g[:4],delta)>=0,'No bounded obstruction to repeated outgoing word')
        for off in offsets:
            first=compose(off[0],source);require(first[:2]==(1,0),'Bounded prefix longitudinal motion')
            longitudinal=max(longitudinal,abs(first[2]))
        if len(orders)==9:
            n=0;phase_now=phase;current=target
            while phase_now!=(-1,-2):
                M3,_,_=replay(from_orders(aa.three_orders(aa.branch(phase_now))))
                current=tuple(compose(row,(*current,(0,1,0),(0,0,1))) for row in M3)
                phase_now=(2*current[0][2]-current[1][2],current[0][2]+current[1][2]-current[2][2]);n+=3
                require(n<=6,'Canonical nine-step alignment')
            s=current[0][2]+current[1][2]-2*current[3][2]
            require(s in (-2,-1,0,1),'Reached nine-step transverse phase');canonical_s.add(s)
        max_prefix=max(max_prefix,len(record['orders']));counts[len(orders)]+=1
    require(dict(counts)=={42:344,9:40} and max_prefix==42,'Finite phase transition classification')
    require(canonical_s=={-2,-1,0,1},'Complete reached nine-step edge inventory')
    for p in ab.PHASES:
        _,_,offsets=replay(from_orders(ab.return_orders(p)))
        require(all(evaluate(off[0][:4],ab.DELTA_AB)==-12 for off in offsets),
                'First coordinate decreases by twelve between all AB word repetitions')
    result={'incoming_transverse_guards':exit_guards,'transition_guards':guards_count,
            'thresholds':thresholds,'incoming_states':384,'return_counts':dict(counts),
            'maximum_prefix':max_prefix,'maximum_longitudinal_offset':str(longitudinal)}
    result.update(terminal_checks(data['terminals']));return result


def terminal_checks(records):
    expected=[(s,r) for s in (-2,-1,0,1) for r in range(5)]
    require([(r['phase_s'],r['base_residue']) for r in records]==expected,'All twenty nine-step terminal families')
    edge_guards=0;terminal_guards=0;tables=[];edge_min={};errors=[];terminal_error=Q(0)
    for s in (-2,-1,0,1):
        source=chart((-1,-2,s));orders,delta=aa.return_spec((-1,-2));M,guards,offsets=replay(from_orders(orders))
        require(tuple(compose(row,source) for row in M)==tuple((a,b,c-2*a) for a,b,c in source[:4]),
                'Exact nine-step longitudinal return')
        for g in guards:
            a,b,c=compose(g,source);slope=b+a*(Q(4,5) if a>=0 else Q(5,6));constant=c+(252*a if a>=0 else -256*a)
            require(slope>=0 and 1000000*slope+constant>=0,'Unbounded buffered nine-step guard');edge_guards+=1
        us=[compose(off[0],source) for off in offsets]
        require(all(a==1 and b==0 for a,b,c in us),'Bounded nine-step first-coordinate offsets')
        edge_min[s]=min(c for a,b,c in us)
    for record,(s,r) in zip(records,expected):
        e,matrix=nine_seed(s,r);seen={};path=[]
        require(record['offset']==e and record['edge_period']==9,'Correct parity-compatible terminal seed')
        for order in record['orders']:
            require(matrix not in seen,'First terminal repeat');seen[matrix]=len(path);path.append(matrix)
            matrix,guards=step(matrix,(5,r),(0,1),order)
            for a,c in guards:
                require(a>=0 and 200000*a+c>=0,'All-base nine-step terminal guard');terminal_guards+=1
        require(matrix in seen,'Nine-step terminal closure');h=seen[matrix];period=len(path)-h
        require(all(tuple(a for a,c in m)==(4,3,2,1) for m in path),'No hidden earlier numeric repeat')
        capture=[c for a,c in path[h]];umax=max(m[0][1] for m in path[h:])
        require((h,period,capture,umax)==(record['arrival'],record['period'],record['capture_offsets'],record['cycle_u_max']),
                'Exact first capture and terminal cycle range')
        require(e+edge_min[s]>umax,'Entire earlier nine-step edge lies outside the terminal cycle')
        tables.append({k:v for k,v in record.items() if k!='orders'})
        errors.append(Q(h)-Q(9,2)*(Q(e)-Q(4*r,5)))
        terminal_error=max(terminal_error,*[abs(Q(c)-Q(a*r,5)) for a,c in path[h]])
    require(tables==json.loads((ROOT/'nine_aab_terminals.json').read_text()),'Runtime nine-step table matches proof')
    return {'nine_edge_guards':edge_guards,'nine_terminal_guards':terminal_guards,
            'nine_clock_correction_range':list(map(str,(min(errors),max(errors)))),
            'nine_first_terminal_error_constant':str(terminal_error)}


def actual_checks():
    steps=0;checks=0
    for phase in incoming()[0]:
        for B in range(1000000,1000005):
            U=41*B//50;U+=(phase[0]+phase[2]-U)%2;q=encode(B,U,phase);out=transfer(B,q)
            for _ in range(out['elapsed']):q=raw_step(B,q)
            require(q==out['target'],'Independent complete triple-corner transition')
            steps+=out['elapsed'];checks+=1
    joins=[];full={}
    for phase in ab.PHASES:
        for B in range(1000000,1000005):
            U=21*B//25;U+=(sum(phase)-U)%2;q=ab.encode(B,U,49*B//100,phase)
            out=complete_from_recurrent_ab(B,q);length=out['transport']['return_length']
            joins.append({'base':B,'phase':phase,'return_length':length,'pair_time':out['pair_time'],
                          'period':out['period'],'first_cycle_state':out['first_cycle_state']})
    # A single macroscopic center need not realize both phase branches in
    # every base. Choose an explicit reached exit of each branch and prepend
    # many actual AB returns; this retains their deterministic residues.
    representatives={}
    for initial,record in transitions().items():representatives.setdefault(record['return_length'],initial)
    for length,phase in representatives.items():
        for B in range(1000000,1000005):
            U=33*B//40;U+=(phase[0]+phase[2]-U)%2;end=encode(B,U,phase)
            q=tuple(a-(B//800)*d for a,d in zip(end,ab.DELTA_AB));out=complete_from_recurrent_ab(B,q)
            require(out['transport']['return_length']==length,'Preserved explicit incoming residue branch')
            full[length,B]=(B,q,out)
    from nine_aaa_checks import raw_capture
    complete=[];complete_steps=0
    require(len(full)==10,'Both outgoing branches in every base residue')
    for B,q,out in full.values():
        h,period,captured,n=raw_capture(B,q)
        require((h,period,captured)==(out['pair_time'],out['period'],out['first_cycle_state']),
                'Independent complete AB-to-terminal first capture')
        complete_steps+=n;complete.append({'base':B,'source':q,'pair_time':h,'period':period,
                                         'first_cycle_state':captured,'return_length':out['transport']['return_length']})
    return {'independent_transitions':checks,'independent_steps':steps,'joined_ab_checks':joins,
            'independent_complete_paths':complete,'independent_complete_steps':complete_steps}


def actual_number_checks(paths):
    chosen={};records=[];steps=0
    for path in paths:chosen.setdefault(path['return_length'],path)
    for path in chosen.values():
        B=path['base'];value=0
        for d in (*path['source'],0,0,0,0,0):value=B*value+d
        seen={}
        while value not in seen:
            seen[value]=len(seen);n=value;digits=[]
            for _ in range(9):n,d=divmod(n,B);digits.append(d)
            require(n==0,'Fixed nine-slot original numerical state');digits.sort();low=high=0
            for d in digits:low=B*low+d
            for d in reversed(digits):high=B*high+d
            value=high-low
        h=seen[value];period=len(seen)-h
        require((h,period)==(path['pair_time']+1,path['period']),'Actual original-number clock across the complete route')
        steps+=len(seen);records.append({'base':B,'source':path['source'],'return_length':path['return_length'],
                                       'actual_time':h,'period':period,'steps':len(seen)})
    require(len(records)==2,'Both complete route types checked numerically')
    return {'independent_number_paths':records,'independent_number_steps':steps}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    data=generate() if args.generate else json.loads((ROOT/'nine_aab_certificates.json').read_text())
    result=verify(data);controls=0
    for kind in ('missing','word','return','terminal'):
        bad=copy.deepcopy(data)
        if kind=='missing':bad['transitions'].pop()
        elif kind=='word':
            i=next(i for i,r in enumerate(bad['transitions']) if r['orders']);bad['transitions'][i]['orders'].pop()
        elif kind=='return':bad['transitions'][0]['return_length']=9
        else:bad['terminals'].pop()
        try:verify(bad)
        except (ArithmeticError,AssertionError,ValueError):controls+=1
        else:raise ArithmeticError('Corrupt triple-corner proof accepted')
    result.update(actual_checks());result.update(actual_number_checks(result['independent_complete_paths']))
    result.update({'status':'passed','corruptions_rejected':controls,
          'scope':'All 384 recurrent AB exits and both complete terminal routes, eventually on compact AB/AAB interiors through all bases; state counts are not input probabilities.'})
    (ROOT/'nine_aab_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f'PASS: {result["incoming_transverse_guards"]} incoming transverse guards and {result["transition_guards"]} transition guards.\n'
            'All 384 AB exits reach an outgoing word within 42 steps: 344 select length 42 and 40 select length 9.\n'
            f'{result["independent_transitions"]} independent complete transitions over {result["independent_steps"]} steps; '
            f'80 complete AB captures; {controls} corruptions rejected.\n'
            f'{result["nine_edge_guards"]} edge and {result["nine_terminal_guards"]} terminal guards cover all twenty nine-step regimes.\n'
            f'Ten complete independent captures: {result["independent_complete_steps"]} further subtraction steps.\n'
            f'Two original-number captures independently verify T=h+1 over {result["independent_number_steps"]} steps.\n')
    (ROOT/'nine_aab_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
