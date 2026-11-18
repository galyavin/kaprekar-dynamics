"""All-base guarded edge returns and finite terminal certificates.

The three-step affine return is verified for every offset in the reachable
finite graph, over unbounded integer bases and the stated longitudinal strip.
Terminal trajectories are checked as affine families in B modulo four.
"""
if not __debug__:
    raise RuntimeError(
        'Proof audits require assertions; disable -O/-OO and PYTHONOPTIMIZE.')

from fractions import Fraction as Q
import json
from pathlib import Path
import z3

from odd_gap import subtract_pair
from seven_boundary import offset_return, phase_data


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def scale(a, c):
    return tuple(c*x for x in a)


def affine_return(parity, state):
    m, s = state
    U, B, ONE = (1,0,0), (0,2,parity), (0,0,1)
    q = [U, (0,1,(parity+m)//2), (-1,2,parity+s)]
    guards = []
    for _ in range(3):
        guards.extend([sub(q[2],ONE), sub(q[0],q[1]), sub(q[1],q[2]),
                       sub(sub(B,ONE),q[0])])
        ds = [q[0], q[1], sub(q[2],ONE), sub(B,ONE),
              sub(sub(B,q[2]),ONE), sub(sub(B,q[1]),ONE), sub(B,q[0])]
        order = sorted(range(7), key=lambda i: (
            -sum(x*y for x,y in zip(ds[i], (80000,50000,1))), i))
        for i,j in zip(order,order[1:]):
            guard = sub(ds[i],ds[j])
            if i > j:
                guard = sub(guard,ONE)
            guards.append(guard)
        q = [sub(ds[order[i]],ds[order[-1-i]]) for i in range(3)]
        # Certify every intermediate longitudinal displacement on the entire strip.
        longitudinal = sub(q[0], U)
        assert longitudinal[:2] == (0,0) and abs(longitudinal[2]) <= 20
    du, mn, sn = sub(q[0],U), sub(scale(q[1],2),B), sub(add(q[0],q[2]),B)
    assert du[:2] == mn[:2] == sn[:2] == (0,0)
    result = (mn[2],sn[2]),du[2]
    assert result == offset_return(parity,state)
    return result, guards


def terminal_family(residue, state, decrement):
    base = 100000+residue
    threshold = (3*base+3)//4+256
    result = []
    for j in range(1,decrement+1):
        m,s = state; u=threshold-j
        start=(u,(base+m)//2,base-u+s)
        q, path, seen = start, [], {}
        while q not in seen:
            seen[q]=len(path);path.append(q)
            # For B -> B+4t and q -> q+(3,2,1)t, each borrowed digit
            # is affine. Its sorted order stays valid for all t>=0 when
            # consecutive slopes are nonincreasing; ties have equal value.
            a,b,c=q
            digits=[a,b,c-1,base-1,base-c-1,base-b-1,base-a]
            slopes=[3,2,1,4,3,2,1]
            order=sorted(range(7),key=lambda i:(-digits[i],i))
            assert all(slopes[i]>=slopes[k] for i,k in zip(order,order[1:]))
            assert [slopes[order[i]]-slopes[order[-1-i]] for i in range(3)] == [3,2,1]
            assert c>0
            q=subtract_pair(base,q)
            assert len(path)<5000
        h=seen[q]; cycle=path[h:]
        result.append({'u_above_center':str(Q(u)-Q(3*base,4)),
                       'tail':h,'period':len(cycle),
                       'cycle_max_u_above_center':str(max(Q(p[0])-Q(3*base,4) for p in cycle)),
                       'max_abs_center_margin':max(max(abs(4*p[0]-3*base),
                                                      abs(2*p[1]-base),
                                                      abs(4*p[2]-base)) for p in path)})
    return result


def audit():
    summaries, all_nodes, all_guards = [], {}, {}
    for parity in (0,1):
        initial=[(m,s) for m in range(-20,21) if m%2==parity
                 for s in range(-20,21)]
        pending=list(initial); nodes={}; guards=set()
        while pending:
            state=pending.pop()
            if state in nodes:
                continue
            result, local_guards=affine_return(parity,state)
            nodes[state]=result;guards.update(local_guards)
            if result[0] not in nodes:
                pending.append(result[0])
            assert len(nodes)<2000
        b,u=z3.Ints('b u');B=2*b+parity
        solver=z3.Solver()
        solver.add(B>=100000,4*u>=3*B+512,6*u<=5*B-768)
        solver.add(z3.Or(*[row[0]*u+row[1]*b+row[2]<0 for row in guards]))
        assert solver.check()==z3.unsat, ('affine edge guards',parity)
        cycles={};maximum_tail=0;micro_shift_bound=0
        for state in initial:
            seen={};path=[];shifts=[];q=state
            while q not in seen:
                seen[q]=len(path);path.append(q)
                q,d=nodes[q];shifts.append(d)
            cut=seen[q];cycle=path[cut:];j=cycle.index(min(cycle))
            cycle=tuple(cycle[j:]+cycle[:j]);delta=sum(shifts[cut:])
            cycles[cycle]=delta
            maximum_tail=max(maximum_tail,cut)
            micro_shift_bound=max(micro_shift_bound,
                                  max(abs(sum(shifts[:j])) for j in range(len(shifts)+1)))
        # Bound intermediate (single-step) first-coordinate excursions in
        # any certified offset return, not merely return endpoints.
        intermediate_bound=0
        for m,s in nodes:
            base=100000+parity;start=(80000,(base+m)//2,base-80000+s);q=start
            for _ in range(3):
                q=subtract_pair(base,q)
                intermediate_bound=max(intermediate_bound,abs(q[0]-start[0]))
        summaries.append({'parity':parity,'initial_states':len(initial),
                          'reachable_states':len(nodes),'guards':len(guards),
                          'unbounded_guard_result':'unsat',
                          'max_transient_blocks':maximum_tail,
                          'micro_cumulative_u_shift_bound':micro_shift_bound,
                          'single_return_intermediate_u_shift_bound':intermediate_bound,
                          'cycles':[{'offsets':c,'u_shift':d,'steps':3*len(c),
                                     'time_per_unit_u':str(Q(3*len(c),-d))}
                                    for c,d in sorted(cycles.items())]})
        all_nodes[parity]=nodes;all_guards[parity]=guards

    terminal_cache={};errors=[]
    for parity,nodes in all_nodes.items():
        for m in range(-20,21):
            if m%2!=parity:
                continue
            for s in range(-20,21):
                phase=phase_data(parity,(m,s));rate=phase['time_per_unit_u']
                assert rate in (Q(3,2),Q(2))
                for residue in range(parity,4,2):
                    key=(residue,phase['cycle_state'],-phase['cycle_u_shift'])
                    if key not in terminal_cache:
                        terminal_cache[key]=terminal_family(*key)
                    for result in terminal_cache[key]:
                        error=(3*phase['transient_blocks']+
                               rate*phase['transient_u_shift']+result['tail']-
                               rate*Q(result['u_above_center']))
                        errors.append(error)
                        assert Q(result['cycle_max_u_above_center'])<128
    assert max(map(abs,errors))<=200
    assert all(2*r['micro_cumulative_u_shift_bound']+
               r['single_return_intermediate_u_shift_bound']<128 for r in summaries)
    terminals=[{'base_residue_mod_4':key[0],'offsets':key[1],
                'decrement':key[2],'families':value} for key,value in sorted(terminal_cache.items())]
    output={'solver':z3.get_version_string(),'offset_graphs':summaries,
            'certified_domain':'B >= 100000, |2v-B|<=20, |u+w-B|<=20, '
                               '3B/4+1000 <= u <= 5B/6-1000',
            'pair_clock_error_range':[str(min(errors)),str(max(errors))],
            'stated_absolute_error_bound':200,
            'terminal_families':terminals,
            'scope':'Conditional edge arrival theorem, not weighted/global seven-digit entry or a population law.'}
    Path(__file__).with_name('seven_boundary_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Both offset graphs certified for all bases in scope; '
          f'{sum(len(x) for x in terminal_cache.values())} terminal families; '
          f'clock error [{min(errors)}, {max(errors)}].')


if __name__=='__main__':
    audit()
