"""Complete exact nine-digit capture on generic strict-core incoming routes.

All integer phases are retained. Finite-base guard failures and leading ties
are rejected; every compact generic route interior is eventually accepted.
This runtime is not a global leading-entry theorem or a phase-weight law.
"""
from fractions import Fraction as Q
import nine_entry as core
import nine_a_facet as af
import nine_b_facet as bf
import nine_aa as aa
import nine_ab as ab
import nine_bb as bb
import nine_bbb as bbb
from nine_aaa import from_recurrent_aa
from nine_aab import complete_from_recurrent_ab
from nine_releases import release
from odd_gap import subtract_pair,valid
from odd_affine import word

FACE_NAMES={row:name for name,row,speed in core.EXIT_FACETS}
ALLOWED={'A':{'A0','A1','B1'},'B':{'A0','A1','B1','B2'},'BB':{'A1','B2'}}


def leading_exit(x,kind):
    delta={'A':af.DELTA_A,'B':bf.DELTA_B,'BB':bb.DELTA_BB}[kind]
    candidates=[(core.evaluate(row,(*x,1))/(-core.evaluate(row[:4],delta)),row)
                for row in bf.FACETS if core.evaluate(row[:4],delta)<0]
    t=min(v for v,row in candidates);rows=[row for v,row in candidates if v==t]
    if t<=0 or len(rows)!=1 or FACE_NAMES.get(rows[0]) not in ALLOWED[kind]:
        raise ValueError('Leading route has a tie or a lower-dimensional exit.')
    return FACE_NAMES[rows[0]],tuple(v+t*d for v,d in zip(x,delta))


def exact_arrival(B,q,first_image=True):
    q=tuple(q)
    if not isinstance(B,int) or B<1000000 or len(q)!=4 or any(not isinstance(v,int) for v in q) or not valid(B,q):
        raise ValueError('Require valid integer nine-slot differences and B>=1000000.')
    initial=q;elapsed=0;minimum=None;stages=[];route=[]
    if first_image:
        if not B>q[0]>q[1]>q[2]>q[3]>0:raise ValueError('Require strict original differences.')
        radius=tuple(sorted((abs(Q(2*v,B)-1) for v in q),reverse=True));x=tuple(core.C(radius))
        q=subtract_pair(B,q);elapsed=1
    else:x=tuple(Q(v,B) for v in q)
    if not core.in_core(x,strict=True):raise ValueError('Leading source is outside the strict core.')

    def retain(v):
        nonlocal minimum
        minimum=v if minimum is None else min(minimum,v)

    def short(label,result):
        nonlocal q,elapsed
        start=q;walk=q
        for _ in range(result['elapsed']):retain(walk[0]);walk=subtract_pair(B,walk)
        if walk!=result['target']:raise ArithmeticError('Finite actual prefix disagrees with the certified transition.')
        q=walk;elapsed+=result['elapsed'];stages.append({'stage':label,'source':start,**result})

    def repeated(label,result,length,delta):
        nonlocal q,elapsed
        start=q;count=result['repetitions'];rec=word(B,q,length)
        actual=tuple(core.evaluate(row,(*q,B,1))-v for row,v in zip(rec['matrix'],q))
        if actual!=delta or tuple(core.evaluate(row[:4],delta) for row in rec['matrix'])!=delta:
            raise ArithmeticError('Repeated word has the wrong affine progression.')
        for off in rec['offsets']:
            v=core.evaluate(off[0],(*q,B,1));slope=core.evaluate(off[0][:4],delta)
            retain(min(v,v+(count-1)*slope))
        if result['target']!=tuple(v+count*d for v,d in zip(q,delta)) or result['elapsed']!=length*count:
            raise ArithmeticError('Incorrect compressed endpoint or elapsed time.')
        q=result['target'];elapsed+=result['elapsed'];stages.append({'stage':label,'source':start,**result})

    def rotate(n):
        nonlocal x
        for _ in range(n):short('rotation',{'elapsed':1,'target':subtract_pair(B,q)});x=core.leading_bulk(x)

    def absorb_facet(kind):
        module=af if kind=='A' else bf;start=q;target=q;n=0;e=module.phase(B,q)
        if abs(e)>512:raise ValueError('Outside the proved bounded first-facet phase inventory.')
        while module.phase(B,target) not in module.PHASES:
            e=module.phase(B,target);key=af.branch(e) if kind=='A' else 'high' if e>=-1 else 'low'
            if word(B,target,3)['orders']!=module.THREE_ORDERS[key]:raise ValueError('Facet absorption guards fail.')
            target=module.three_formula(B,target);n+=3
            if n>(390 if kind=='A' else 513):raise ArithmeticError('Proved facet absorption bound failed.')
        short(kind+' absorption',{'elapsed':n,'target':target})

    initial_exit=core.limiting_exit(x)
    if len(initial_exit['faces'])!=1:raise ValueError('Leading bulk exit is tied.')
    first=initial_exit['faces'][0];route.append(first);x=initial_exit['target']
    repeated('bulk',core.bulk_exit(B,q),3,core.DELTA)
    mode=first[0];rotate((int(first[1])+1)%3 if mode=='A' else int(first[1]));absorb_facet(mode)
    tail=None
    # The proved route tree is acyclic: B may enter BB, AB or A; A enters
    # AA or AB; BB enters BBB or AB. No numerical search cap replaces closure.
    while tail is None:
        if mode in ('A','B'):
            selected,target=leading_exit(x,mode);route.append(mode+':'+selected)
            module=af if mode=='A' else bf
            repeated(mode+' facet',module.facet_exit(B,q),12 if mode=='A' else 6,
                     af.DELTA_A if mode=='A' else bf.DELTA_B);x=target
            if mode=='A':
                if selected=='A0':rotate(1);mode='AA'
                elif selected=='A1':mode='AA'
                else:mode='AB'
            elif selected=='A0':
                short('BA0 release',release(B,q,'BA0'));x=core.leading_bulk(x);mode='A'
            elif selected=='A1':rotate(2);mode='AB'
            else:rotate(0 if selected=='B1' else 2);mode='BB'
        elif mode=='BB':
            p=bb.decode(B,q)[2];short('BB absorption',bb.absorb(B,q,max_returns=bb.phase_entry_bound(p)+1))
            selected,target=leading_exit(x,'BB');route.append('BB:'+selected)
            repeated('BB corner',bb.corner_exit(B,q),12,bb.DELTA_BB);x=target
            if selected=='B2':mode='BBB'
            else:
                short('BBA1 release',release(B,q,'BBA1'));x=core.leading_bulk(core.leading_bulk(x));mode='AB'
        elif mode=='AA':
            p=aa.decode(B,q)[2];short('AA absorption',aa.absorb(B,q,max_returns=aa.phase_entry_bound(p)+1))
            route.append('AA');tail=from_recurrent_aa(B,q)
        elif mode=='AB':
            p=ab.decode(B,q)[2];short('AB absorption',ab.absorb(B,q,max_returns=ab.phase_entry_bound(p)+1))
            route.append('AB');tail=complete_from_recurrent_ab(B,q)
        elif mode=='BBB':route.append('BBB');tail=bbb.exact_arrival(B,q)
        else:raise ArithmeticError('Unclassified reached region.')
    cycle=tail['first_cycle_state'];umax=cycle[0]
    for _ in range(tail['period']):umax=max(umax,cycle[0]);cycle=subtract_pair(B,cycle)
    if cycle!=tail['first_cycle_state'] or minimum is None or minimum<=umax:
        raise ValueError('Missing first-capture exclusion for an earlier route segment.')
    h=elapsed+tail['pair_time']
    return {'source':initial,'first_image':first_image,'pair_time':h,'actual_time':h+1,
            'period':tail['period'],'first_cycle_state':tail['first_cycle_state'],
            'route':tuple(route),'stages':stages,'tail':tail,'prefix_minimum':minimum,'cycle_maximum':umax}
