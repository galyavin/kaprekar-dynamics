"""Exact finite A1=A2=B1 transition from every recurrent AB exit.

The 384 reached phase states are a complete geometric exit inventory, not
equiprobable input events. Both outgoing routes have all-base terminal
closure; forty-two-step outputs compose with nine_aaa.from_recurrent_aa.
"""
import json
from functools import lru_cache
from pathlib import Path
import nine_ab as ab
import nine_aa as aa
from nine_entry import evaluate
from odd_affine import word
from odd_gap import subtract_pair


@lru_cache(None)
def transitions():
    records=json.loads(Path(__file__).with_name('nine_aab_certificates.json').read_text())['transitions']
    return {tuple(r['initial']):r for r in records}


def decode(B,q):
    u,v,w,z=q
    return u,(2*u-v-B,u+v-w-B,u+v-2*z-B)


def encode(B,U,phase):
    p,r,s=phase
    if (U-p-s)%2:raise ValueError('AAB chart requires U=p+s modulo two.')
    return U,2*U-B-p,3*U-2*B-p-r,(3*U-2*B-p-s)//2


def transfer(B,q):
    q=tuple(q);_,phase=decode(B,q);record=transitions().get(phase)
    if record is None:raise ValueError('Outside the complete AB exit inventory.')
    orders=tuple(map(tuple,record['orders']));rec=word(B,q,len(orders))
    if rec['orders']!=orders:raise ValueError('Outside the AAB transition guards.')
    target=tuple(evaluate(row,(*q,B,1)) for row in rec['matrix'])
    orders,_=aa.return_spec(aa.decode(B,target)[2]);test=word(B,target,len(orders))
    if test['orders']!=orders or len(orders)!=record['return_length']:
        raise ValueError('Missing outgoing recurrent AA word.')
    return {'elapsed':len(record['orders']),'target':target,'return_length':len(orders)}


def from_recurrent_ab(B,q):
    """Compose the full AB segment and its bounded AAB transition."""
    out=ab.corner_exit(B,q);bridge=transfer(B,out['target'])
    return {'elapsed':out['elapsed']+bridge['elapsed'],'target':bridge['target'],
            'return_length':bridge['return_length'],'corner':out,'bridge':bridge}


@lru_cache(None)
def nine_terminals():
    records=json.loads(Path(__file__).with_name('nine_aab_terminals.json').read_text())
    return {(r['phase_s'],r['base_residue']):r for r in records}


def nine_pair_arrival(B,q):
    """First capture on the reached nine-step AAB edges, in every base."""
    from fractions import Fraction as Q
    q=tuple(q);prefix=[]
    for _ in range(3):
        U,(p,r,s)=decode(B,q)
        if (p,r)==(-1,-2):break
        if (p,r) not in aa.PHASE_CYCLES[1]:raise ValueError('Outside the recurrent nine-step phase cycle.')
        for _ in range(3):prefix.append(q);q=subtract_pair(B,q)
    else:raise ValueError('Canonical nine-step phase was not reached.')
    if s not in (-2,-1,0,1) or not Q(4*B,5)+252<=U<=Q(5*B,6)-256 or B<1000000:
        raise ValueError('Outside the reached, buffered AAB edge.')
    b,r=divmod(B,5);record=nine_terminals()[s,r];e=record['offset'];count=(U-4*b-e)//2
    if count<0 or 4*b+e+2*count!=U:raise ValueError('AAB edge parity or buffer failed.')
    if min([q[0],*[x[0] for x in prefix]])<=4*b+record['cycle_u_max']:
        raise ValueError('Missing nine-step prefix exclusion.')
    return {'pair_time':len(prefix)+9*count+record['arrival'],'period':record['period'],
            'first_cycle_state':tuple(k*b+c for k,c in zip((4,3,2,1),record['capture_offsets'])),
            'phase_s':s,'prefix':len(prefix),'edge_repetitions':count}


def complete_from_recurrent_ab(B,q):
    """Complete first capture after either AB-to-AAB branch, with guards."""
    from nine_aaa import from_recurrent_aa
    result=from_recurrent_ab(B,q);target=result['target']
    tail=from_recurrent_aa(B,target) if result['return_length']==42 else nine_pair_arrival(B,target)
    # Verify exclusion for the earlier AB word and the finite bridge. The
    # first coordinate decreases by 12 between successive AB returns.
    end=result['corner']['target'];start=tuple(a-d for a,d in zip(end,ab.DELTA_AB));rec=word(B,start,24)
    minimum=min(evaluate(m[0],(*start,B,1)) for m in rec['offsets'])
    bridge=word(B,end,result['bridge']['elapsed'])
    minimum=min([minimum,*[evaluate(m[0],(*end,B,1)) for m in bridge['offsets']]])
    cyc=tail['first_cycle_state'];umax=cyc[0]
    for _ in range(tail['period']):umax=max(umax,cyc[0]);cyc=subtract_pair(B,cyc)
    if cyc!=tail['first_cycle_state'] or minimum<=umax:raise ValueError('Missing AB prefix first-capture exclusion.')
    h=result['elapsed']+tail['pair_time']
    return {'pair_time':h,'actual_time':h+1,'period':tail['period'],
            'first_cycle_state':tail['first_cycle_state'],'transport':result,'tail':tail}
