"""All-base first pair capture on the complete finite AA exit inventory.

The input may already be cyclic. Original-number arrival therefore requires
the caller's transient-prefix proof; no unconditional h+1 conversion is made.
"""
import json
from pathlib import Path
from functools import lru_cache
from fractions import Fraction as Q
from nine_entry import evaluate
from odd_affine import word

THRESHOLDS=((-2,-2,27,30),(-2,0,27,30),(-2,2,27,30),(-2,4,31,30),
            (-1,-2,2,5),(-1,4,6,5),(0,-2,27,30),(0,0,27,30),(0,2,27,30),
            (0,4,31,30),(2,-2,29,30),(2,0,29,30),(2,2,29,30),(2,4,31,30),
            (3,2,5,5),(4,2,31,30),(4,4,31,30))
INCOMING=tuple((-p,-p-r,e) for p,r,m,d in THRESHOLDS for e in range(m-d,m))
CANONICAL={32:(-3,-5,1),37:(-2,-2,-1)}
DECREMENTS={32:2,37:3}
# For these eight incoming phases, first capture occurs inside the final
# preceding AA word. Offsets are relative to (U,2U-B,3U-2B,B-U) at its exit.
EARLY_CAPTURE={
    (2,4,-1):(5,(1,0,3,-2)),(2,2,1):(2,(2,2,2,1)),
    (0,-2,-1):(5,(1,0,1,2)),(0,-2,1):(2,(0,-2,-4,-1)),
    (-2,0,-1):(5,(1,0,-1,-2)),(-2,-2,3):(2,(0,-2,-2,3)),
    (-2,-4,-1):(5,(1,0,-1,2)),(-2,-4,3):(2,(-2,-6,-6,1)),
}
ORDERS={
    'low':(4,5,0,1,6,2,7,8,3),
    'middle':(4,0,5,1,6,2,7,8,3),
    'high':(4,0,5,1,6,2,7,3,8),
}


def encode(B,U,phase):
    a,b,c=phase
    return U,2*U-B+a,3*U-2*B+b,B-U+c


def decode(B,q):
    u,v,w,z=q
    return u,(v-2*u+B,w-3*u+2*B,z+u-B)


def branch(phase):
    c=phase[2]
    return 'low' if c<=-2 else 'middle' if c<=0 else 'high'


def phase_step(phase):
    a,b,c=phase;m=min(c,1)
    return (2*m+max(-c-1,0)-max(c-1,0),
            3*m+min(0,-c-1)+a+1,a-b-m),-m


@lru_cache(None)
def phase_path(initial):
    """Finite, exhaustively verified incoming graph; no all-phase claim."""
    if initial not in INCOMING:raise ValueError('Outside the complete AA exit inventory.')
    phase=initial;seen={};states=[];shifts=[0]
    while phase not in seen:
        seen[phase]=len(states);states.append(phase);phase,du=phase_step(phase);shifts.append(shifts[-1]+du)
        if len(states)>1000:raise ArithmeticError('Finite AA exit graph failed to close.')
    h=seen[phase];period=len(states)-h;drift=shifts[-1]-shifts[h]
    if drift==0:return {'states':tuple(states),'shifts':tuple(shifts),'prefix':h,'period':period,'stationary':True}
    if period not in CANONICAL or drift!=-DECREMENTS[period]:raise ArithmeticError('Unknown drifting terminal phase.')
    index=states.index(CANONICAL[period],h)
    return {'states':tuple(states),'shifts':tuple(shifts),'prefix':index,'period':period,'stationary':False}


def edge_orders(period):
    p=CANONICAL[period];orders=[]
    for _ in range(period):orders.append(ORDERS[branch(p)]);p,_=phase_step(p)
    if p!=CANONICAL[period]:raise ArithmeticError('Canonical terminal phase did not close.')
    return tuple(orders)


@lru_cache(None)
def terminal_table():
    records=json.loads(Path(__file__).with_name('nine_aaa_terminals.json').read_text())
    return {(r['edge_period'],r['base_residue'],r['offset']):r for r in records}


def exact_pair_arrival(B,q):
    q=tuple(q)
    if not isinstance(B,int) or B<1000000 or any(not isinstance(x,int) for x in q):
        raise ValueError('Require integer differences and B>=1000000.')
    U,p=decode(B,q);path=phase_path(p);prefix=path['prefix'];period=path['period']
    # Check every actual prefix/finite cycle word; phase-only discovery is not
    # a substitute for guards at the caller's finite base.
    for i,phase in enumerate(path['states']):
        state=encode(B,U+path['shifts'][i],phase);rec=word(B,state,1)
        if rec['orders']!=(ORDERS[branch(phase)],):raise ValueError('Outside the AAA strip guards.')
        target,du=phase_step(phase)
        if tuple(evaluate(row,(*state,B,1)) for row in rec['matrix'])!=encode(B,U+path['shifts'][i]+du,target):
            raise ArithmeticError('AAA phase formula failed.')
    if path['stationary']:
        state=encode(B,U+path['shifts'][prefix],path['states'][prefix])
        return {'pair_time':prefix,'period':period,'first_cycle_state':state,
                'phase_period':period,'prefix':prefix,'edge_repetitions':0}
    V=U+path['shifts'][prefix];b,r=divmod(B,5);d=DECREMENTS[period]
    if not Q(4*B,5)+252<=V<=Q(5*B,6)-256:raise ValueError('Outside the certified terminal-edge strip.')
    count=(V-4*b-256)//d
    if count<0:raise ValueError('Missing edge buffer.')
    e=V-4*b-d*count;terminal=terminal_table()[period,r,e]
    # Every prefix first coordinate and every edge interior lies above the
    # eventual terminal cycle. Edge offsets are checked in the proof replay.
    if min(U+s for s in path['shifts'][:prefix+1])<=4*b+terminal['cycle_u_max']:
        raise ValueError('Missing prefix exclusion from the eventual cycle.')
    h=prefix+period*count+terminal['arrival']
    captured=tuple(k*b+c for k,c in zip((4,3,2,1),terminal['capture_offsets']))
    return {'pair_time':h,'period':terminal['period'],'first_cycle_state':captured,
            'phase_period':period,'prefix':prefix,'edge_repetitions':count,
            'terminal_offset':e,'base_residue':r}


def from_recurrent_aa(B,q):
    """Exact first pair capture from a recurrent AA word, including its seam.

    Finite-base guard failures are rejected. On every compact relative AA
    interior, with a recurrent phase, this succeeds eventually through all
    bases. It makes no assertion about an earlier caller-owned prefix.
    """
    import nine_aa as aa
    q=tuple(q);out=aa.corner_exit(B,q);target=out['target'];U,p=decode(B,target)
    tail=exact_pair_arrival(B,target);h=out['elapsed']+tail['pair_time'];captured=tail['first_cycle_state']
    correction=0
    if p in EARLY_CAPTURE:
        correction,offset=EARLY_CAPTURE[p]
        if tail['pair_time']!=0:raise ArithmeticError('Early capture must join an already cyclic exit.')
        h-=correction;base=encode(B,U,(0,0,0));captured=tuple(a+b for a,b in zip(base,offset))
    elif tail['phase_period'] in (32,37):
        # Within each repeated AA word its first coordinate decreases by a
        # positive fixed amount between repeats. The final word is therefore
        # enough to bound the whole preceding segment.
        orders,delta=aa.return_spec(aa.decode(B,q)[2]);start=tuple(a-d for a,d in zip(target,delta))
        rec=word(B,start,len(orders));minimum=min(evaluate(m[0],(*start,B,1)) for m in rec['offsets'])
        r=tail['base_residue'];e=tail['terminal_offset'];record=terminal_table()[tail['phase_period'],r,e]
        if minimum<=4*(B//5)+record['cycle_u_max']:raise ValueError('Missing AA-segment terminal exclusion.')
    return {'pair_time':h,'actual_time':h+1,'period':tail['period'],'first_cycle_state':captured,
            'corner':out,'tail':tail,'seam_correction':correction}
