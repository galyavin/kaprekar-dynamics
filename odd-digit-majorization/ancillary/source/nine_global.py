"""Nine-digit first-core entry and transfer of the established joint kernel.

The evaluator has an explicit resource cap and rejects exceptional leading
boundaries. A cap or rejected finite-base guard is never a nonentry proof.
Global almost-everywhere validity is supplied separately by the entry audit
and analytical proof, not by this search routine.
"""
from fractions import Fraction as Q
from nine_entry import H,C,inverse_C,in_core_radius,evaluate


def leading_target(x,max_steps=200000):
    x=tuple(map(Q,x))
    if len(x)!=4 or not 1>x[0]>x[1]>x[2]>x[3]>0:
        raise ValueError('Require strict normalized original differences.')
    r=tuple(sorted((abs(2*v-1) for v in x),reverse=True))
    if min(r)<=0 or len(set(r))!=4:
        raise ValueError('Initial folded source lies on a leading boundary.')
    steps=0;seen=set()
    while not in_core_radius(r,strict=True):
        if in_core_radius(r):raise ValueError('Leading entry lies on a core boundary.')
        if r in seen:raise ValueError('Exceptional leading periodic orbit outside the strict core.')
        if steps>=max_steps:raise ValueError('Leading evaluation cap; no nonentry conclusion.')
        seen.add(r)
        raw=(r[0],*(abs(a+b-1) for a,b in zip(r,r[1:])))
        if min(raw)<=0 or len(set(raw))!=4:
            raise ValueError('Leading prefix lies on a sorting or absolute-value boundary.')
        r=H(r);steps+=1
    return C(r),1+steps


def limit_kernel(x,residue,max_steps=200000):
    from nine_collar import route_table
    if residue not in range(5):raise ValueError('Base residue must be 0,...,4.')
    y,prefix=leading_target(x,max_steps);r=inverse_C(y);domains,table=route_table()
    routes=[tuple(s['route']) for s in domains if min(evaluate(g,(*r,1)) for g in s['domain'])>0]
    if len(routes)!=1:raise ValueError('Leading core route lies on an exit tie.')
    record=table[routes[0],residue]
    atoms=[{'weight':Q(atom['weight']),'time':evaluate(tuple(map(Q,atom['clock'])),(*y,1)),
            'terminal':tuple(evaluate(tuple(map(Q,g)),(*y,1)) for g in atom['terminal'])}
           for atom in record['atoms']]
    return {'source':tuple(map(Q,x)),'core_source':y,'prefix_steps':prefix,
            'route':routes[0],'atoms':atoms}


def exact_arrival(B,q,max_steps=200000):
    from odd_gap import subtract_pair
    from nine_core import exact_arrival as core_arrival
    q=tuple(q)
    if not isinstance(B,int) or B<1000000 or any(not isinstance(v,int) for v in q):
        raise ValueError('Require integer differences and B>=1000000.')
    _,prefix=leading_target(tuple(Q(v,B) for v in q),max_steps);target=q
    for _ in range(prefix):target=subtract_pair(B,target)
    tail=core_arrival(B,target,first_image=False)
    if tail['pair_time']<=0:raise ArithmeticError('Core target must be transient.')
    return {'source':q,'prefix_steps':prefix,'core_source':target,'core_result':tail,
            'pair_time':prefix+tail['pair_time'],'actual_time':prefix+tail['actual_time'],
            'period':tail['period'],'first_cycle_state':tail['first_cycle_state'],
            'first_cyclic_number_differences':subtract_pair(B,tail['first_cycle_state']),
            'route':tail['route']}
