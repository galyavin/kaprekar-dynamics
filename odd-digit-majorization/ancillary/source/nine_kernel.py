"""Exact lattice-weighted leading summaries on the 33 strict core routes.

Each branch retains an affine integer image lattice, not its saturation.
All sorting and stopping choices are checked on whole relative core faces.
The original-input counting argument is separate from this finite compiler.
"""
import json
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil
from pathlib import Path
import nine_entry as core
import nine_a_facet as af
import nine_b_facet as bf
import nine_aa as aa
import nine_ab as ab
import nine_bb as bb
import nine_aaa as aaa
from nine_aab import transitions
from nine_releases import records
from nine_releases_checks import face_vertices
from nine_b_routes import I,L,CM,flow,multiply
from nine_b_facet_checks import from_orders,require
from nine_phase_checks import replay
from nine_lattice import (identity,matvec,canonical,apply_affine,constant_form,
                          return_split)
from odd_affine import word

ROOT=Path(__file__).resolve().parent
ROWS={name:row for name,row,speed in core.EXIT_FACETS}
ZERO=(0,0,0,0,0)
ENDPOINT=tuple((0,0,0,0,Q(k,5)) for k in (4,3,2,1))
PHASE_ROWS={'A':(af.A2,),'B':(bf.B0,),'AA':(af.A2,aa.A1),
            'AB':(af.A2,bb.B1),'BB':(bf.B0,bb.B1)}
MODULES={'A':af,'B':bf,'AA':aa,'AB':ab,'BB':bb}
STATS=defaultdict(int)


@lru_cache(None)
def faces(active):
    vertices=face_vertices(active)
    require(bool(vertices),'Nonempty relative leading face')
    return vertices


@lru_cache(None)
def reconstructed(orders):return replay(from_orders(orders))


def point(T,center):return matvec(T,(*center,1))[:4]


def active_at(x):return tuple(g for g in bf.FACETS if core.evaluate(g,(*x,1))==0)


@lru_cache(None)
def active_matrix(T,center):
    x=point(T,center)
    require(min(core.evaluate(g,(*x,1)) for g in bf.FACETS)>=0,'Current leading point stays in the core')
    active=active_at(x)
    require(active==tuple(g for g in bf.FACETS if core.compose(g,T)==ZERO),
            'Face identities hold on the whole route, not only at a chosen center')
    return active


def guard_constant(state,g,P,r):
    return constant_form(state,g[:5],P,r)+g[5]


@lru_cache(None)
def guard_geometry(g,active):
    values=tuple(core.evaluate(g[:5],(*v,1)) for v in faces(active))
    require(min(values)>=0,('Negative homogeneous guard',g,active))
    return not any(values)


@lru_cache(None)
def prefix_leading(orders,active):
    M,_,_=reconstructed(orders);R=I
    for _ in range(len(orders)%3):R=multiply(L,R)
    for v in faces(active):
        require(matvec(tuple(row[:5] for row in M),(*v,1))==matvec(R,(*v,1))[:4],
                'Finite prefix agrees with its leading rotation on the whole face')


@lru_cache(None)
def check_prefix(state,orders,active,P,r):
    M,guards,_=reconstructed(orders)
    for g in guards:
        if guard_geometry(g,active):require(guard_constant(state,g,P,r)>=0,('Negative exact prefix guard',g,state))
        STATS['prefix_guards']+=1
    # Every finite prefix used here has the known leading core rotation.
    prefix_leading(orders,active)
    return apply_affine(state,M,P,r)


def prefix(branch,orders,center,P,r):
    state,T,clock=branch;orders=tuple(map(tuple,orders))
    state=check_prefix(state,orders,active_matrix(T,center),P,r)
    for _ in range(len(orders)%3):T=multiply(L,T)
    return state,T,clock


def near_point(state,x,P,r):
    """Choose a lattice point only to discover a rotation's sorting word."""
    b=10**10;B=P*b+r;a,c,basis=state;q=[v*b+w for v,w in zip(a,c)]
    for j in range(len(basis[0])-1,-1,-1):
        pivot=max(i for i,line in enumerate(basis) if line[j])
        n=round((B*x[pivot]-q[pivot])/basis[pivot][j])
        for i,line in enumerate(basis):q[i]+=n*line[j]
    require(max(abs(Q(v,B)-w) for v,w in zip(q,x))<Q(1,100000),
            'Rotation discovery point approximates its intended leading face')
    return B,tuple(q)


def rotate(branch,n,center,P,r):
    if not n:return branch
    B,q=near_point(branch[0],point(branch[1],center),P,r)
    orders=word(B,q,n)['orders']
    return prefix(branch,orders,center,P,r)


def phase(state,kind,P,r):
    values=tuple(constant_form(state,row,P,r) for row in PHASE_ROWS[kind])
    return values[0] if len(values)==1 else values


def absorb(branch,kind,center,P,r):
    module=MODULES[kind];seen=set()
    while (p:=phase(branch[0],kind,P,r)) not in module.PHASES:
        require(p not in seen,('Nonabsorbing incoming phase',kind,p));seen.add(p)
        if kind=='A':orders=af.THREE_ORDERS[af.branch(p)]
        elif kind=='B':orders=bf.THREE_ORDERS['high' if p>=-1 else 'low']
        elif kind=='BB':orders=bb.THREE_ORDERS[bb.branch(p)]
        else:orders=module.three_orders(module.branch(p))
        branch=prefix(branch,orders,center,P,r)
    return branch


@lru_cache(None)
def stopping_threshold(state,orders,active,selected,delta,P,r,strict=False):
    """Derive the exact last-word threshold from every sorting guard."""
    M,guards,_=(core.BULK_MATRIX,core.BULK_GUARDS,core.BULK_OFFSETS) if strict else reconstructed(orders)
    # Establish the full return matrix on this exact affine lattice.
    a,c,basis=state
    linear=tuple(g[:4] for g in M)
    require(multiply(linear,basis)==basis,'Repeated word fixes each incoming lattice vector')
    require(tuple(core.evaluate(g[:4],a)+g[4]*P for g in M)==a,
            'Repeated word fixes the incoming base coefficient')
    require(tuple(core.evaluate(g[:4],c)+g[4]*r+g[5] for g in M)==tuple(v+d for v,d in zip(c,delta)),
            'Repeated word is the stated pointwise translation')
    expected=canonical(a,tuple(v+d for v,d in zip(c,delta)),basis)
    require(apply_affine(state,M,P,r)==expected,'Repeated word translates the exact incoming lattice')
    # The matrix must fix every displacement, not just preserve its lattice.
    require(matvec(tuple(g[:4] for g in M),delta)==delta,'Repeated drift is invariant')
    speed=-core.evaluate(selected[:4],delta);require(speed>0,'Decreasing selected facet')
    fixed=tuple(g for g in active if core.evaluate(g[:4],delta)==0)
    outgoing=tuple(dict.fromkeys((*fixed,selected)))
    thresholds=[]
    h_a=core.evaluate(selected[:4],a)+selected[4]*P
    h_c=core.evaluate(selected[:4],c)+selected[4]*r
    h_basis=core.compose(selected[:4],basis)
    for g in guards:
        if guard_geometry(g,outgoing):
            slope=core.evaluate(g[:4],delta)
            if slope<0:
                lam=Q(-slope,speed)
                require(core.compose(g[:4],basis)==tuple(lam*v for v in h_basis),
                        'Stopping guard depends only on the selected remainder')
                require(core.evaluate(g[:4],a)+g[4]*P==lam*h_a,'Stopping guard has no hidden base slope')
                c0=core.evaluate(g[:4],c)+g[4]*r+g[5]-lam*h_c
                thresholds.append(ceil(-c0/lam))
            else:
                require(slope==0 and guard_constant(state,g,P,r)>=0,
                        ('Unaccounted exact stopping guard',g,slope))
        STATS['stopping_guards']+=1
    require(bool(thresholds),'A finite stopping threshold was found')
    return max(thresholds)


def repeated(branch,kind,selected,center,P,r):
    state,T,clock=branch;row=ROWS[selected];p=None
    if kind=='bulk':orders=(core.ORDER,)*3;delta=core.DELTA
    else:
        p=phase(state,kind,P,r)
        if kind=='A':orders=af.twelve_orders(p);delta=af.DELTA_A
        elif kind=='B':orders=bf.six_orders(p);delta=bf.DELTA_B
        elif kind=='BB':orders=bb.twelve_orders(p);delta=bb.DELTA_BB
        elif kind=='AB':orders=ab.return_orders(p);delta=ab.DELTA_AB
        elif kind=='AA':orders,delta=aa.return_spec(p)
        else:raise ArithmeticError('Unknown return kind')
    active=active_matrix(T,center)
    # Full prefix guards also prove the initial word is valid on the current
    # relative face, including sources released from a lower-dimensional one.
    check_prefix(state,orders,active,P,r)
    m=stopping_threshold(state,orders,active,row,delta,P,r,kind=='bulk')
    speed=-core.evaluate(row[:4],delta)
    Tnew=multiply(flow(row,speed,delta),T)
    clocknew=tuple(c+Q(len(orders),speed)*v for c,v in zip(clock,core.compose(row,T)))
    outputs=[]
    for e,w,new in return_split(state,row,speed,m,delta,P,r):
        outputs.append((w,(new,Tnew,clocknew)))
        STATS['weighted_transfers']+=1
    return outputs


def release(branch,kind,center,P,r):
    rows=(bf.B0,aa.A0) if kind=='BA0' else (bf.B0,bb.B1,aa.A1)
    p=tuple(constant_form(branch[0],g,P,r) for g in rows)
    record=records().get((kind,p));require(record is not None,'Complete reached release inventory')
    return prefix(branch,record['orders'],center,P,r)


def leaf(branch,kind):
    state,T,clock=branch
    if kind=='stationary':return clock,T[:4]
    rate={'32':Q(16),'37':Q(37,3),'9':Q(9,2),'BBB':Q(1)}[kind]
    distance=tuple(v-(Q(4,5) if i==4 else 0) for i,v in enumerate(T[0]))
    return tuple(c+rate*d for c,d in zip(clock,distance)),ENDPOINT


def compile_route(route,center,P,r):
    """Return exact (clock row, terminal rows) probabilities for one route."""
    source=((0,)*4,(0,)*4,identity(4));population={(source,I,ZERO):Q(1)}
    first=route[0];mode=first[0];turn=(int(first[1])+1)%3 if mode=='A' else int(first[1])

    def operation(pop,fn):
        out=defaultdict(Q)
        for b,w in pop.items():
            for p,new in fn(b):out[new]+=w*p
        require(sum(out.values())==sum(pop.values()),'Exact probability conservation')
        STATS['peak_population']=max(STATS['peak_population'],len(out))
        return dict(out)

    def deterministic(pop,fn):return operation(pop,lambda b:((Q(1),fn(b)),))

    population=operation(population,lambda b:repeated(b,'bulk',first,center,P,r))
    population=deterministic(population,lambda b:absorb(rotate(b,turn,center,P,r),mode,center,P,r))
    for stage in route[1:-1]:
        kind,selected=stage.split(':');require(kind==mode,'Consistent leading route')
        population=deterministic(population,lambda b:absorb(b,mode,center,P,r))
        population=operation(population,lambda b:repeated(b,kind,selected,center,P,r))
        if mode=='A':
            turns=1 if selected=='A0' else 0;mode='AB' if selected=='B1' else 'AA'
            population=deterministic(population,lambda b:rotate(b,turns,center,P,r))
        elif mode=='B' and selected=='A0':
            population=deterministic(population,lambda b:release(b,'BA0',center,P,r));mode='A'
        elif mode=='BB' and selected=='A1':
            population=deterministic(population,lambda b:release(b,'BBA1',center,P,r));mode='AB'
        else:
            turns=2 if mode=='B' and selected in ('A1','B2') else 0
            mode='BBB' if mode=='BB' else 'AB' if selected=='A1' else 'BB'
            population=deterministic(population,lambda b:rotate(b,turns,center,P,r))
    require(mode==route[-1],'Final leading route channel')
    atoms=defaultdict(Q)
    if mode=='BBB':
        for b,w in population.items():atoms[leaf(b,'BBB')]+=w
    else:
        population=deterministic(population,lambda b:absorb(b,mode,center,P,r))
        if mode=='AB':
            population=operation(population,lambda b:repeated(b,'AB','A1',center,P,r))
            onward=defaultdict(Q)
            for b,w in population.items():
                p=tuple(constant_form(b[0],row,P,r) for row in (af.A2,aa.A1,bb.B1))
                rec=transitions().get(p);require(rec is not None,'Complete AB exit inventory')
                b=prefix(b,rec['orders'],center,P,r)
                if rec['return_length']==9:atoms[leaf(b,'9')]+=w
                else:
                    require(rec['return_length']==42 and phase(b[0],'AA',P,r) in aa.PHASE_CYCLES[0],
                            'AAB transfer reaches the stated AA return')
                    onward[b]+=w
            population=dict(onward)
        population=operation(population,lambda b:repeated(b,'AA','A0',center,P,r))
        for b,w in population.items():
            p=tuple(constant_form(b[0],row,P,r) for row in ((-2,1,0,0,1),(-3,0,1,0,2),aa.A0))
            path=aaa.phase_path(p)
            kind='stationary' if path['stationary'] else str(path['period'])
            atoms[leaf(b,kind)]+=w
    require(sum(atoms.values())==1,'Complete route law has mass one')
    return dict(atoms)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--route',type=int);parser.add_argument('--residue',type=int)
    args=parser.parse_args();centers=json.loads((ROOT/'nine_core_test_points.json').read_text())['centers'];out=[]
    for i,rec in enumerate(centers):
        if args.route is not None and args.route!=i:continue
        center=core.C(tuple(map(Q,rec['radius'])))
        for r in range(5):
            if args.residue is not None and args.residue!=r:continue
            atoms=compile_route(tuple(rec['route']),center,5,r)
            out.append({'route':rec['route'],'base_residue':r,'atoms':[
                {'weight':str(w),'clock':list(map(str,key[0])),
                 'terminal':[list(map(str,row)) for row in key[1]]} for key,w in sorted(atoms.items())]})
            print(i,r,'/'.join(rec['route']),len(atoms),'atoms',dict(STATS),flush=True)
    if args.route is None and args.residue is None:
        (ROOT/'nine_kernel_results.json').write_text(json.dumps({'modulus':5,'kernels':out,'audit':dict(STATS)},indent=2)+'\n')


if __name__=='__main__':main()
