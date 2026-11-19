from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
from lp import family,vertex_solver,canonical,anchor_cut,dot
from faces import spec_of
from boundary_all import ORBIT_DATA

def leading(q):
    k=len(q)
    digits=sorted((Q(1),)+q+tuple(1-x for x in q),reverse=True)
    return tuple(digits[i]-digits[-1-i] for i in range(k))


def as_gap(q):
    return tuple(a-b for a,b in zip((Q(1),)+q,q+(Q(0),)))


def rotations(orbit):
    return tuple(orbit[i:]+orbit[:i] for i in range(len(orbit)))


def two_cone(g):
    a,t,u,v,w,z=g
    return (t==z==a and
            ((u==a and v>=a and 3*a+2*v<=1) or
             (v==a and w>=a and 3*a+2*w<=1)))


def exceptional_points():
    return {as_gap(tuple(Q(x,d) for x in phase))
            for name,d,orbit,_,_ in ORBIT_DATA[:2] for phase in orbit}


def verify_boundary_and_x(records,k):
    if k!=5:
        return
    special=exceptional_points()
    boundary_specs=set()
    known={}
    for key,d,orbit,_,_ in ORBIT_DATA:
        normalized=tuple(tuple(Q(x,d) for x in q) for q in orbit)
        for word in rotations(normalized):
            known[tuple(spec_of(q) for q in word)]=(key,word)
    top=((Q(4,5),Q(3,5),Q(2,5),Q(1,5),Q(1,5)),
         (Q(4,5),Q(3,5),Q(3,5),Q(2,5),Q(1,5)))
    for word in rotations(top):
        known[tuple(spec_of(q) for q in word)]=('two-top',word)
    xcounts=Counter()
    for word,verts,q,phase in records:
        if min(phase[0])==0:
            leading_word=tuple(spec_of(tuple(sum(g[i:]) for i in range(1,6)))
                               for g in phase)
            name,witness=known[leading_word]
            boundary_specs.add(name)
            assert min(phase[0])==0
            zero_slots=[i for i,g in enumerate(phase[0]) if not g]
            assert all(all(v[i]==0 for i in zero_slots) for v in verts)
            if name!='three-open':
                assert verts==[as_gap(witness[0])]
            else:
                assert len(verts)==2
                # Check the closure against the two explicitly checked endpoint orbits.
                endpoint_fibers=[]
                for label,den,orbit,_,_ in ORBIT_DATA[3:]:
                    endpoint_fibers.append({as_gap(tuple(Q(x,den) for x in q)) for q in orbit})
                assert all(len(set(verts)&p)==1 for p in endpoint_fibers)
        h=len(word)
        labels,eq,strict=family(word,k)
        inequalities=tuple({canonical(f) for f in strict}
                           | {tuple(int(i==j) for j in range(6)) for i in range(6)})
        a=tuple(-c for c in anchor_cut(6))
        xvert=vertex_solver(eq,inequalities+(a,),6)
        if not xvert:
            continue
        mid=tuple(sum(v[i] for v in xvert)/len(xvert) for i in range(6))
        if not all(dot(f,mid)>0 for f in strict):
            continue
        for lab in labels:
            phases=tuple(tuple(g[i] for i in lab) for g in xvert)
            if len(phases)==1 and phases[0] in special:
                continue
            # Each phase lies in one convex orientation of the closed Two triangle.
            assert (all(two_cone(g) and g[2]==g[0] for g in phases)
                    or all(two_cone(g) and g[3]==g[0] for g in phases)),phases
        xcounts[h]+=1
    assert boundary_specs=={'fixed','two','three-open','three-left','three-right','two-top'}
    print('positive-anchor boundary cells',sorted(boundary_specs),
          'exceptional intersections',sum(xcounts.values()))


def evaluate(k,scope,file):
    m=k+1
    families=json.loads(Path(file).read_text())
    records=[]
    periodic_counts=Counter()
    for r in families:
        orbit=tuple((a,tuple(map(tuple,b)),c) for a,b,c in r['faces'])
        if orbit!=min(orbit[i:]+orbit[:i] for i in range(len(orbit))):
            continue
        q=tuple(map(Q,r['q']))
        g=tuple(a-b for a,b in zip((Q(1),)+q,q+(Q(0),)))
        assert spec_of(q)==orbit[0]
        labels,eq,strict=family(orbit,k)
        if scope=='positive':
            strict.append(tuple(int(i==0) for i in range(m)))
        else:
            eq.append(tuple(int(i==0) for i in range(m)))
        eq=tuple(sorted({canonical(f,True) for f in eq if any(f)}))
        inequalities=tuple(sorted({canonical(f) for f in strict}
                                  | {tuple(int(i==j) for j in range(m)) for i in range(m)}))
        verts=vertex_solver(eq,inequalities,m)
        assert all(dot(f,g)==0 for f in eq)
        assert all(dot(f,g)>0 for f in strict)
        assert verts
        phase=tuple(tuple(g[i] for i in l) for l in labels)
        x=q
        for gphase, spec in zip(phase,orbit):
            assert spec_of(x)==spec and as_gap(x)==gphase
            x=leading(x)
        assert x==q
        per=next(i for i in range(1,len(phase)+1) if phase[i%len(phase)]==g)
        periodic_counts[per]+=1
        records.append((orbit,verts,q,phase[:per]))
    print('scope',scope,'families',len(records),'witness minimal period counts',sorted(periodic_counts.items()))
    if scope=='positive':
        verify_boundary_and_x(records,k)
        for orbit,verts,q,phase in records:
            if any(g==0 for g in phase[0]):
                print('boundary orbit', 'v=', len(verts), 'g=',phase[0],'phase=',phase)
    return records

if __name__=='__main__':
    k,scope,file=sys.argv[1:]
    evaluate(int(k),scope,file)
