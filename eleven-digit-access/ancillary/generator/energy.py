"""Exact range of centered squared-gap energy on each periodic family."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path

from geometry import as_gap, evaluate, leading
from lp import dot, family, reduce_rows, sub

SPHERE = Q(1,30)


def energy(g):
    return dot(g,g)-Q(1,len(g))


def affine_projection(vertices):
    origin,*rest=vertices
    basis=[]
    m=len(origin)
    for point in rest:
        d=sub(point,origin)
        rows,pivots=reduce_rows(basis+[d],m)
        if len(pivots)>len(basis):
            basis=rows[:len(pivots)]
    if not basis:
        return origin
    center=(Q(1,m),)*m
    normal=sub(center,origin)
    gram=[tuple(dot(a,b) for b in basis)+(dot(a,normal),) for a in basis]
    rows,pivots=reduce_rows(gram,len(basis))
    assert len(pivots)==len(basis)
    return tuple(a+sum(row[-1]*d[i] for row,d in zip(rows,basis))
                 for i,a in enumerate(origin))


def extremizers(verts,halfspaces):
    # Every face of a bounded polytope is an intersection of halfspaces'
    # equality sets. The minimizer is the orthogonal projection onto the
    # affine hull of the unique face containing it in relative interior.
    universe=frozenset(range(len(verts)))
    facets={frozenset(i for i,v in enumerate(verts) if dot(f,v)==0)
            for f in halfspaces}-{frozenset(),universe}
    seen={universe}
    stack=[universe]
    minimum=None
    while stack:
        inds=stack.pop()
        p=affine_projection([verts[i] for i in sorted(inds)])
        if all(dot(f,p)>=0 for f in halfspaces):
            value=energy(p)
            if minimum is None or value<minimum[0]:
                minimum=value,p
        for facet in facets:
            section=inds & facet
            if section and section not in seen:
                seen.add(section)
                stack.append(section)
    assert minimum is not None
    maximum=max(map(energy,verts))
    return minimum,maximum,[v for v in verts if energy(v)==maximum]


def classify(k,scope,path):
    m=k+1
    records=evaluate(k,scope,path)
    result=[]
    for word,vertices,q,phase in records:
        _,equal,strict=family(word,k)
        if scope=='positive':
            strict.append((1,)+(0,)*k)
        halfspaces=list(strict)+[tuple(int(i==j) for i in range(m)) for j in range(m)]
        (lo,minimizer),hi,maximizers=extremizers(vertices,halfspaces)
        g=as_gap(q)
        assert lo<=energy(g)<=hi
        for p in [minimizer,*maximizers]:
            assert all(dot(f,p)==0 for f in equal)
            assert min(p)>=0
        in_family=lambda p:all(dot(f,p)>0 for f in strict)
        lo_in=in_family(minimizer)
        hi_in=any(in_family(p) for p in maximizers)
        high=hi>SPHERE or (hi==SPHERE and hi_in)
        place='zero-anchor' if scope=='zero' else ('boundary' if min(g)==0 else 'interior')
        row={'faces':word,'q_witness':list(map(str,q)),
             'return_length':len(word),'witness_period':len(phase),
             'location':place,'energy_inf':str(lo),'energy_sup':str(hi),
             'inf_attained':lo_in,'sup_attained':hi_in,
             'in_ball':lo<SPHERE,'meets_high_energy':high,
             'maximizing_gaps':[list(map(str,v)) for v in maximizers]}
        result.append(row)
        if high and place=='interior':
            print('high energy interior',row,flush=True)
    print('counts',Counter((x['location'],x['in_ball'],x['meets_high_energy']) for x in result))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('positive',type=Path)
    p.add_argument('zero',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    records=classify(5,'positive',a.positive)+classify(5,'zero',a.zero)
    records.sort(key=lambda r:(not r['in_ball'],r['location'],r['faces']))
    a.output.write_text(json.dumps(records,indent=2)+'\n')
