from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
from functools import lru_cache, reduce
import json
import sys
from faces import layout, permutation, spec_of
from work import w_phase, dot

def scale(c,row):
    return tuple(c*x for x in row)

def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))

def eq_guard(label,spec):
    middle,groups,endpoint=spec
    m=len(label)
    one=(1,)*m
    q=tuple(tuple(int(i in label[j:]) for i in range(m)) for j in range(1,m))
    radii=[[] for _ in groups]
    equal=[]
    for f,(j,s) in zip(q,layout(spec)):
        size=sub(scale(2,f),one)
        if s:
            radii[j].append(scale(s,size))
        else:
            equal.append(size)
    top=[r[0] for r in radii]
    equal.extend(sub(f,r[0]) for r in radii for f in r[1:])
    strict=[sub(a,b) for a,b in zip(top,top[1:])]
    if top:
        strict.append(top[-1])
        if endpoint:
            equal.append(sub(one,top[0]))
        else:
            strict.append(sub(one,top[0]))
    return equal,strict

def canonical(row,sign_free=False):
    d=reduce(gcd,map(abs,row),0)
    if not d:
        return row
    if sign_free and next(x for x in row if x)<0:
        d=-d
    return tuple(x//d for x in row)

def reduce_rows(rows,m):
    rows=[tuple(map(Q,r)) for r in rows]
    pivots=[]
    for matrixj in range(m):
        n=len(pivots)
        loc=next((i for i in range(n,len(rows)) if rows[i][matrixj]),None)
        if loc is None:
            continue
        rows[n],rows[loc]=rows[loc],rows[n]
        pivots.append(matrixj)
        rows[n]=scale(1/rows[n][matrixj],rows[n])
        for i,row in enumerate(rows):
            if i!=n and row[matrixj]:
                rows[i]=sub(row,scale(row[matrixj],rows[n]))
    return rows,pivots

def vertex_solver(eq,inequalities,m):
    # All constraints except normalization are homogeneous in gap space.
    equal,pivots=reduce_rows([r+(0,) for r in eq]+[(1,)*m+(1,)],m)
    assert all(not r[-1] for r in equal[len(pivots):])
    free=[i for i in range(m) if i not in pivots]
    shape=[None]*m
    for i in free:
        shape[i]=tuple(int(i==j) for j in free)+(0,)
    for i,row in zip(pivots,equal):
        shape[i]=tuple(-row[j] for j in free)+(row[-1],)
    def transform(constraint):
        row=tuple(sum(c*shape[i][j] for i,c in enumerate(constraint))
                  for j in range(len(free)+1))
        if any(row[:-1]):
            v=next(x for x in row if x)
            return scale(1/abs(Q(v)),row)
        return row
    base=tuple(sorted({transform(tuple(int(i==j) for i in range(m)))
                        for j in range(m)}))
    cuts=sorted({transform(r) for r in inequalities}-set(base))
    present=[]
    for group in combinations(base,len(free)):
        sol,piv=reduce_rows([r[:-1]+(-r[-1],) for r in group],len(free))
        if len(piv)<len(free):
            continue
        root=tuple(r[-1] for r in sol)+(1,)
        if all(dot(r,root)>=0 for r in base):
            present.append(root)
    active_rows=list(base)
    vertices={v:frozenset(i for i,r in enumerate(base) if dot(r,v)==0)
              for v in present}
    for cut in cuts:
        values={v:dot(cut,v) for v in vertices}
        n=len(active_rows)
        next_vertices={v:active | {n} if values[v]==0 else active
                       for v,active in vertices.items() if values[v]>=0}
        plus=[v for v,x in values.items() if x>0]
        minus=[v for v,x in values.items() if x<0]
        @lru_cache(maxsize=None)
        def edge(common):
            if len(common)<len(free)-1:
                return False
            _,piv=reduce_rows([active_rows[i][:-1] for i in sorted(common)],len(free))
            return len(piv)>=len(free)-1
        for a in plus:
            for b in minus:
                common=vertices[a] & vertices[b]
                if edge(common):
                    t=values[a]/(values[a]-values[b])
                    point=tuple(x+t*(y-x) for x,y in zip(a,b))
                    next_vertices[point]=common | {n}
        vertices=next_vertices
        active_rows.append(cut)
    return sorted(tuple(dot(r,v) for r in shape) for v in vertices)


def family(orbit,k):
    m=k+1
    label=tuple(range(m))
    equal=[]
    strict=[]
    labels=[]
    for spec in orbit:
        labels.append(label)
        eq,guards=eq_guard(label,spec)
        equal.extend(eq)
        strict.extend(guards)
        p=permutation(k,spec)
        label=tuple(label[i] for i in p)
    equal.extend(tuple(int(i==a)-int(i==b) for i in range(m)) for a,b in enumerate(label))
    return labels,equal,strict


def anchor_cut(m):
    return (1-m,)+(1,)*(m-1)


def normal(g):
    return max(abs(x-Q(1,len(g))) for x in g)


def work_orbit(g,orbit,labels):
    phases=[tuple(g[i] for i in l) for l in labels]
    return sum(w_phase(a,b,spec) for a,b,spec in zip(phases,phases[1:]+phases[:1],orbit))


def run(records,k,ball_words=None):
    m=k+1
    rows=[]
    winner=None
    hcounts=Counter()
    exceptional={}
    det_cache={}
    ball_ratios=[]
    if ball_words is not None:
        records=sorted(records,key=lambda row: json.dumps(row['faces']) not in ball_words)
    for r in records:
        orbit=tuple((a,tuple(map(tuple,b)),c) for a,b,c in r['faces'])
        if orbit!=min(orbit[i:]+orbit[:i] for i in range(len(orbit))):
            continue
        h=len(orbit)
        labels,eq,strict=family(orbit,k)
        eq=tuple(sorted({canonical(f,True) for f in eq if any(f)}))
        inequalities=tuple(sorted({canonical(f) for f in strict}
                                  | {tuple(int(i==j) for j in range(m)) for i in range(m)}))
        def vertices_for(extra):
            key=(eq,inequalities+(extra,))
            if key not in det_cache:
                det_cache[key]=vertex_solver(*key,m)
            return det_cache[key]
        cut=anchor_cut(m)
        verts=vertices_for(cut)
        if verts and max(v[0] for v in verts)>0 and min(v[0] for v in verts)<Q(1,m):
            # The open family intersects the strict anchor slab, since all
            # these vertices are closure points of a nonempty open face.
            hcounts[h]+=1
            values=[]
            for v in verts:
                w=work_orbit(v,orbit,labels)
                n=normal(v)
                if n:
                    values.append((w/(h*n),v,w))
            z=max(values)
            if winner is None or z[0]>winner[0]:
                winner=(*z,orbit)
                print('new upper ratio',winner,flush=True)
            if ball_words is not None and json.dumps(orbit) in ball_words:
                ball_ratios.append(z[0])
            rows.append(dict(faces=orbit,ratio=str(z[0]),point=list(map(str,z[1])),
                             work=str(z[2]),vertices=len(verts)))
        verts=vertices_for(scale(-1,cut))
        if verts and max(v[0] for v in verts)>=Q(1,m):
            # Vertices on a strict guard need not lie in this family.
            centroid=tuple(sum(v[i] for v in verts)/len(verts) for i in range(m))
            if all(dot(f,centroid)>0 for f in strict):
                exceptional[orbit]=verts
    if ball_words is not None:
        assert len(ball_ratios)==185 and max(ball_ratios)==Q(-4,7)
        print('in-ball W cells',len(ball_ratios),'certified ratio',max(ball_ratios),flush=True)
    print('ordinary counts' ,hcounts,'max ratio',winner,'exceptional counts',Counter(map(len,exceptional)),flush=True)
    for orbit,vertices in exceptional.items():
        print('exceptional',orbit,vertices,flush=True)
    return rows

if __name__=='__main__':
    from pathlib import Path
    rec=json.loads(Path(sys.argv[2]).read_text())
    ball_words=None
    if len(sys.argv)>4:
        energy_rows=json.loads(Path(sys.argv[4]).read_text())
        ball_words={json.dumps(row['faces']) for row in energy_rows if row['in_ball']}
    rows=run(rec,int(sys.argv[1]),ball_words)
    Path(sys.argv[3]).write_text(json.dumps(rows,indent=2)+'\n')
