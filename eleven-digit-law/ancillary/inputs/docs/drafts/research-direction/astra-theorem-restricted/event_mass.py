"""Independently integrate the three-step row-181 event on its sixty strata.

Use Task 145's ordinary-difference Poly clipper, an independent inverse tree,
pulling triangulation and Dirichlet moments by polarization. No Task 130/141
inverse, geometry, triangulation or quadrature code is called.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import fcntl
import json
import shutil
import subprocess
import tempfile
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import lcm
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
ACCESS = ROOT/'kaprekar-investigation/eleven_access'
sys.path.insert(0, str(ACCESS/'restricted_theorem_check'))
from cover_exact import Poly

ZERO = (0,)*6
ONE = (0,0,0,0,0,1)
ID = tuple(tuple(int(i==j) for j in range(6)) for i in range(5))
MASKS = tuple(product((False,True), repeat=5))


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def scale(a,k):
    return tuple(x*k for x in a)


def gaps(q):
    return tuple(sub(a,b) for a,b in zip((ONE,*q),(*q,ZERO)))


def restrict(P, rows):
    for row in rows:
        P = P.clip(row)
        if P is None:
            return None
    return P


def low_chart(q):
    lows = []
    previous = ZERO
    for row in q:
        previous = sub(sub(ONE,row),previous)
        lows.append(previous)
    return lows


def paths(P,q,depth):
    lows = low_chart(q)
    cuts = (lows[0], *(sub(b,a) for a,b in zip(lows,lows[1:])), sub(ONE,scale(lows[-1],2)))
    P = restrict(P,cuts)
    if P is None:
        return
    for mask in MASKS:
        original = tuple([sub(ONE,v) for v,k in zip(lows,mask) if k]
                         +[v for v,k in reversed(list(zip(lows,mask))) if not k])
        if depth == 1:
            yield P,original
        else:
            yield from paths(P,original,depth-1)


def triangulate(P):
    facets = tuple(frozenset(f) for f in P.faces)
    @lru_cache(None)
    def pull(face, dim):
        if len(face) == dim+1:
            return (tuple(sorted(face)),)
        v = min(face)
        candidates = {f & face for f in facets if not face <= f} - {frozenset()}
        covered = [f for f in candidates if v not in f and all(not f<h for h in candidates)]
        assert covered
        return tuple((v,*s) for f in covered for s in pull(f,dim-1))
    return tuple(tuple(P.vertices[i] for i in s) for s in pull(frozenset(range(len(P.vertices))),5))


def det(A):
    a = [list(row) for row in A]
    sign, previous = 1,1
    for i in range(len(a)-1):
        j = next((k for k in range(i,len(a)) if a[k][i]),None)
        if j is None:
            return 0
        if j != i:
            a[j],a[i]=a[i],a[j]
            sign=-sign
        pivot=a[i][i]
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):
                n=a[j][k]*pivot-a[j][i]*a[i][k]
                a[j][k],r=divmod(n,previous)
                assert not r
        previous=pivot
    return sign*a[-1][-1]


@lru_cache(None)
def prepared(P):
    result = []
    for simplex in triangulate(P):
        scale = lcm(*(x.denominator for v in simplex for x in v))
        vertices = [tuple(int(scale*x) for x in v)+(scale,) for v in simplex]
        jac = abs(det(vertices))
        assert jac
        result.append((vertices,scale,jac))
    return result


@lru_cache(maxsize=100000)
def h6(args):
    coeff = [1,0,0,0,0,0,0]
    for a in args:
        for i in range(1,7):
            coeff[i] += a*coeff[i-1]
    return coeff[6]


def measure(P,q):
    # Jacobian is 1 for each inverse. Polarization of the six gap factors
    # cancels the factor 11! in the density against the Dirichlet denominator.
    assert abs(det((*q,ONE))) == 1
    rows = gaps(q)
    result = Q(0)
    for vertices,den,jac in prepared(P):
        values = [tuple(sum(a*b for a,b in zip(row,v)) for v in vertices) for row in rows]
        assert min(min(row) for row in values) >= 0
        sums = [((0,)*6,1)]
        for v in values:
            sums += [(tuple(x+y for x,y in zip(u,v)),-sign) for u,sign in sums]
        subtotal = sum(sign*h6(tuple(sorted(u))) for u,sign in sums)
        result += Q(jac*subtotal,den**12)
    return result


def check(S,largest):
    verts=tuple(tuple(Q(j<k) for j in range(5)) for k in range(6))
    P=Poly(5,verts,tuple(tuple(j for j in range(6) if j!=i) for i in range(6)))
    g=gaps(ID)
    rest=set(range(6))-set(S)
    assert largest in rest
    cuts=[sub(scale(g[i],8),ONE) for i in (0,5)]
    cuts.append(sub(scale(ONE,3),scale(tuple(map(sum,zip(*(g[i] for i in S)))),8)))
    cuts.extend(sub(g[j],g[i]) for i in S for j in rest)
    cuts.append(sub(scale(g[largest],4),ONE))
    cuts.extend(sub(g[largest],g[j]) for j in rest-{largest})
    P=restrict(P,cuts)
    if P is None:
        return Q(0),0
    weight=Q(0)
    copies=0
    for poly,q in paths(P,ID,3):
        weight+=measure(poly,q)
        copies+=1
    return weight,copies


def one_case(indices):
    S,j=tuple(indices[:3]),indices[-1]
    weight,paths_count=check(S,j)
    line=re.search(r'^'+re.escape(str(S))+f' {j}: '+r'(\S+) \((\d+) paths,',
                   (ACCESS/'e1c_row181/verification.log').read_text(), re.M)
    expected = (Q(line[1]),int(line[2])) if line else (0,0)
    assert (weight,paths_count)==expected, (S,j,weight,paths_count,expected)
    h6.cache_clear()
    prepared.cache_clear()
    return {'three_smallest':S,'largest':j,'paths':paths_count,'mass':str(weight)}


def all_cases():
    trials=[(*S,j) for S in combinations(range(6),3) for j in range(6) if j not in S]
    with ProcessPoolExecutor(max_workers=6) as pool:
        cells=list(pool.map(one_case,trials))
    assert len(cells)==60
    weight=sum(Q(row['mass']) for row in cells)
    count=sum(row['paths'] for row in cells)
    assert count==31648 and sum(Q(row['mass'])>0 for row in cells)==23
    print(json.dumps({'total':str(weight),'paths':count,'cells':cells},indent=2))


def locked_snapshot():
    files=[Path(__file__),ACCESS/'e1c_row181/verification.log',
           *(ACCESS/'restricted_theorem_check'/n for n in ('cover_exact.py','exact_model.py'))]
    with tempfile.TemporaryDirectory(prefix='astra148-event-') as tmp:
        for source in files:
            target=Path(tmp)/source.relative_to(ROOT)
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,target)
        with (ROOT/'.task-monitor/compute.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            print('Exclusive lock acquired, integrating all 60 strata',file=sys.stderr,flush=True)
            subprocess.run([sys.executable,'-B',str(Path(__file__).relative_to(ROOT)),'--private'],
                           cwd=tmp,check=True)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('--all',action='store_true')
    args.add_argument('--private',action='store_true',help=argparse.SUPPRESS)
    args.add_argument('indices',nargs='*',type=int,default=[2,3,4,5])
    options=args.parse_args()
    if options.private:
        all_cases()
    elif options.all:
        locked_snapshot()
    else:
        assert len(options.indices)==4
        print(one_case(options.indices))
