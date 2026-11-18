"""Independent finite checks; no finite experiment certifies global descent."""
from itertools import permutations, product
from pathlib import Path
from random import Random
import hashlib
import json
import subprocess
import tempfile
import z3
import sympy as sp
from odd_gap import subtract_pair, gaps
from nine_clock_reduction import (dispersion, terminal_residual, return_cover,
                                  terminal_tube_count_bound, require)
from nine_zero_return_audit import NETWORK, KEYS, symbolic_cover, accept

ROOT = Path(__file__).resolve().parent


def sorted_network(xs, network):
    xs = list(xs)
    for i,j in network:
        xs[i],xs[j] = max(xs[i],xs[j]),min(xs[i],xs[j])
    return xs


def formula(base,q):
    hi = [q[i] if q[i+1]>0 else q[i]-1 if q[i]>0 else base-1
          for i in range(3)]+[q[3]-1 if q[3]>0 else base-1]
    ds = sorted_network(hi+[base-q[3]-1,base-q[2]-1,base-q[1]-1,base-q[0]],NETWORK)
    return (base-1-ds[7],ds[0]-ds[6],ds[1]-ds[5],ds[2]-ds[4])


def whole_integer(base,digits):
    ds = sorted(digits)
    lo=hi=0
    for a,b in zip(ds,reversed(ds)):
        lo=base*lo+a
        hi=base*hi+b
    out=[]
    value=hi-lo
    for _ in range(9):
        value,r=divmod(value,base)
        out.append(r)
    require(value==0,'Nine output digits')
    out.sort(reverse=True)
    return tuple(out[i]-out[-1-i] for i in range(4))


def counterfamily():
    s=sp.Symbol('s')
    B=4*s+2
    qs=[(3*s+1,2*s+1,s+1,s-1),(3*s+3,2*s+1,2*s,s),
        (3*s+2,2*s+4,s+1,1),(4*s+1,3*s,2*s+1,s+3),
        (4*s,3*s,2*s-2,s-2),(4*s-1,3*s+3,2*s+2,s+2),
        (4*s-2,3*s+1,2*s+2,s)]
    qs=[tuple(map(sp.expand,q)) for q in qs]
    checks=0
    def nonnegative(f):
        nonlocal checks
        f=sp.expand(f)
        require(sp.Poly(f,s).degree()<=1 and f.coeff(s)>=0 and f.subs(s,10)>=0,
                'Affine guard throughout s>=10')
        checks+=1
    values=[]
    for q in qs:
        for f in (B-1-q[0],q[0]-q[1],q[1]-q[2],q[2]-q[3],q[3]-1):
            nonnegative(f)
        gs=sorted(gaps(B,q),key=lambda x:x.subs(s,10),reverse=True)
        for a,b in zip(gs,gs[1:]):
            nonnegative(a-b)
        values.append(sp.expand(sum((4-2*i)*g for i,g in enumerate(gs))))
    for q,nq in zip(qs,qs[1:]):
        u,v,w,x=q
        ds=sorted([u,v,w,x-1,B-1,B-x-1,B-w-1,B-v-1,B-u],
                  key=lambda x:x.subs(s,10),reverse=True)
        for a,b in zip(ds,ds[1:]):
            nonnegative(a-b)
        require(all(sp.expand(ds[i]-ds[8-i]-nq[i])==0 for i in range(4)),
                'Exact symbolic six-step path')
    for g in values:
        nonnegative(g-values[0])
    require(values==[4*s-2,4*s+6,4*s+12,4*s+14,4*s+4,4*s+6,4*s-2],
            'Dispersion path')
    # The initial sorted gaps are (2,s-1,s,s,s+1), so min(M,2m)=4.
    require(sp.expand(B-4*2-4)==4*s-10,'Unbounded source residual')
    return dict(base=str(B),range='integer s>=10',path=[[str(x) for x in q] for q in qs],
                dispersion=[str(x) for x in values],residual='4*s-10',
                affine_guards=checks)


def main():
    if not __debug__:
        raise RuntimeError('Optimized Python is not accepted.')
    audit=json.loads((ROOT/'nine_zero_return_results.json').read_text())
    statements=json.loads((ROOT/'nine_zero_return_obligations.json').read_text())
    accept(audit['records'])
    require(tuple(r['key'] for r in statements)==KEYS,'Obligation keys')
    for r,s in zip(audit['records'],statements):
        require(hashlib.sha256(s['negated_lemma_smtlib'].encode()).hexdigest()==r['sha256'],
                'Saved statement digest')
    network_count=0
    for xs in permutations(range(8)):
        require(sorted_network(xs,NETWORK)==sorted(xs,reverse=True),'Eight-value sorting')
        network_count+=1
    for xs in product((0,1),repeat=8):
        require(sorted_network(xs,NETWORK)==sorted(xs,reverse=True),'Sorting ties')
        network_count+=1
    numeric=0
    for B in range(2,26):
        for u in range(1,B):
            for v in range(u+1):
                for w in range(v+1):
                    for x in range(w+1):
                        q=(u,v,w,x)
                        target=subtract_pair(B,q)
                        require(formula(B,q)==target,'General borrow encoding')
                        require(0<=dispersion(B,q)<=4*B,'Dispersion range')
                        require(terminal_residual(B,q)>=0,'Residual nonnegativity')
                        if x:
                            source_g=sorted(gaps(B,q),reverse=True)
                            target_g=sorted(gaps(B,target),reverse=True)
                            require(all(sum(target_g[:r])<=sum(source_g[:r])+1
                                        for r in range(1,5)),'Ordered partial sums')
                        numeric+=1
    rng=Random(551144)
    whole=0
    for B in (40,41,128,129,1000,1001,10**6,10**6+1,10**12):
        for _ in range(600):
            ds=sorted([rng.randrange(B) for i in range(9)],reverse=True)
            q=tuple(ds[i]-ds[-1-i] for i in range(4))
            require(q[0]>0,'Nonconstant generated input')
            require(whole_integer(B,ds)==formula(B,q)==subtract_pair(B,q),
                    'Independent whole-integer subtraction')
            if q[3]:
                require(dispersion(B,subtract_pair(B,q))<=dispersion(B,q)+8,
                        'Large-base one-step bound')
            whole+=1
    cover_comparisons=0
    for _ in range(800):
        B=rng.choice((40,41,42,43,100,101,10**6,10**6+1))
        rank=rng.choice((2,3))
        q=tuple(sorted([rng.randrange(1,B) for j in range(rank)],reverse=True))+((0,)*(4-rank))
        value=z3.is_true(z3.simplify(symbolic_cover(B,q,rank)))
        require(value==return_cover(B,q),'Symbolic versus numeric cover')
        cover_comparisons+=1
    # Exact finite count of all gap tubes, including all simplex boundaries.
    tube_cases=0
    for B in range(2,26):
        for C in (0,1,2,5):
            count=0
            for a in range(B+1):
                for b in range(B-a+1):
                    for c in range(B-a-b+1):
                        for d in range(B-a-b-c+1):
                            g=(a,b,c,d,B-a-b-c-d);m=min(g)
                            count+=B-4*m-min(max(g),2*m)<=C
            require(count<=terminal_tube_count_bound(B,C),'Terminal tube count')
            tube_cases+=1
    # Exact counterexamples retained even if a longer descent lemma succeeds.
    q=(797,533,353,174);B=1061;path=[q]
    for _ in range(6):
        path.append(subtract_pair(B,path[-1]))
    require(terminal_residual(B,q)>100 and dispersion(B,path[-1])>=dispersion(B,q),
            'False fixed six-step decrement')
    q=(85,57,29,27);B=114;path2=[q]
    for _ in range(6):
        path2.append(subtract_pair(B,path2[-1]))
    require(terminal_residual(B,q)>100 and
            all(dispersion(B,x)>=dispersion(B,q) for x in path2),
            'False six-step lower-hit claim')
    q=(21,21,11,0);path3=[q]
    for _ in range(5):
        path3.append(subtract_pair(43,path3[-1]))
    require(path3[-1][3]==0,'False odd five-step clearance')
    with tempfile.TemporaryDirectory(prefix='nine-clock-') as tmp:
        binary=str(Path(tmp)/'holdout')
        subprocess.run(['c++','-std=c++17','-O3',str(ROOT/'nine_clock_holdout.cpp'),'-o',binary],check=True)
        holdouts=[]
        for mode,bases in (('zero',range(97,105)),('graph',(128,129))):
            for B in bases:
                r=json.loads(subprocess.run([binary,str(B),mode],check=True,capture_output=True,text=True).stdout)
                holdouts.append(r)
                print(r,flush=True)
    family=counterfamily()
    result=dict(network_checks=network_count,complete_small_states=numeric,
                whole_integer_checks=whole,cover_comparisons=cover_comparisons,
                tube_count_cases=tube_cases,holdouts=holdouts,
                six_step_counterfamily=family,
                negative_controls=[dict(base=1061,path=path),dict(base=114,path=path2),
                                   dict(base=43,path=path3)],
                scope='Independent finite checks of accepted reduction; finite lower-hit diagnostics do not prove its unbounded version.')
    (ROOT/'nine_clock_reduction_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: sorting, borrow encoding, independent whole-integer updates, sparse covers and complete finite graphs.',flush=True)


if __name__=='__main__':
    main()
