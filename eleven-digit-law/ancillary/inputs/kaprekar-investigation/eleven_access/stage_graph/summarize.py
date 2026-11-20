"""Read the frozen held-out outputs and score the predeclared observables."""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys

GRID = 2**38
BASES = range(12)
SIZES = (12000, 48000)


def get_actual(file):
    rows = [tuple(map(int,line.split())) for line in file.read_text().splitlines()]
    assert all(len(r)==23 for r in rows)
    return rows


def get_graph(file):
    rows = [tuple([*(int(y) for y in f[:8]), *(Q(y) for y in f[8:10]),
                   *(int(y) for y in f[10:])])
            for f in (line.split() for line in file.read_text().splitlines())]
    assert all(len(r)==15 for r in rows)
    return rows


def means(rows, source, actual):
    # Statistics on the complete no-loss subset only.
    if actual:
        rows = [(Q(r[0],source),r[8:13]) for r in rows if
                r[0]>=0 and r[2]>=0 and r[3]<0]
    else:
        rows = [(r[8],r[10:15]) for r in rows if r[0]==0 and r[2]<0]
    n = len(rows)
    stats = [sum((t for t,x in rows),Q(0))/n]
    stats += [sum((t*x[i] for t,x in rows),Q(0))/(source*n) for i in range(5)]
    return n, stats, sorted(t for t,x in rows)


def ks(xs,ys):
    i=j=score=0
    n,m=len(xs),len(ys)
    while i<n or j<m:
        t=xs[i] if j==m or (i<n and xs[i]<=ys[j]) else ys[j]
        while i<n and xs[i]==t: i+=1
        while j<m and ys[j]==t: j+=1
        score=max(score,abs(i*m-j*n))
    return Q(score,n*m)


def quantiles(sorted_vals):
    if not sorted_vals:
        return None
    return [float(sorted_vals[(len(sorted_vals)-1)*k//100]) for k in (10,50,90,99)]


def compare(a,g,base):
    na, ma, ta = means(a,base,True)
    ng, mg, tg = means(g,GRID,False)
    md = abs(ma[0]-mg[0])
    mixed = [abs(x-y) for x,y in zip(ma[1:],mg[1:])]
    cdf = ks(ta,tg)
    terminal = max(abs(6*r[13+i]-(5-i)*base) for r in a
                   if r[0]>=0 and r[2]>=0 and r[3]<0 for i in range(5))
    nm, nk = ((Q(17,1000),Q(75,1000)) if base < 48000 else (Q(8,1000),Q(55,1000)))
    failures = [r for r in a if r[3]>=0]
    residual = [r for r in g if r[0]]
    result = {
        "actual_n":len(a), "graph_n":len(g), "actual_bulk":na, "graph_bulk":ng,
        "actual_no_entry":sum(r[2]<0 for r in a), "actual_censored":sum(r[0]<0 for r in a),
        "actual_failures":len(failures),
        "actual_failure_mechanisms":Counter(r[6] for r in failures),
        "actual_failure_rows":failures,
        "graph_failures":sum(r[2]>=0 for r in g),
        "graph_failure_stage":Counter(r[2]+1 for r in g if r[2]>=0),
        "graph_residual":Counter(r[0] for r in residual),
        "graph_residual_examples":residual[:5],
        "graph_complete_stage_counts":Counter(r[4] for r in g if not r[0]),
        "graph_rank_preserving_counts":Counter(r[7] for r in g if not r[0]),
        "graph_maximum_bounded_steps":max(r[5] for r in g),
        "bulk_mean_actual":float(ma[0]), "bulk_mean_graph":float(mg[0]),
        "bulk_quantiles_actual":quantiles(ta), "bulk_quantiles_graph":quantiles(tg),
        "mean_discrepancy":float(md), "cdf_discrepancy":float(cdf),
        "mixed_moment_discrepancies":list(map(float,mixed)),
        "terminal_centered_coordinate_max":terminal,
        "pass":(ng*1000>=995*len(g) and max(md,*mixed)<=nm and cdf<=nk and terminal<=2520)
    }
    return result, (ma,ta),(mg,tg)


def paired_comparison(a,b):
    diff = [abs(v-w) for v,w in zip(a[0],b[0])]
    dist = ks(a[1],b[1])
    return {"mean_discrepancy":float(diff[0]), "max_joint_moment_discrepancy":float(max(diff[1:])),
            "ks":float(dist), "pass":max(diff)<=Q(6,1000) and dist<=Q(55,1000)}


def box_test(base,row):
    q = row[8:13]
    a,b,c,d,e = q
    x = [Q(k,base) for k in (a,b,c,2*(b-e)-base,base-c-d)]
    u=(940800,738087,579763,484,2491)
    s=(446,365,342,280,449)
    v=(590,-307,-447,127,-296)
    theta=Q(121049,10**6)+sum(Q(vv,ss)*(xx-Q(uu,10**6))
                            for vv,ss,xx,uu in zip(v,s,x,u))
    return min(Q(1),max(Q(0),500*(Q(row[0],base)-theta+Q(1,1000))))


def main():
    directory=Path(sys.argv[1])
    if "--followup" in sys.argv:
        result={}
        for base in (48001,48007):
            a=get_actual(directory/f"actual-uniform-pi-{base}.tsv")
            g=get_graph(directory/f"graph-uniform-pi-{base}.tsv")
            assert len(a)==len(g)==5000
            result[base]=compare(a,g,base)[0]
        print(json.dumps(result,default=str,indent=2))
        return
    allow_partial="--partial" in sys.argv
    by_source={}
    for actual_mode, graph_mode, actual_count, graph_count in (
            ("original","entry",8000,3000),("uniform-pi","uniform-pi",1000,1000)):
        cases, raw_actual, raw_graph = {},{},{}
        for start,r in itertools.product(SIZES, BASES):
            base=start+r
            f=directory/f"actual-{actual_mode}-{base}.tsv"
            h=directory/f"graph-{graph_mode}-{base}.tsv"
            if allow_partial and (not f.exists() or not h.exists()): continue
            a,g=get_actual(f),get_graph(h)
            if allow_partial and (len(a)!=actual_count or len(g)!=graph_count): continue
            assert len(a)==actual_count and len(g)==graph_count
            cases[base],raw_actual[base],raw_graph[base]=compare(a,g,base)
        six_shifts = {48000+r: {"actual":paired_comparison(raw_actual[48000+r],raw_actual[48006+r]),
                                "graph":paired_comparison(raw_graph[48000+r],raw_graph[48006+r])}
                      for r in range(6) if 48006+r in raw_actual and 48000+r in raw_actual}
        pooled={}
        for start in SIZES:
            v=[c for base,c in cases.items() if start<=base<start+12]
            if not v: continue
            count=sum(c["actual_n"] for c in v)
            losses=sum(c["actual_failures"] for c in v)
            mechanisms=Counter()
            for c in v: mechanisms.update(c["actual_failure_mechanisms"])
            pooled[start]={"actual_draws":count,"actual_failures":losses,
                           "actual_failure_fraction":losses/count,"mechanisms":mechanisms,
                           "actual_no_entry":sum(c["actual_no_entry"] for c in v),
                           "graph_failure_count":sum(c["graph_failures"] for c in v),
                           "graph_draws":sum(c["graph_n"] for c in v),
                           "max_mean_discrepancy":max(c["mean_discrepancy"] for c in v),
                           "max_cdf_discrepancy":max(c["cdf_discrepancy"] for c in v),
                           "gates_passed":sum(c["pass"] for c in v)}
        by_source[actual_mode]={"bases":cases,"six_shifts":six_shifts,"pooled":pooled}
    boxes={}
    for base in (1200005,1200011):
        f=directory/f"actual-box-{base}.tsv"
        if allow_partial and not f.exists(): continue
        a=get_actual(f)
        if allow_partial and len(a)!=24: continue
        assert len(a)==24
        test=[box_test(base,r) for r in a]
        modes=Counter(map(str,test))
        early=sum(r[3]>=0 and 1<=r[6]<=4 for r in a)
        boxes[base]={"count":len(a),"frozen_test":modes,"first_loss_interlevel":early,
                     "minimum_failure_count":min(r[5] for r in a),
                     "maximum_failure_count":max(r[5] for r in a),
                     "pass":early==24 and all(t in (0,1) for t in test),"rows":a}
    print(json.dumps({"sources":by_source,"box":boxes},default=str,indent=2))


if __name__=="__main__":
    main()

