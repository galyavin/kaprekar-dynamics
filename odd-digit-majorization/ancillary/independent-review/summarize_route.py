"""Pool the distinguishing-route runs per residue class; chi-square vs the Section 4.4 table;
homogeneity of the two base families (same B mod 5, opposite parity)."""
import json, math
from scipy.stats import chi2 as C2
R=json.load(open('sim_route_120010.json'))+json.load(open('sim_route_290015.json'))
tot=0; df=0
print("beta  n       freq (7/2, 9/2, 50/9, 37/6)             predicted                               z                      chi2(3)  p")
for b in range(5):
    rs=[r for r in R if r['beta']==b]
    cnt=[sum(r['counts'][i] for r in rs) for i in range(4)]; n=sum(cnt); p=rs[0]['pred']
    z=[(c/n-q)/math.sqrt(q*(1-q)/n) for c,q in zip(cnt,p)]
    ch=sum((c-n*q)**2/(n*q) for c,q in zip(cnt,p)); tot+=ch; df+=3
    print(f"{b}  {n}  {[round(c/n,4) for c in cnt]}  {[round(q,4) for q in p]}  {[round(t,2) for t in z]}  {ch:6.2f}  {C2.sf(ch,3):.3f}")
print(f"all classes: chi2={tot:.2f} df={df} p={C2.sf(tot,df):.3f}")
hc=0;hd=0
for b in range(5):
    a,c=[r for r in R if r['beta']==b]
    na,nc=a['n'],c['n']
    for x,y in zip(a['counts'],c['counts']):
        t=x+y; ea=t*na/(na+nc); ec=t*nc/(na+nc); hc+=(x-ea)**2/ea+(y-ec)**2/ec
    hd+=3
print(f"homogeneity B=120010+b vs 290015+b (opposite parity, different mod 3): chi2={hc:.2f} df={hd} p={C2.sf(hc,hd):.3f}")
print("max clock error:",max(r['max_time_err'] for r in R)," max terminal error:",max(r['max_terminal_err'] for r in R),
      " unmatched:",sum(r['unmatched'] for r in R)," not paired:",sum(r['notpaired'] for r in R))
# per-branch homogeneity (16 inverse branches) pooled over classes, chi2 test vs table per class
bc=0;bd=0
for r in R:
    for bi,c in r['perbranch'].items():
        n=sum(c)
        if n<50: continue
        bc+=sum((ci-n*q)**2/(n*q) for ci,q in zip(c,r['pred'])); bd+=3
print(f"per-branch goodness of fit (16 branches x 10 bases): chi2={bc:.1f} df={bd} p={C2.sf(bc,bd):.3f}")
