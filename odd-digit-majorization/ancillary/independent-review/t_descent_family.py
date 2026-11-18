"""Targeted test of eq. (descent) around the slow family B=34t-9, q=(26t-12,20t-9,14t-6,8t-3),
where nine steps are needed: random t in [51,1e6], offsets in [-40,40]^5."""
import random
from km import F, G, gaps
random.seed(99)
chk=viol=0; hist=[0]*10
for _ in range(100000):
    t=random.randint(51,10**6)
    B=34*t-9+random.randint(-40,40)
    q=sorted([26*t-12+random.randint(-40,40),20*t-9+random.randint(-40,40),14*t-6+random.randint(-40,40),8*t-3+random.randint(-40,40)],reverse=True)
    q=tuple(q)
    if not (B-1>=q[0] and q[3]>=1): continue
    g=gaps(q,B); m=min(g); M=max(g)
    if B-4*m-min(M,2*m)<=100: continue
    G0=G(q,B); p=q; first=0; pos=True
    for i in range(1,10):
        p=F(p,B)
        if p[3]<=0: pos=False; break
        if not first and G(p,B)<G0: first=i
    if not pos: continue
    chk+=1
    if not first: viol+=1; print("VIOLATION",B,q)
    else: hist[first]+=1
print(f"family neighbourhood: checked={chk} violations={viol} first-descent-step histogram (1..9)={hist[1:]}")
