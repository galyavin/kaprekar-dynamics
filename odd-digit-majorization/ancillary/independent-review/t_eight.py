"""Sharpness of the nine-step horizon: the family B=34t-9, q=(26t-12,20t-9,14t-6,8t-3)
(taken from the authors' control) checked with the independent model."""
from km import F, G, gaps
for t in (4, 10, 100, 1000, 10**5):
    B=34*t-9; q=(26*t-12,20*t-9,14*t-6,8*t-3); g=gaps(q,B)
    m=min(g); M=max(g); P=B-4*m-min(M,2*m)
    path=[q]
    for _ in range(9): path.append(F(path[-1],B))
    G0=G(q,B); first=[i for i in range(1,10) if G(path[i],B)<G0]
    print(f"t={t}: B={B} P={P} all positive={all(p[3]>0 for p in path)} first descent step={first[0] if first else None}")
