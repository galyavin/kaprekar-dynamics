"""Sensitivity of implication (iv): which hypotheses are load-bearing."""
import time, z3
from zmodel import *
def base_r(s):
    r=[z3.Real(f"r{i}") for i in range(4)]; s.add(1>r[0],r[0]>=r[1],r[1]>=r[2],r[2]>=r[3],r[3]>=0); return r
def its(r,n):
    out=[r]
    for _ in range(n): out.append(H_sym(out[-1]))
    return out
for label,mult,wb in (("min D >= 0 (no gap bound), W<=eta",0,True),("min D >= 100 eta, no W bound",100,False),("min D >= eta/10, W<=eta",z3.Q(1,10),True)):
    s=z3.Solver(); r=base_r(s); eta=z3.Real("eta"); s.add(eta>0)
    if wb: s.add(W_sym(r)<=eta)
    for t in D_sym(r): s.add(t>=mult*eta)
    it=its(r,9); s.add(z3.Not(in_core(it[9])))
    s.add(z3.Or([it[9][i]-it[6][i]!=it[6][i]-it[3][i] for i in range(4)]))
    t=time.time(); res=s.check(); print(f"(iv) mutated: {label:45s} {res}  {time.time()-t:.1f}s",flush=True)
    if res==z3.sat:
        m=s.model(); print("    r =",[m.eval(x) for x in r],"eta =",m[eta])
