"""Cross-check C simulator against the Python column-subtraction model."""
import random, subprocess
from km import F
random.seed(7)
lines=[]; exp=[]
for _ in range(300):
    B=random.choice([10,37,101,1000,99991])
    q=sorted([random.randrange(B) for _ in range(4)],reverse=True)
    lines.append(f"{B} {q[0]} {q[1]} {q[2]} {q[3]}")
    # python orbit
    seen={}; t=tuple(q); i=0
    while t not in seen:
        seen[t]=i; t=F(t,B); i+=1
    mu=seen[t]; lam=i-mu
    # successor of first cyclic
    c=tuple(q)
    for _ in range(mu): c=F(c,B)
    exp.append((mu,lam)+F(c,B))
out=subprocess.run(["./kap9","orbits"],input="\n".join(lines)+"\n",capture_output=True,text=True).stdout.split("\n")
bad=0
for l,e,o in zip(lines,exp,out):
    v=tuple(map(int,o.split()[:6]))
    if v!=e: bad+=1; print("MISMATCH",l,e,v)
print("C vs Python orbits: 300 cases, mismatches:",bad)
