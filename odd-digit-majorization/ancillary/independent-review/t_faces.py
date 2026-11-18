"""Core facets in y=(u,v',w,z) coordinates; which slacks decrease under the bulk drift;
compare with the six faces A_i,B_i of Section 4.1; R-action on the slacks; Table 1 slack bookkeeping."""
from fractions import Fraction as Fr
from km import CORE
import itertools
# r = Cinv(y): a=2u-1, b=2v-2u+1, c=2w-2v+2u-1, d=2z-2w+2v-2u+1   -> as affine forms (coef_u,coef_v,coef_w,coef_z,const)
ra=[(2,0,0,0,-1),(-2,2,0,0,1),(2,-2,2,0,-1),(-2,2,-2,2,1)]
def row_in_y(row):
    out=[0]*5
    for i in range(4):
        for j in range(5): out[j]+=row[i]*ra[i][j]
    out[4]+=row[4]
    return tuple(out)
faces={'A0':(1,0,0,1,-1),'B0':(0,1,1,0,-1),'A1':(1,1,-1,0,-1),'B1':(1,1,0,-2,-1),'A2':(2,-1,0,0,-1),'B2':(2,-2,2,-1,-1)}
def norm(f):
    from math import gcd
    g=0
    for t in f: g=gcd(g,abs(t))
    return tuple(t//g for t in f)
inv={norm(v):k for k,v in faces.items()}
delta=(-3,-2,-1,-1)
print("core row (r-coords) -> y-coords (u,v',w,z,1), normalised; drift l.delta; matching face")
dec=[]
for row in CORE:
    y=row_in_y(row); n=norm(y)
    ld=sum(n[i]*delta[i] for i in range(4))
    name=inv.get(n,'-')
    print(f"{str(row):22s} -> {str(n):24s} l.delta={ld:3d}  {name}")
    if ld<0: dec.append(name)
print("decreasing facets:",dec)
# R action on slack linear forms: (f o R)(y)
def compR(f):
    u,v,w,z,c=f
    # R(y)=(u,u-z,v-z,v-w): f(Ry)=f_u*u+f_v*(u-z)+f_w*(v-z)+f_z*(v-w)+c
    return (u+v, w+z, -z, -v-w, c)
for k,f in faces.items():
    g=compR(f); print(f"{k} o R = {inv.get(norm(g),g)}")
# Table 1 slack changes
tab={'A2':(12,(-6,-12,-6,-4)),'B0':(6,(-6,-1,1,-2)),'BB':(12,(-12,-2,2,-7)),'AA42':(42,(-10,-20,-30,-20)),'AA9':(9,(-2,-4,-6,-3)),'AB':(24,(-12,-24,-12,-18))}
print("Table 1 slack changes per return (A0,B0,A1,B1,A2,B2):")
for k,(l,d) in tab.items():
    print(f"  {k:5s} l={l:2d} ", {f: sum(faces[f][i]*d[i] for i in range(4)) for f in faces})
# R fixes the bulk drift and v
def Rlin(x): u,v,w,z=x; return (u,u-z,v-z,v-w)
print("R(delta)=",Rlin(delta)," R(4,3,2,1)=",Rlin((4,3,2,1)))
# 24 d1 + 5 delta_AB = 18 d0
d1=(4,8,7,6); dab=(-12,-24,-12,-18); print("24d1+5dAB =",tuple(24*a+5*b for a,b in zip(d1,dab))," 18d0 =",tuple(18*a for a in (2,4,6,3)))
print("(A1)lin d1 =",sum(faces['A1'][i]*d1[i] for i in range(4)),"(A1)lin dAB =",sum(faces['A1'][i]*dab[i] for i in range(4)))
# Q(4/5) on all six faces
q=(Fr(4,5),Fr(3,5),Fr(2,5),Fr(1,5)); print("slacks at Q(4/5):",{k:sum(f[i]*q[i] for i in range(4))+f[4] for k,f in faces.items()})
# the distinguishing radius
r=(Fr(18,25),Fr(14,25),Fr(6,25),Fr(3,25))
from km import Cmap, core_slacks
y=Cmap(r); print("y=C(r)=",y," core slacks:",core_slacks(r))
print("face slacks at y:",{k:sum(f[i]*y[i] for i in range(4))+f[4] for k,f in faces.items()})
u,v,w,z=y; L=Fr(7,4)*u-v-z/4-Fr(3,4); s=-u/4+v/2-z/4+Fr(3,4); print("L =",L," s =",s, float(s))
