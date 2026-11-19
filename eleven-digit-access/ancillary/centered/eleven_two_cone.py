"""Primitive two-fold return near the whole eleven-slot leading two-cone.

Its clock proof and acceptance records are in eleven_two_cone.md and the
separate audit. Both deep guards are required at each return.
"""
from fractions import Fraction as Q


def offsets(v):
    a,w,p,s,r=v
    return (-a,-2*a+6*p,-3*a+6*s,6*w-4*a,a+6*r)


def coordinates(x):
    a=-x[0]
    return (a,Q(x[3]+4*a,6),Q(x[1]+2*a,6),Q(x[2]+3*a,6),Q(x[4]-a,6))


def guards(v):
    a,w,p,s,r=v
    return w-a+p<=-1 and a-2*w-2*max(0,r+1)>=-1


def step(v):
    a,w,p,s,r=v
    C=max(0,-a+2*s+1)
    A=max(-1,-r)
    J=max(0,-r-1)-max(0,r-1)
    K=max(0,r+1)
    H=max(0,w-J-p-1)
    b=s-p-C
    D=max(A-1,b)
    E=max(A,b-1)-max(-b-1,-A)
    F=min(A,b-1)+J+1
    target=(a-6*D,w+J-p-4*D,E-2*D,F-3*D+H,-p-K+D-H)
    return target,dict(C=C,H=H,D=D,b=b,current=r+J+K+H,
                       gap3_change=-6*(C+max(0,b-A-1)-H))


def reduced_phases():
    """Finite transverse superset forced by three failed native-gap hits."""
    return tuple((p,p+b,r) for p in range(-2,3) for b in range(-4,0)
                 for r in range(-b-4,1))


def reduced(v):
    a,w,p,s,r=v
    b=s-p
    return a>=-27 and w<=2 and -2<=p<=2 and -4<=b<=-1 and -b-4<=r<=0


def fourth_gap(v):
    a,w,p,s,r=v
    return 6*w-5*a-6*r


def lower_hit(v,cap):
    start=fourth_gap(v)
    for h in range(1,cap+1):
        if not guards(v):
            return dict(kind='guard_exit',time=h-1,state=v)
        v,_=step(v)
        if fourth_gap(v)>start:
            return dict(kind='hit',time=h,state=v)
    return dict(kind='no_hit',time=cap,state=v)
