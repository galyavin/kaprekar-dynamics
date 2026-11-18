"""Exact integer row kernels and affine-lattice return transport.

All reductions use tracked unimodular column operations. This module supplies
arithmetic, not a probability theorem without the accompanying lattice-count
and guarded-route arguments.
"""
from fractions import Fraction as Q
from math import gcd
from functools import lru_cache


def identity(n):return tuple(tuple(int(i==j) for j in range(n)) for i in range(n))
def multiply(a,b):return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)
def matvec(a,x):return tuple(sum(v*w for v,w in zip(row,x)) for row in a)


def egcd(a,b):
    r0,r1=a,b;s0,s1=1,0;t0,t1=0,1
    while r1:
        q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1;t0,t1=t1,t0-q*t1
    if r0<0:r0,s0,t0=-r0,-s0,-t0
    return r0,s0,t0


def columns(a,i,j,u,v,w,z):
    for row in a:
        x,y=row[i],row[j];row[i]=u*x+v*y;row[j]=w*x+z*y


def inverse(a):
    n=len(a);rows=[list(map(Q,row))+list(map(Q,identity(n)[i])) for i,row in enumerate(a)]
    for i in range(n):
        j=next((j for j in range(i,n) if rows[j][i]),None)
        if j is None:raise ArithmeticError('Singular reduction matrix.')
        rows[i],rows[j]=rows[j],rows[i];p=rows[i][i];rows[i]=[x/p for x in rows[i]]
        for j in range(n):
            if j!=i:
                p=rows[j][i];rows[j]=[x-p*y for x,y in zip(rows[j],rows[i])]
    return tuple(tuple(row[n:]) for row in rows)


def unimodular_check(u):
    inv=inverse(u)
    if any(x.denominator!=1 for row in inv for x in row) or multiply(u,inv)!=identity(len(u)):
        raise ArithmeticError('Reduction is not integer unimodular.')
    return tuple(tuple(map(int,row)) for row in inv)


@lru_cache(None)
def row_kernel(row):
    """Return g,U with row*U=(g,0,...), U unimodular, g positive."""
    row=tuple(map(int,row));n=len(row);u=[list(r) for r in identity(n)];a=[list(row)]
    for j in range(1,n):
        x,y=a[0][0],a[0][j]
        if not y:continue
        g,s,t=egcd(x,y);columns(a,0,j,s,t,-y//g,x//g);columns(u,0,j,s,t,-y//g,x//g)
    if a[0][0]<0:
        for r in a:r[0]*=-1
        for r in u:r[0]*=-1
    u=tuple(map(tuple,u));g=a[0][0]
    if g<=0 or multiply((row,),u)!=((g,*(0 for _ in range(n-1))),):raise ArithmeticError('Row gcd identity failed.')
    unimodular_check(u);return g,u


@lru_cache(None)
def column_basis(matrix):
    """Column Hermite reduction, retaining the entire integer image lattice."""
    matrix=tuple(tuple(map(int,row)) for row in matrix);m=len(matrix);n=len(matrix[0])
    a=[list(row) for row in matrix];u=[list(row) for row in identity(n)];k=n
    for i in range(m-1,-1,-1):
        if k==0:break
        k-=1
        for j in range(k-1,-1,-1):
            if a[i][j]:
                x,y=a[i][k],a[i][j];g,s,t=egcd(x,y)
                columns(a,k,j,s,t,-y//g,x//g);columns(u,k,j,s,t,-y//g,x//g)
        if a[i][k]<0:
            for row in a:row[k]*=-1
            for row in u:row[k]*=-1
        pivot=a[i][k]
        if not pivot:k+=1;continue
        for j in range(k+1,n):
            q=a[i][j]//pivot
            for row in a:row[j]-=q*row[k]
            for row in u:row[j]-=q*row[k]
    u=tuple(map(tuple,u));a=tuple(map(tuple,a));unimodular_check(u)
    if multiply(matrix,u)!=a or any(x for row in a for x in row[:k]):raise ArithmeticError('Column image identity failed.')
    return tuple(row[k:] for row in a)


def reduce_vector(vector,basis):
    out=list(vector)
    for j in range(len(basis[0])-1,-1,-1):
        pivot=max(i for i,row in enumerate(basis) if row[j]);v=basis[pivot][j]
        if v<=0:raise ArithmeticError('Nonpositive column pivot.')
        q=out[pivot]//v
        for i,row in enumerate(basis):out[i]-=q*row[j]
    return tuple(out)


def canonical(a,c,basis):
    basis=column_basis(basis)
    return reduce_vector(a,basis),reduce_vector(c,basis),basis


def apply_affine(state,matrix,P,r):
    a,c,basis=state
    new_a=tuple(sum(row[i]*a[i] for i in range(4))+row[4]*P for row in matrix)
    new_c=tuple(sum(row[i]*c[i] for i in range(4))+row[4]*r+row[5] for row in matrix)
    new_basis=multiply(tuple(row[:4] for row in matrix),basis)
    return canonical(new_a,new_c,new_basis)


def constant_form(state,row,P,r):
    a,c,basis=state
    if any(x for x in multiply((tuple(row[:4]),),basis)[0]) or sum(x*y for x,y in zip(row[:4],a))+row[4]*P:
        raise ArithmeticError('Requested phase is not constant on the incoming affine lattice.')
    return sum(x*y for x,y in zip(row[:4],c))+row[4]*r


def return_split(state,row,speed,threshold,delta,P,r):
    """All exact outgoing remainder lattices, each with its residue weight."""
    a,c,basis=state;h=tuple(row[:4]);coeff=multiply((h,),basis)[0]
    g,u=row_kernel((*coeff,-speed));ha=sum(x*y for x,y in zip(h,a))+row[4]*P
    hc=sum(x*y for x,y in zip(h,c))+row[4]*r
    if ha%g:raise ArithmeticError(('Base progression needs refinement',P,r,g,ha))
    combined=tuple((*line,d) for line,d in zip(basis,delta));image=multiply(combined,u)
    generator=tuple(line[0] for line in image);new_basis=tuple(line[1:] for line in image)
    new_basis=column_basis(new_basis);out=[]
    for e in range(threshold-speed,threshold):
        if (e-hc)%g:continue
        new_a=tuple(x-v*(ha//g) for x,v in zip(a,generator))
        new_c=tuple(x+v*((e-hc)//g) for x,v in zip(c,generator))
        new=canonical(new_a,new_c,new_basis)
        if constant_form(new,row,P,r)!=e:raise ArithmeticError('Outgoing affine lattice has the wrong remainder.')
        out.append((e,Q(g,speed),new))
    if sum(w for e,w,new in out)!=1:raise ArithmeticError('Outgoing residue weights do not sum to one.')
    return out
