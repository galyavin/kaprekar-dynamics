from core import *
EPS=F(1,2**80)
def raw(r): return [r[0]]+[r[i]+r[i+1]-1 for i in range(4)]
def rawoff(d): return [d[0]]+[d[i]+d[i+1] for i in range(4)]
def T_num(ref,d):
    """exact local homogeneous map of H at ref (derived from H itself by epsilon-scaling); d must keep ref's order (ties nonincreasing)"""
    base=Hmap(ref)
    out=[]
    for e in (EPS,EPS/2**20):
        r=[ref[i]+e*d[i] for i in range(5)]
        assert all(r[i]>=r[i+1] for i in range(4)), "input not sorted"
        h=Hmap(r)
        out.append([(h[i]-base[i])/e for i in range(5)])
    assert out[0]==out[1],"epsilon not small enough"
    return out[0]
def T_rule(ref,d):
    """README rule: signed offsets, groups keep order, max/min in tied groups, abs at zero"""
    R=raw(ref); ro=rawoff(d)
    idx=sorted(range(5),key=lambda i:-abs(R[i]))
    out=[];i=0
    zeros=[j for j in idx if R[j]==0]
    nz=[j for j in idx if R[j]!=0]
    k=0
    while k<len(nz):
        m=k
        while m<len(nz) and abs(R[nz[m]])==abs(R[nz[k]]): m+=1
        vals=sorted([ (1 if R[j]>0 else -1)*ro[j] for j in nz[k:m]],reverse=True)
        out+=vals;k=m
    out+=sorted([abs(ro[j]) for j in zeros],reverse=True)
    return out
def rand_dir(ref,rng,scale=1):
    d=[F(rng.randint(-1000,1000),1000)*scale for _ in range(5)]
    # enforce nonincreasing within tied groups of ref
    for i in range(4):
        pass
    i=0
    while i<5:
        j=i
        while j+1<5 and ref[j+1]==ref[i]: j+=1
        if j>i: 
            seg=sorted(d[i:j+1],reverse=True); d[i:j+1]=seg
        i=j+1
    return d
