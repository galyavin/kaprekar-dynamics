"""Finite-prefix zero-access obligations and independent leading branches.

Only the listed odd-base caps are certified by odd_zero_access_audit.py.
No cap for larger odd slot counts or smaller odd bases is inferred.
"""


def symbolic_positive_prefix(k,depth,minimum_half):
    import z3
    a=z3.Int('a');B=2*a+1;constraints=[a>=minimum_half]
    q=z3.Ints(' '.join(f'q0_{j}' for j in range(k)));states=[q]
    for t in range(depth+1):
        constraints.extend((q[0]<B,q[-1]>=1,*[u>=v for u,v in zip(q,q[1:])]))
        if t==depth:break
        # The central B-1 is maximal; sort the other 2k borrowed digits.
        ds=[*q[:-1],q[-1]-1,*[B-v-1 for v in reversed(q[1:])],B-q[0]]
        for i in range(1,len(ds)):
            for j in range(i,0,-1):
                x,y=z3.Ints(f's{t}_{i}_{j}_hi s{t}_{i}_{j}_lo')
                constraints.extend((x==z3.If(ds[j-1]>=ds[j],ds[j-1],ds[j]),
                                    y==ds[j-1]+ds[j]-x))
                ds[j-1],ds[j]=x,y
        ds=[B-1,*ds];nxt=z3.Ints(' '.join(f'q{t+1}_{j}' for j in range(k)))
        constraints.extend(nxt[j]==ds[j]-ds[-1-j] for j in range(k))
        q=nxt;states.append(q)
    constraints.append(z3.Or(*[z3.And(q[i]==a,q[i+1]==a) for i in range(1,k-1)]))
    return a,B,states,constraints


def folded_branch(base,q):
    """Return M,c with global leading image M*q+c*base.

    A signed permutation gives the folded radii; a unit lower-bidiagonal
    matrix gives their consecutive sums. Thus every branch is unimodular.
    """
    k=len(q);radii=[]
    for i,v in enumerate(q):
        sign=1 if 2*v>=base else -1
        radii.append((max(v,base-v),i,sign))
    radii.sort(key=lambda t:(-t[0],t[1]))
    rows=[];offset=[]
    for _,i,sign in radii:
        rows.append(tuple(sign*int(i==j) for j in range(k)))
        offset.append(int(sign<0))
    M=[rows[0]];c=[offset[0]]
    for j in range(1,k):
        M.append(tuple(a+b for a,b in zip(rows[j-1],rows[j])))
        c.append(offset[j-1]+offset[j]-1)
    return tuple(M),tuple(c)


def positive_first_zero_cap(base,k):
    """An upper bound if a positive source ever reaches zero; 0 means never."""
    if k==1:return 0
    if base%2==0 and base>=8:return 1
    if base%2==1:
        if k==2:return 0
        if k==3 and base>=17:return 2
        if k==4 and base>=33:return 3
    raise ValueError('No certified cap for these parameters.')
