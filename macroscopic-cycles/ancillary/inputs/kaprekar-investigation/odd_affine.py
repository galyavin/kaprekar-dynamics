"""Small exact affine-word constructor for guarded odd-slot return identities.

Only positive last differences are encoded. It does not assert that a word
applies outside its guards or give global trajectory/arrival conclusions.
"""
from odd_gap import subtract_pair,valid


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def evaluate(row,base,q):
    return sum(a*b for a,b in zip(row,(*q,base,1)))


def word(base,start,length):
    k=len(start);dim=k+2
    matrix=tuple(tuple(int(i==j) for j in range(dim)) for i in range(k))
    one=(0,)*(dim-1)+(1,);B=(0,)*k+(1,0)
    guards,orders,offsets=[],[],[]
    for _ in range(length):
        q=tuple(evaluate(r,base,start) for r in matrix)
        if not valid(base,q) or q[-1]<=0:
            raise ValueError('Require positive last differences throughout the word.')
        offsets.append(matrix)
        guards.extend([sub(matrix[-1],one),sub(sub(B,one),matrix[0]),
                       *[sub(matrix[i],matrix[i+1]) for i in range(k-1)]])
        ds=[*matrix[:-1],sub(matrix[-1],one),sub(B,one),
            *[sub(sub(B,r),one) for r in reversed(matrix[1:])],sub(B,matrix[0])]
        order=sorted(range(2*k+1),key=lambda i:(-evaluate(ds[i],base,start),i))
        guards.extend(sub(sub(ds[i],ds[j]),one) if i>j else sub(ds[i],ds[j])
                      for i,j in zip(order,order[1:]))
        matrix=tuple(sub(ds[order[i]],ds[order[-1-i]]) for i in range(k))
        if tuple(evaluate(r,base,start) for r in matrix)!=subtract_pair(base,q):
            raise ArithmeticError('Affine word disagrees with actual digit columns.')
        orders.append(tuple(order))
    return {'matrix':matrix,'guards':tuple(guards),'orders':tuple(orders),'offsets':tuple(offsets)}


def progression(base,start,certificate):
    """Maximal repetitions on an affine progression when A*delta=delta."""
    M=certificate['matrix'];k=len(start)
    delta=tuple(evaluate(r,base,start)-start[i] for i,r in enumerate(M))
    if tuple(sum(r[i]*delta[i] for i in range(k)) for r in M)!=delta:
        return 1,delta
    limits=[]
    for g in certificate['guards']:
        value=evaluate(g,base,start)
        if value<0:
            raise ValueError('Starting point violates the word guards.')
        slope=sum(g[i]*delta[i] for i in range(k))
        if slope<0:
            limits.append(value//(-slope)+1)
    return min(limits) if limits else 1,delta
