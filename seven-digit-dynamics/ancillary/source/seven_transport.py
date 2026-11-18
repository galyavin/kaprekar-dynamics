"""Exact guarded affine return exploration for positive seven-slot states.

Every accelerated word is checked along its whole affine progression. The
routine does not assert global termination, coverage, or a limiting law. Its
saved segments permit an exact first-cycle lookup, including skipped states.
"""
from odd_gap import subtract_pair, valid

IDENTITY = ((1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0))
B_ROW, ONE = (0,0,0,1,0),(0,0,0,0,1)


def plus(a,b):
    return tuple(x+y for x,y in zip(a,b))


def minus(a,b):
    return tuple(x-y for x,y in zip(a,b))


def evaluate(row,base,q):
    return sum(x*y for x,y in zip(row,(*q,base,1)))


def word_step(base, start, matrix):
    """Extend a word by the uniquely tie-broken actual sorting order."""
    u,v,w = matrix
    q = tuple(evaluate(row,base,start) for row in matrix)
    if not valid(base,q) or q[2] <= 0:
        raise ValueError('Positive last difference is required for this word.')
    guards = [minus(w,ONE),minus(u,v),minus(v,w),minus(minus(B_ROW,ONE),u)]
    digits = [u,v,minus(w,ONE),minus(B_ROW,ONE),
              minus(minus(B_ROW,w),ONE),minus(minus(B_ROW,v),ONE),minus(B_ROW,u)]
    order = sorted(range(7),key=lambda i:(-evaluate(digits[i],base,start),i))
    for i,j in zip(order,order[1:]):
        row = minus(digits[i],digits[j])
        if i>j:
            row = minus(row,ONE)
        guards.append(row)
    result = tuple(minus(digits[order[i]],digits[order[-1-i]]) for i in range(3))
    if tuple(evaluate(row,base,start) for row in result) != subtract_pair(base,q):
        raise ArithmeticError('Affine word disagrees with digit subtraction.')
    return result,tuple(guards),tuple(order)


def guarded_return(base, q, max_word=24):
    """Find a repeatable word; otherwise return one exact step.

    A*delta=delta proves the entire arithmetic progression, even if A is not I.
    Affine guard slopes bound the exact maximal number of valid repetitions.
    """
    matrix, guards, offsets, orders = IDENTITY,[],[],[]
    best = None
    for length in range(1,max_word+1):
        offsets.append(matrix)
        matrix, local, order = word_step(base,q,matrix)
        guards.extend(local)
        orders.append(order)
        image = tuple(evaluate(row,base,q) for row in matrix)
        delta = tuple(a-b for a,b in zip(image,q))
        if not any(delta):
            return {'length':length,'repeats':1,'delta':delta,'matrix':matrix,
                    'guards':tuple(guards),'offsets':tuple(offsets),'orders':tuple(orders),'periodic':True}
        if tuple(sum(row[i]*delta[i] for i in range(3)) for row in matrix) != delta:
            continue
        limits = []
        for row in guards:
            value = evaluate(row,base,q)
            if value < 0:
                raise ArithmeticError('Invalid initial affine guard.')
            slope = sum(row[i]*delta[i] for i in range(3))
            if slope < 0:
                limits.append(value//(-slope)+1)
        if not limits:
            raise ArithmeticError('A nonzero bounded-state progression must exit a guard.')
        repeats = min(limits)
        if repeats >= 2:
            return {'length':length,'repeats':repeats,'delta':delta,'matrix':matrix,
                    'guards':tuple(guards),'offsets':tuple(offsets),'orders':tuple(orders),'periodic':False}
        if best is None:
            best = {'length':length,'repeats':1,'delta':delta,'matrix':matrix,
                    'guards':tuple(guards),'offsets':tuple(offsets),'orders':tuple(orders),'periodic':False}
    if best:
        return best
    matrix, guards, order = word_step(base,q,IDENTITY)
    image = tuple(evaluate(row,base,q) for row in matrix)
    return {'length':1,'repeats':1,'delta':tuple(a-b for a,b in zip(image,q)),
            'matrix':matrix,'guards':guards,'offsets':(IDENTITY,),'orders':(order,),'periodic':False}


def trace(base,start,max_segments=10000,max_word=24):
    """Exact pair arrival plus guarded compressed segments, if cap suffices."""
    q,elapsed,seen,segments = start,0,{},[]
    while q not in seen:
        if len(segments)>=max_segments:
            raise RuntimeError('Segment cap reached; no conclusion about termination.')
        seen[q] = elapsed
        word = guarded_return(base,q,max_word)
        segments.append({'start':q,'time':elapsed,**word})
        q = tuple(x+word['repeats']*d for x,d in zip(q,word['delta']))
        elapsed += word['length']*word['repeats']
        if word['periodic']:
            break
    # At a repeated anchor (or exact periodic word), reconstruct its actual cycle.
    cycle,p = [],q
    while p not in cycle:
        cycle.append(p)
        p = subtract_pair(base,p)
        if len(cycle)>10000:
            raise RuntimeError('Cycle cap reached; no first-arrival claim returned.')
    if p != q:
        raise ArithmeticError('Repeated anchor is not on the reconstructed cycle.')
    arrival = elapsed
    for segment in segments:
        for j,matrix in enumerate(segment['offsets']):
            origin = tuple(evaluate(row,base,segment['start']) for row in matrix)
            velocity = tuple(sum(row[i]*segment['delta'][i] for i in range(3)) for row in matrix)
            for target in cycle:
                index = None
                compatible = True
                for x,d,y in zip(origin,velocity,target):
                    if d==0:
                        if x!=y:
                            compatible=False;break
                    else:
                        k,rem=divmod(y-x,d)
                        if rem or (index is not None and index!=k):
                            compatible=False;break
                        index=k
                if index is None:
                    index=0
                if compatible and 0<=index<segment['repeats']:
                    arrival=min(arrival,segment['time']+index*segment['length']+j)
    return {'arrival':arrival,'cycle':cycle,'segments':segments}
