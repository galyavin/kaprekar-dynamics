"""Sanity: column subtraction vs paper's example and vs the borrow formula (eq. borrow)."""
from km import *
import itertools
# Example in Section 2
d = kaprekar_digits([8,6,5,3,0], 10); print("86530 ->", d, qvec(d))
assert d == [8,2,9,6,2] and qvec(d) == (7,6)
print("F(8,3) =", F((8,3),10,5), "gaps src", gaps((8,3),10), "image", gaps((7,6),10))
# borrow formula: all q with q_k>=1, and ranks r<k (copies of B-1), exhaustively
def borrow_digits(q, B, n):
    k=n//2; r=max([i+1 for i in range(k) if q[i]>0], default=0)
    if r==0: return [0]*n
    lead=list(q[:r-1])+[q[r-1]-1]
    tail=[B-1-q[i] for i in range(r-1,0,-1)]+[B-q[0]]
    return lead+[B-1]*(n-2*r)+tail
for n in (3,5,7,9):
    k=n//2
    for B in range(2,14 if n<9 else 11):
        for q in itertools.combinations_with_replacement(range(B-1,-1,-1),k):
            dd=kaprekar_digits(realize(q,n),B)
            assert dd==borrow_digits(q,B,n),(n,B,q,dd)
            # realization independence: random other strings with same q give same image
print("borrow formula == column subtraction: n=3,5,7 (B<14), n=9 (B<11), all q incl. lower ranks")
# independence of realization: all strings at n=5,B=6 and n=9,B=4
for n,B in ((5,6),(7,4),(9,3)):
    for s in itertools.product(range(B),repeat=n):
        q=qvec(s); assert qvec(kaprekar_digits(list(s),B))==F(q,B,n)
print("F well defined on q (all strings checked for (5,6),(7,4),(9,3))")
for n in (3,5,7,9,11,13):
    net=batcher_network(n); assert check_network(net,n); print("batcher",n,len(net),"comparators OK")
