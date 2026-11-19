from fractions import Fraction as F
import random
ONE=F(1)
def fold(q): return sorted((abs(2*x-1) for x in q), reverse=True)
def C(r): return [ (1+r[0])/2 ]+[ (r[i]+r[i+1])/2 for i in range(4)]
def lead(q): return C(fold(q))
def Hmap(r):
    raw=[r[0]]+[r[i]+r[i+1]-1 for i in range(4)]
    return sorted((abs(x) for x in raw), reverse=True)
def gaps(q): return [1-q[0]]+[q[i]-q[i+1] for i in range(4)]+[q[4]]
def q_from_gaps(g):
    q=[];s=F(1)
    for i in range(5):
        s-=g[i]; q.append(s)
    return q
def energy(q): return sum((x-F(1,6))**2 for x in gaps(q))
def K_sorted(q):
    # independent formula from odd_global_majorization: sort {1,q,1-q}
    s=sorted([ONE]+list(q)+[1-x for x in q], reverse=True)
    return [s[j-1]-s[11-j] for j in range(1,6)]
def orig_q(digits,B):
    s=sorted(digits, reverse=True)
    return [F(s[j-1]-s[11-j],B) for j in range(1,6)]
