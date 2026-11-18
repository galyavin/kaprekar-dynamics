"""Brute-force check of the textual gate characterisations in Appendix C (First access)."""
import itertools
from km import F
bad=0; tot=0
for B in list(range(40,66))+[101,102]:
    a=B//2
    for q in itertools.combinations_with_replacement(range(B-1,0,-1),4):
        u,v,w,z=q; tot+=1
        q1=F(q,B)
        zero = q1[3]==0
        if B%2: gate = sum(t==a for t in (v,w,z))>=2
        else: gate = (u==a and v>=a-1 and z<a)
        if zero!=gate: bad+=1; print("gate mismatch",B,q,q1) if bad<5 else None
        if q1[0]!=max(u-1+(u==v), B-z): bad+=1; print("first-coordinate mismatch",B,q,q1) if bad<5 else None
print(f"gate characterisation (odd: two of v',w,z equal a; even: u=a,v'>=a-1,z<a) and first-coordinate formula: {tot} positive sources over B=40..65,101,102; mismatches={bad}")
