"""Reconstruct the 'Exit' column of Table 1 from the core geometry alone:
for each canonical face (or face intersection) and its return translation delta,
which decreasing slacks can be the strict first hit (min h/s) from the relative interior."""
import z3
from km import CORE
faces={'A0':(1,0,0,1,-1),'B0':(0,1,1,0,-1),'A1':(1,1,-1,0,-1),'B1':(1,1,0,-2,-1),'A2':(2,-1,0,0,-1),'B2':(2,-2,2,-1,-1)}
ra=[(2,0,0,0,-1),(-2,2,0,0,1),(2,-2,2,0,-1),(-2,2,-2,2,1)]
def row_in_y(row):
    out=[0]*5
    for i in range(4):
        for j in range(5): out[j]+=row[i]*ra[i][j]
    out[4]+=row[4]; return out
rows=[row_in_y(r) for r in CORE]
y=z3.Reals('u v w z')
ev=lambda f: sum(f[i]*y[i] for i in range(4))+f[4]
cases=[('A2',['A2'],(-6,-12,-6,-4)),('B0',['B0'],(-6,-1,1,-2)),('BB',['B0','B1'],(-12,-2,2,-7)),
       ('AA(42)',['A1','A2'],(-10,-20,-30,-20)),('AA(9)',['A1','A2'],(-2,-4,-6,-3)),('AB',['A2','B1'],(-12,-24,-12,-18))]
facename={}
for k,f in faces.items():
    for r,yr in zip(CORE,rows):
        from math import gcd
        g=0
        for t in yr: g=gcd(g,abs(t))
        if tuple(t//g for t in yr)==f: facename[tuple(r)]=k
for name,zero,d in cases:
    s_of={k:-sum(faces[k][i]*d[i] for i in range(4)) for k in faces}
    decr=[k for k in faces if s_of[k]>0]
    reach=[]
    for G in decr:
        s=z3.Solver()
        for r,yr in zip(CORE,rows):
            fn=facename.get(tuple(r))
            if fn in zero: s.add(ev(yr)==0)
            else: s.add(ev(yr)>0)
        for H in decr:
            if H!=G: s.add(s_of[H]*ev(faces[G]) < s_of[G]*ev(faces[H]))
        if s.check()==z3.sat: reach.append(G)
    print(f"{name:7s} zero={zero} decreasing (s): { {k:s_of[k] for k in decr} }  strict first exits: {reach}")
