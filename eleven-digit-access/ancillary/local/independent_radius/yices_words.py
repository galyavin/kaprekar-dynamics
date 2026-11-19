import subprocess
from localmap import *
R3=[[F(5,7),F(4,7),F(1,7),F(1,7),F(1,7)],[F(3,4),F(1,2),F(1,4),F(1,4),F(0)],[F(3,4),F(3,4),F(1,4),F(1,4),F(0)]]
def comb(lam): return [sum(l*p[i] for l,p in zip(lam,R3)) for i in range(5)]
def smt(x): return f"(/ {x.numerator} {x.denominator})" if x.denominator!=1 else (str(x) if x>=0 else f"(- {-x})")
def mx(a,b): return f"(ite (>= {a} {b}) {a} {b})"
def mn(a,b): return f"(ite (>= {a} {b}) {b} {a})"
def ab(a): return f"(ite (>= {a} 0) {a} (- {a}))"
def sortdesc(ts):
    if len(ts)<=1: return ts
    if len(ts)==2: return [mx(*ts),mn(*ts)]
    raise NotImplementedError
def step_terms(ref,d):
    """terms of local offset map at ref applied to offset terms d (list of 5 smt strings)"""
    R=raw(ref)
    ro=[d[0]]+[f"(+ {d[i]} {d[i+1]})" for i in range(4)]
    idx=sorted(range(5),key=lambda i:-abs(R[i]))
    zs=[i for i in idx if R[i]==0]; nz=[i for i in idx if R[i]!=0]
    out=[];k=0
    while k<len(nz):
        m=k
        while m<len(nz) and abs(R[nz[m]])==abs(R[nz[k]]): m+=1
        out+=sortdesc([ro[j] if R[j]>0 else f"(- {ro[j]})" for j in nz[k:m]]);k=m
    out+=sortdesc([ab(ro[j]) for j in zs])
    return out
def run(name,asserts,decls,expect):
    L=["(set-logic QF_LRA)"]+[f"(declare-const {v} Real)" for v in decls]+asserts+["(check-sat)"]
    open("query_tmp.smt2","w").write("\n".join(L)+"\n")
    res=subprocess.run(["yices-smt2","query_tmp.smt2"],capture_output=True,text=True).stdout.strip()
    print(f"{name:55s} -> {res}  (expected {expect})")
    return res
V=["a","b","c","d","h"]
def normT(o):
    a,b,c,d,h=o
    return [f"(- {c} {d})", f"(+ {a} {d} {h})", f"(- (- (- (* 2 {a}) {c}) {d}) {h})"]
def FT(n):
    u,v,w=n; L=f"(+ (- {u} {v}) {w})"
    wp=f"(ite (>= {w} 0) {w} 0)"; Lp=f"(ite (>= {L} 0) {L} 0)"
    return [ab(L), f"(- (+ {ab(v)} {wp}) {Lp})", f"(- {L} {ab(w)})"]
def neq_any(A,B): return "(assert (or "+" ".join(f"(not (= {x} {y}))" for x,y in zip(A,B))+"))"
def ordered(ref,vs): return [f"(assert (>= {vs[i]} {vs[i+1]}))" for i in range(4) if ref[i]==ref[i+1]]
# Q1 ray, one step: F_S
f=F(1,10); rs=[F(2,3)+2*f,F(2,3)-2*f,F(1,3),F(1,3),2*f]
o=step_terms(rs,V)
nin=[f"(+ a b)","c","d","(- h a)"]; nout=[f"(+ {o[0]} {o[1]})",o[2],o[3],f"(- {o[4]} {o[0]})"]
FS=[ "(- (- d) (- h a))" ]
v_,c_,d_,z_=nin
FS=[f"(- (- {d_}) {z_})", mx(v_,f"(- (- {c_}) {d_})"), mn(v_,f"(- (- {c_}) {d_})"), f"(- (- {v_}) {c_})"]
run("Q1 ray one-step normal map == F_S",ordered(rs,V)+[neq_any(nout,FS)],V,"unsat")
# Q2-4 triangle
for name,lam in [("Q2 triangle interior",[F(1,3)]*3),("Q3 triangle zero-radius edge (zero entry)",[F(1,2),F(0),F(1,2)]),("Q4 triangle equal-radius edge (tie 2,-3)",[F(1,2),F(1,2),F(0)])]:
    rho=comb(lam); x=V; r=rho
    for k in range(3):
        x=step_terms(r,x); r=Hmap(r)
    run(name+": 3-step normals == F_T",ordered(rho,V)+[neq_any(normT(x),FT(normT(V)))],V,"unsat")
    # sat control: deliberately wrong map (replace |L| by L)
    u,v,w=normT(V); bad=FT(normT(V)); bad[0]=f"(+ (- {u} {v}) {w})"
    run(name+" [control: |L| replaced by L]",ordered(rho,V)+[neq_any(normT(x),bad)],V,"sat")
# Q6 sphere
rA=[F(4,5),F(3,5),F(3,5),F(2,5),F(1,5)]
x=V; r=rA
for k in range(3): x=step_terms(r,x); r=Hmap(r)
ell=lambda t: f"(+ (- {t[2]} {t[1]}) (* 2 {t[3]}) {t[4]})"
asserts=ordered(rA,V)+[f"(assert (<= {V[0]} 0))",f"(assert (>= {ell(V)} 0))",f"(assert (<= {x[0]} 0))",f"(assert (>= {ell(x)} 0))"]
asserts+=[f"(assert (and (<= {v} 1) (>= {v} (- 1))))" for v in V]
asserts+=["(assert (or "+" ".join(f"(= {v} 1) (= {v} (- 1))" for v in V)+"))"]
run("Q6 sphere anchor+energy at both sections, ||x||=1",asserts,V,"unsat")
asserts2=[a for a in asserts if ell(x) not in a]
run("Q6 control: drop energy condition at 2nd section",asserts2,V,"sat")
