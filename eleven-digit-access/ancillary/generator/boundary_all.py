"""Unbounded integer local escape checks at the non-Two boundary orbits."""
import argparse
from collections import defaultdict
from pathlib import Path
import subprocess
import z3

ORBIT_DATA = (
    ('fixed', 3, ((2,1,1,1,1),), ((1,4),), 2),
    ('two', 5, ((4,2,1,1,1),(4,3,3,3,2)), ((2,4),(1,3)), 4),
    ('three-open', 28, ((24,16,12,12,7),(24,17,9,4,4),(24,20,20,15,8)),
     ((2,3),(3,4),(1,2)), 9),
    ('three-left', 14, ((12,8,6,6,3),(12,9,5,2,2),(12,10,10,7,4)),
     ((2,3),(3,4),(1,2)), 9),
    ('three-right', 14, ((12,8,6,6,5),(12,7,3,2,2),(12,10,10,9,4)),
     ((2,3),(3,4),(1,2)), 9),
)

def sorted_network(values):
    out=[]
    for value in values:
        for i,x in enumerate(out):
            out[i]=z3.If(x<=value,x,value)
            value=z3.If(x<=value,value,x)
        out.append(value)
    return out

def symbolic_step(B, den, phase, q, successor):
    groups=defaultdict(list)
    groups[den].append(B-1)
    for i,(a,d) in enumerate(zip(phase,q)):
        groups[a].append(d-(i==len(q)-1))
        groups[den-a].append(B-d-(i!=0))
    macro=[]
    digits=[]
    for a,entries in sorted(groups.items(),reverse=True):
        macro.extend([a]*len(entries))
        digits.extend(sorted_network(entries)[::-1])
    assert tuple(macro[i]-macro[-1-i] for i in range(len(q)))==successor
    return [digits[i]-digits[-1-i] for i in range(len(q))]

def prove(name, den, orbit, small, start, steps, out):
    s=z3.SolverFor('QF_LIA')
    s.set(timeout=15000)
    B=z3.Int('B')
    allq=[z3.Ints(' '.join(f'q{j}_{i}' for i in range(5))) for j in range(steps+1)]
    labels=tuple((i+start)%len(orbit) for i in range(steps+1))
    source=orbit[start]
    s.add(*[a>=b for i,(a,b) in enumerate(zip(allq[0],allq[0][1:]))
            if source[i]==source[i+1]])
    def spread(j):
        a,b=small[labels[j]]
        return allq[j][a]-allq[j][b]
    for j in range(steps):
        target=symbolic_step(B,den,orbit[labels[j]],allq[j],orbit[labels[j+1]])
        s.add(*[a==b for a,b in zip(allq[j+1],target)])
        s.add(spread(j+1)<=spread(0))
        status=s.check()
        if status==z3.unsat:
            smt='(set-logic QF_LIA)\n'+s.to_smt2()
            trace=subprocess.run(['yices-smt2'],input=smt,text=True,capture_output=True,check=True)
            assert trace.stdout.strip()=='unsat',trace
            (out/f'{name}-phase-{start}.smt2').write_text(smt)
            print(f'family={name} phase={start} gain within {j+1} updates',flush=True)
            return
    raise ValueError((name,start,steps,status))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    a=p.parse_args()
    a.output.mkdir(exist_ok=True)
    for name,den,orbit,small,steps in ORBIT_DATA:
        for start in range(len(orbit)):
            prove(name,den,orbit,small,start,steps,a.output)
