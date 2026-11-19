from collections import Counter
from fractions import Fraction
from pathlib import Path
import fcntl
import gzip
import json
import subprocess
import sys
import z3
from faces import layout, permutation, spec_of, templates

K = int(sys.argv[1])
m = K+1
scope = sys.argv[2]
OUT = Path(sys.argv[3])
ID = tuple(range(m))
v = z3.Reals(' '.join(f'g{i}' for i in range(m)))

def translate(g):
    return tuple(sum(g[i:]) for i in range(1,m))

def guards(q,spec):
    middle, groups, endpoint = spec
    radius = [[] for _ in groups]
    middle_guard=[]
    for x,(i,s) in zip(q,layout(spec)):
        if s==0:
            middle_guard.append(2*x==1)
        else:
            radius[i].append(s*(2*x-1))
    reps=[r[0] for r in radius]
    conditions = middle_guard + [x==r[0] for r in radius for x in r[1:]]
    if reps:
        conditions.extend([reps[i]>reps[i+1] for i in range(len(reps)-1)])
        conditions.append(reps[-1]>0)
        conditions.append(reps[0]==1 if endpoint else reps[0]<1)
    return z3.And(conditions)

specs={}
for spec in templates(K):
    p=permutation(K,spec)
    # On a periodic orbit the outer anchor never decreases, hence is fixed.
    if p is not None and p[0]==0:
        specs[spec]=p

def new_solver():
    s=z3.SolverFor('QF_LRA')
    s.add(sum(v)==1,*[a>=0 for a in v])
    if scope=='positive':
        s.add(v[0]>0)
    if scope=='zero':
        s.add(v[0]==0)
    return s

s=new_solver()
guards_by_label={}
stats=Counter()
families={}
commands=[]
statuses=[]
clauses={}

def push():
    s.push()
    commands.append('(push 1)\n')

def pop():
    s.pop()
    commands.append('(pop 1)\n')

def add(*assertions):
    s.add(*assertions)
    for a in assertions:
        key=a.get_id()
        if key not in clauses:
            # Retain the AST to prevent Z3 from reusing its identifier.
            clauses[key]=(a,f'c{len(clauses)}')
        commands.append(f'(assert {clauses[key][1]})\n')

def check():
    res=s.check()
    assert res in (z3.sat,z3.unsat),res
    stats[str(res)]+=1
    statuses.append(str(res))
    commands.append('(check-sat)\n')
    return res

def guard_dict(label):
    if label not in guards_by_label:
        q=translate(tuple(v[i] for i in label))
        guards_by_label[label]={r:guards(q,r) for r in specs}
    return guards_by_label[label]

def val(model,expr):
    x=model.eval(expr)
    return Fraction(x.numerator_as_long(),x.denominator_as_long())

def terminal(label,path,word):
    j=path.index(label)
    for phase in range(j,len(path)):
        if phase and path[phase] in path[j:phase]:
            continue
        push()
        if phase:
            add(*[v[i]==v[l] for i,l in enumerate(path[phase])])
        if check()==z3.sat:
            model = s.model()
            q=tuple(val(model,expr) for expr in translate(v))
            assert spec_of(q)==word[0]
            cycle=word[j:]
            key=tuple(cycle[phase-j:]+cycle[:phase-j])
            if key not in families:
                families[key] = (q,tuple(path[j:]),tuple(word))
        pop()
    stats['terminals']+=1

def visit(label,path,word):
    stats['nodes']+=1
    if label in path:
        terminal(label,path,word)
        return
    gd=guard_dict(label)
    push()
    add(z3.Or(*gd.values()))
    while check()==z3.sat:
        model=s.model()
        q=translate(tuple(v[i] for i in label))
        concrete=tuple(val(model,x) for x in q)
        r=spec_of(concrete)
        if r not in gd:
            raise ArithmeticError((r,concrete))
        push()
        add(gd[r])
        visit(tuple(label[i] for i in specs[r]),path+[label],word+[r])
        pop()
        add(z3.Not(gd[r]))
    pop()

def main():
    visit(ID,[],[])
    record=[]
    for cycle,(q,path,word) in sorted(families.items()):
        record.append({'faces':cycle,'q':list(map(str,q)),
                       'labels':path,'word':word})
    OUT.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Z3',z3.get_version_string(),flush=True)
    print('stats',dict(stats),'family words',len(record),'counts',Counter(map(len,families)),
          'specs',len(specs),flush=True)
    prefix = ['(set-logic QF_LRA)\n']
    prefix.extend(f'(declare-fun {x} () Real)\n' for x in v)
    prefix.extend(f'(define-fun {name} () Bool {expr.sexpr()})\n' for expr,name in clauses.values())
    prefix.extend(f'(assert {a.sexpr()})\n' for a in s.assertions())
    smt=''.join(prefix + commands)
    with gzip.open(OUT.with_suffix('.smt2.gz'),'wt') as fh:
        fh.write(smt)
    yices=subprocess.run(['yices-smt2','--incremental'],input=smt,text=True,capture_output=True,check=True)
    assert yices.stdout.split()==statuses, yices.stderr
    OUT.with_suffix('.statuses').write_text('\n'.join(statuses)+'\n')
    print(f'Yices and Z3 agreed on all {len(statuses)} queries',flush=True)

if __name__=='__main__':
    with open(sys.argv[4]) as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        main()
