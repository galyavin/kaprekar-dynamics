"""Rebuild a rank-two exploratory witness as an all-lane radial certificate."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from exact import (HERE, INPUT, add, constant, coordinate, density, evaluate, gaps,
                   in_pi, integrate, multiply, scale, subtract)

ANCHOR = (114876005150,79117929028,50179400739,26486570318,17111906707)
GRID = 137438953472
CENTRE = tuple(Q(x,GRID) for x in ANCHOR)
RADIUS = Q(1,1000)


def paired_verifier(work):
    source = (INPUT/'stage_graph/lottery.cpp').read_text()
    signature = 'void phase_guards(const Phase &p, const Chart &c) {'
    winner = 'void certify(const Context &ctx, const Stage &st, const Chart &c) {'
    assert source.count(signature)==source.count(winner)==1
    source = source.replace(signature,signature+'''
        // At the balanced center the physical cross-level rows are positive.
        // Inside-pair ordering rows have zero explicit base term.
        assert(p.paired);
        for (const auto &a : p.end) assert(a.b==0);
        for (const auto &a : p.guards) assert(a.b>=0);
    ''')
    source = source.replace(winner,winner+'\n        assert(st.exit.b==0);')
    (work/'paired_lottery.cpp').write_text(source)
    shutil.copyfile(INPUT/'stage_graph/graph.cpp',work/'graph.cpp')
    subprocess.run(['clang++','-std=c++17','-O3','-Wall','-Wextra','-Wno-return-type',
                    'paired_lottery.cpp','-o','paired_lottery'],cwd=work,check=True)


def read_panel(work, residue):
    base = 192000+residue
    differences = [str(v+(5-i)*GRID) for i,v in enumerate(ANCHOR)]
    command = ['./paired_lottery',str(base),str(6*GRID),*differences]
    output = subprocess.check_output(command,cwd=work,timeout=60)
    stage = defaultdict(lambda: [Q(0)]*6)
    shapes = defaultdict(set)
    leaves = []
    total = Q(0)
    for line in output.decode().splitlines():
        r = line.split()
        if r[0]=='base':
            radius = Q(r[7]); assert radius>RADIUS
        elif r[0]=='clock':
            w = Q(r[8]); total += w
            form = [Q(v)/int(r[7]) for v in r[1:7]]
            assert form[5]==0
            leaves.append(dict(row=form,weight=w))
        else:
            assert r[0]=='edge', r
            k,n,r0,r1,failures,s,index = map(int,r[1:8])
            assert failures==0
            z = Q(r[20])*Q(n,s*int(r[19]))
            for i,v in enumerate(r[13:19]):
                stage[k][i] += z*int(v)
            shapes[k].add((r0,r1,n,s,index))
    assert total==1
    assert all(sum(c[i] for c in stage.values())==sum(x['row'][i]*x['weight'] for x in leaves)
               for i in range(6))
    return dict(radius=radius,stages=stage,shapes={k:sorted(v) for k,v in shapes.items()},
                leaves=leaves, sha256=sha256(output).hexdigest())


def integration(laws):
    t = coordinate(0)
    x = (t, *(multiply(t,coordinate(i)) for i in range(1,5)))
    q = [scale(add(constant(5-j), x[j]),Q(1,6)) for j in range(5)]
    bounds = [(Q(1,128), CENTRE[0]),
              *[((c-RADIUS)/CENTRE[0],(c+RADIUS)/CENTRE[0]) for c in CENTRE[1:]]]
    verts = list(product(*bounds))
    exp = laws[3]['stages']
    p = laws[5]['stages']
    for k in range(1,5):
        assert exp[k]==p[k]
    assert set(exp)==set(p)==set(range(1,6))
    delta = [a-b for a,b in zip(exp[5],p[5])]
    signed_form = constant(delta[5])
    for f,k in zip(x,delta):
        signed_form = add(signed_form,scale(f,k))
    assert min(evaluate(signed_form,v) for v in verts)>0
    g = gaps(q)
    lows = [g[0],g[1],add(g[0],g[2]),add(g[1],g[3]),add(add(g[0],g[2]),g[4])]
    for v in verts:
        assert in_pi(tuple(evaluate(p,v) for p in q))
        f = [evaluate(p,v) for p in lows]
        assert 0<f[0]<f[1]<f[2]<f[3]<f[4]<Q(1,2)
        for law in laws.values():
            for leaf in law['leaves']:
                assert 0<sum(a*evaluate(p,v) for a,p in zip(leaf['row'],x))<1
    midpoint = tuple((lo+hi)/2 for lo,hi in bounds)
    inverse = []
    for mask in product((0,1),repeat=5):
        pre = [subtract(constant(1),p) for p,b in zip(lows,mask) if b]
        pre += [p for p,b in reversed(list(zip(lows,mask))) if not b]
        if not in_pi(tuple(evaluate(p,midpoint) for p in pre)):
            inverse.append(pre)
    assert len(inverse)==31
    jac = {(4,0,0,0,0):Q(1,6**5)}
    out = {}
    for name, domains in [('direct',[q]),('first_entry_one',inverse)]:
        poly = {}
        for domain in domains:
            poly = add(poly,multiply(jac,density(domain)))
        out[name] = dict(mass=integrate(poly,bounds),
                         contribution=integrate(multiply(poly,signed_form),bounds))
    mass = sum(q['mass'] for q in out.values())
    contribution = sum(q['contribution'] for q in out.values())
    return dict(centre=CENTRE,bounds=bounds, stages_one_through_four_equal=True,
                signed_stage5_form=delta, cohorts=out,mass=mass,contribution=contribution,
                conditional=contribution/mass)


def main():
    with tempfile.TemporaryDirectory(prefix='s1-large-cell-') as directory:
        work = Path(directory)
        paired_verifier(work)
        with ThreadPoolExecutor(max_workers=2) as pool:
            r = (3,5)
            laws = dict(zip(r,pool.map(lambda i:read_panel(work,i),r)))
    print(json.dumps(dict(laws=laws,cone=integration(laws)),indent=2,default=str))


if __name__=='__main__':
    main()
