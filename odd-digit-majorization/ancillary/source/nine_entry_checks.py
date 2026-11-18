"""Exact leading-core geometry, guard certificates and independent segments.

Uses SymPy and Z3 from .venv-research. Z3 only finds rational witnesses;
acceptance replays their nonnegative coefficient identities with Fractions.
Neither sampled entry nor a timed bulk segment is a terminal-arrival theorem.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations,product
from math import factorial,gcd,lcm,prod
from pathlib import Path
from random import Random
import copy,json

import sympy as sy
import z3

from nine_entry import (CORE_FACETS,RETURN_MATRIX,BRANCH_GUARDS,DELTA,
    BULK_MATRIX,BULK_GUARDS,BULK_OFFSETS,EXIT_FACETS,H,C,inverse_C,
    in_core_radius,in_core,inverse_leading,leading_bulk,bulk_exit,
    limiting_exit,first_image_segment,evaluate,compose)
from nine_phase_checks import raw_step
from odd_gap import leading_pair


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def normalized(row):
    den=lcm(*(Q(x).denominator for x in row))
    values=[int(Q(x)*den) for x in row];g=gcd(*values)
    return tuple(x//g for x in values) if g else (0,)*len(row)


def core_inventory():
    matrix=tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    guards=set()
    for _ in range(3):
        guards.update(normalized(compose(g,matrix)) for g in BRANCH_GUARDS)
        matrix=tuple(compose(r,matrix) for r in RETURN_MATRIX)
    require(matrix==tuple(tuple(int(i==j) for j in range(5)) for i in range(5)),
            'third return is identity')
    require(set(CORE_FACETS)<=guards,'facet inventory')
    return sorted(guards)


def witness(row,domain):
    variables=z3.Reals(' '.join('c'+str(i) for i in range(len(domain))))
    solver=z3.Solver();solver.set(timeout=10000)
    solver.add(*[v>=0 for v in variables])
    solver.add(*[sum(v*int(r[j]) for v,r in zip(variables,domain))==int(row[j])
                 for j in range(len(row))])
    require(solver.check()==z3.sat,'rational witness search failed')
    model=solver.model();cs=[]
    for v in variables:
        x=model.eval(v,model_completion=True)
        cs.append(Q(x.numerator_as_long(),x.denominator_as_long()))
    return {'row':list(row),'coefficients':list(map(str,cs))}


def verify_witnesses(records,expected,domain):
    require([tuple(r['row']) for r in records]==list(expected),'complete witness inventory')
    for record in records:
        cs=list(map(Q,record['coefficients']))
        require(len(cs)==len(domain) and all(c>=0 for c in cs),'coefficient signs')
        require(tuple(sum(c*r[j] for c,r in zip(cs,domain)) for j in range(5))
                ==tuple(record['row']),'exact rational identity')
    return len(records)


def guard_checks():
    original=core_inventory();domain=(*CORE_FACETS,(0,0,0,0,1))
    reduction=[witness(row,domain) for row in original]
    count=verify_witnesses(reduction,original,domain)
    # Pull the three-word homogeneous guards into radius coordinates.
    change=((Q(1,2),0,0,0,Q(1,2)),(Q(1,2),Q(1,2),0,0,0),
            (0,Q(1,2),Q(1,2),0,0),(0,0,Q(1,2),Q(1,2),0),(0,0,0,0,1))
    word_rows=sorted(set(normalized(compose(g[:5],change)) for g in BULK_GUARDS))
    forward=[witness(row,domain) for row in word_rows]
    count+=verify_witnesses(forward,word_rows,domain)
    reverse_domain=(*word_rows,(0,0,0,0,1))
    backward=[witness(row,reverse_domain) for row in CORE_FACETS]
    count+=verify_witnesses(backward,CORE_FACETS,reverse_domain)
    expected=tuple(tuple(int(i==j) if j<5 else DELTA[i] for j in range(6)) for i in range(4))
    require(BULK_MATRIX==expected,'exact three-step displacement')
    slopes=[evaluate(g[:4],DELTA) for g in BULK_GUARDS]
    error=1+max(Q(abs(g[-1]),-s) for g,s in zip(BULK_GUARDS,slopes) if s<0)
    # The transformed facets have six decreasing slopes, exactly as declared.
    inverse=((2,0,0,0,-1),(-2,2,0,0,1),(2,-2,2,0,-1),
             (-2,2,-2,2,1),(0,0,0,0,1))
    shrinking=set()
    for row in CORE_FACETS:
        qrow=normalized(compose(row,inverse));slope=evaluate(qrow[:4],DELTA)
        if slope<0:shrinking.add((qrow,-slope))
    require(shrinking=={(row,speed) for key,row,speed in EXIT_FACETS},'six exits')
    controls=0
    for mutation in ('delete','coefficient','row'):
        bad=copy.deepcopy(reduction)
        if mutation=='delete':bad.pop()
        elif mutation=='coefficient':bad[0]['coefficients'][0]=str(Q(bad[0]['coefficients'][0])+1)
        else:bad[0]['row'][0]+=1
        try:verify_witnesses(bad,original,domain)
        except ArithmeticError:controls+=1
        else:raise ArithmeticError('Corrupted rational proof accepted')
    return {'core_rows':original,'core_reduction':reduction,'word_rows':word_rows,
            'core_implies_word':forward,'word_implies_core':backward,
            'rational_identities':count,'corruption_controls':controls,
            'bulk_guards':BULK_GUARDS,'bulk_repetition_error_bound':str(error)}


def geometry(domain=CORE_FACETS):
    vertices=set()
    for rows in combinations(domain,4):
        matrix=sy.Matrix([r[:4] for r in rows])
        if not matrix.det():continue
        v=tuple(Q(x) for x in matrix.inv()*sy.Matrix([-r[4] for r in rows]))
        if all(evaluate(row,(*v,1))>=0 for row in domain):vertices.add(v)
    vertices=sorted(vertices)
    faces=[frozenset(i for i,v in enumerate(vertices) if evaluate(r,(*v,1))==0)
           for r in domain]
    @lru_cache(None)
    def dimension(face):
        ids=sorted(face);v=vertices[ids[0]]
        return sy.Matrix([[a-b for a,b in zip(vertices[i],v)] for i in ids[1:]]).rank() if len(ids)>1 else 0
    if domain==CORE_FACETS:
        require(all(dimension(face)==3 for face in faces),'all twelve faces are facets')
    faces=[face for face in faces if face and dimension(face)==3]
    @lru_cache(None)
    def triangulate(face):
        ids=sorted(face);d=dimension(face)
        if len(ids)==d+1:return [tuple(ids)]
        v=ids[0]
        subs={face&f for f in faces if face&f and v not in face&f and dimension(face&f)==d-1}
        return [(v,*t) for f in sorted(subs,key=lambda x:tuple(sorted(x))) for t in triangulate(f)]
    simplices=triangulate(frozenset(range(len(vertices))))
    if domain==CORE_FACETS:
        require((len(vertices),len(simplices))==(16,24),'geometry inventory')
    return vertices,simplices


def integrate(vertices,simplices,transform):
    values=[]
    for v in vertices:
        x=transform(v);g=(1-x[0],*[a-b for a,b in zip(x,x[1:])],x[-1])
        require(min(g)>=0,'density support');values.append(g)
    total=Q(0)
    for simplex in simplices:
        v=vertices[simplex[0]]
        determinant=abs(Q(sy.Matrix([[a-b for a,b in zip(vertices[i],v)] for i in simplex[1:]]).det()))
        require(determinant>0,'nondegenerate simplex')
        polynomial={(0,)*5:Q(1)}
        for j in range(5):
            output={}
            for power,c in polynomial.items():
                for k,i in enumerate(simplex):
                    if not values[i][j]:continue
                    p=tuple(e+int(k==l) for l,e in enumerate(power))
                    output[p]=output.get(p,Q(0))+c*values[i][j]
            polynomial=output
        # Degree five in dimension four: the beta denominator 9! cancels rho9's 9!.
        total+=determinant*sum(c*prod(factorial(p) for p in power) for power,c in polynomial.items())
    return total/16  # all coordinate maps here have absolute determinant 1/16


def compositions(total,parts):
    if parts==1:
        yield (total,);return
    for first in range(total+1):
        for rest in compositions(total-first,parts-1):yield (first,*rest)


def simplex_rule():
    """An independently checked 21-node degree-five rule in dimension four."""
    nodes=[]
    for i in range(3):
        denominator=9-2*i
        coefficient=Q((-1)**i*denominator**5,16*factorial(i)*factorial(9-i))
        for beta in compositions(2-i,5):
            nodes.append((tuple(Q(2*b+1,denominator) for b in beta),coefficient))
    count=0
    for degree in range(6):
        for powers in compositions(degree,4):
            value=sum(c*prod(x**p for x,p in zip(lam[1:],powers)) for lam,c in nodes)
            expected=Q(prod(factorial(p) for p in powers),factorial(4+degree))
            require(value==expected,'independent simplex moment');count+=1
    require((len(nodes),count)==(21,126),'cubature inventory')
    return nodes


def density(x):
    return factorial(9)*prod((1-x[0],*[a-b for a,b in zip(x,x[1:])],x[-1]))


def folded_density(r):
    """Density after sorting |2x_i-1|; checked against all sixteen branches."""
    a,b,c,d=r
    polynomial=(a*b*b*c-2*a*b*c+a*c*c-a*c*d*d+a*d*d
                -b*b*c*c+b*b*c-b*b*d*d+2*b*c*d*d-c*d*d)
    return factorial(9)*(a-1)*polynomial/64


def cubature(vertices,simplices,nodes,fn):
    total=Q(0)
    for simplex in simplices:
        v=vertices[simplex[0]]
        determinant=abs(Q(sy.Matrix([[a-b for a,b in zip(vertices[i],v)] for i in simplex[1:]]).det()))
        for lam,c in nodes:
            x=tuple(sum(t*vertices[j][i] for t,j in zip(lam,simplex)) for i in range(4))
            total+=determinant*c*fn(x)
    return total


def exit_partitions(nodes):
    change=((Q(1,2),0,0,0,Q(1,2)),(Q(1,2),Q(1,2),0,0,0),
            (0,Q(1,2),Q(1,2),0,0),(0,0,Q(1,2),Q(1,2),0),(0,0,0,0,1))
    pieces=[]
    for key,row,speed in EXIT_FACETS:
        extra=[]
        for other,other_row,other_speed in EXIT_FACETS:
            if key==other:continue
            comparison=tuple(speed*x-other_speed*y for x,y in zip(other_row,row))
            extra.append(normalized(compose(comparison,change)))
        domain=tuple(sorted(set((*CORE_FACETS,*extra))))
        vertices,simplices=geometry(domain)
        direct=cubature(vertices,simplices,nodes,lambda r:density(C(r))/16)
        preimage=cubature(vertices,simplices,nodes,folded_density)
        require(direct>0 and preimage>0,'each exit has positive input weight')
        center=tuple(sum(v[i] for v in vertices)/len(vertices) for i in range(4))
        require(limiting_exit(C(center))['faces']==(key,),'strict exit representative')
        pieces.append({'face':key,'domain':domain,'vertices':[[str(x) for x in v] for v in vertices],
                       'simplices':simplices,'radius_center':list(map(str,center)),
                       'direct_mass':str(direct),'one_step_mass':str(preimage)})
        print('exit',key,'one-step mass',float(preimage),flush=True)
    return pieces


def segment_checks(error):
    rng=Random(9162026);checks=[];steps=0;faces=Counter();skipped=0
    for B in (100000,100001,100002,100003,100004):
        accepted=0
        for _ in range(500):
            ds=sorted(rng.randrange(B) for _ in range(9))
            q=tuple(ds[-1-i]-ds[i] for i in range(4))
            try:record=first_image_segment(B,q)
            except ValueError:skipped+=1;continue
            y=record['entry']
            if not in_core(y,B,strict=True):continue
            p=q
            for _ in range(record['elapsed']):p=raw_step(B,p)
            require(p==record['target'],'independent exact segment')
            require(any(evaluate(g,(*p,B,1))<0 for g in BULK_GUARDS),'maximal strict exit')
            asymptotic=limiting_exit(tuple(Q(v,B) for v in y))
            require(abs(Q(record['elapsed']-1)-B*asymptotic['elapsed_over_B'])<=3*error,'bulk clock error')
            leading=record['leading_entry']
            require(max(abs(v-B*w) for v,w in zip(y,leading))<=2,'first-prefix error')
            target=limiting_exit(leading)
            require(abs(Q(record['elapsed'])-B*target['elapsed_over_B'])<=15+3*error,'prefix clock error')
            faces.update(target['faces']);steps+=record['elapsed'];accepted+=1
            checks.append({'base':B,'source':q,'elapsed':record['elapsed'],'target':p,'leading_exit_faces':target['faces']})
            if accepted==4:break
        require(accepted==4,'insufficient buffered entry samples')
    return {'independent_segments':checks,'independent_digit_steps':steps,
            'sampled_face_counts':dict(faces),'rejected_or_outside_samples':skipped}


def main():
    root=Path(__file__).resolve().parent
    guards=guard_checks();vertices,simplices=geometry()
    direct=integrate(vertices,simplices,C)
    branches=[]
    for signs in product((-1,1),repeat=4):
        mass=integrate(vertices,simplices,lambda r:inverse_leading(r,signs))
        branches.append({'signs':signs,'mass':str(mass)})
    predecessor=sum(Q(r['mass']) for r in branches)
    require(direct==Q(2601988443328429,299774320312320000),'direct input mass')
    require(predecessor==Q(1973075435497331,24981193359360000),'one-step input mass')
    # Independent finite polynomial identities verify the folded density and
    # a second integration method before using it on the six exit pieces.
    variables=sy.symbols('a b c d');expanded=0
    for signs in product((-1,1),repeat=4):
        x=([(1+r)/2 for r,s in zip(variables,signs) if s==1]
           +[(1-r)/2 for r,s in reversed(list(zip(variables,signs))) if s==-1])
        expanded+=density(x)/16
    require(sy.expand(expanded-folded_density(variables))==0,'folded density polynomial identity')
    nodes=simplex_rule()
    require(cubature(vertices,simplices,nodes,lambda r:density(C(r))/16)==direct,'independent direct integral')
    require(cubature(vertices,simplices,nodes,folded_density)==predecessor,'independent preimage integral')
    pieces=exit_partitions(nodes)
    require(sum(Q(p['direct_mass']) for p in pieces)==direct,'direct exit partition')
    require(sum(Q(p['one_step_mass']) for p in pieces)==predecessor,'preimage exit partition')
    # Independent checks of all inverse branches and the actual return.
    centers=[tuple(sum(vertices[j][i] for j in simplex)/5 for i in range(4)) for simplex in simplices]
    algebra=0
    for r in [*vertices,*centers]:
        require(in_core_radius(H(r)) and H(H(H(r)))==r,'leading core invariance')
        require(leading_pair(C(r))==C(H(r)),'conjugacy')
        for signs in product((-1,1),repeat=4):
            require(leading_pair(inverse_leading(r,signs))==C(r),'sixteen inverse branches');algebra+=1
    numeric=segment_checks(Q(guards['bulk_repetition_error_bound']))
    certificate={**guards,'vertices':[[str(x) for x in v] for v in vertices],
                 'simplices':simplices,'direct_mass':str(direct),'inverse_branch_masses':branches,
                 'one_step_mass':str(predecessor),'exit_partitions':pieces}
    (root/'nine_entry_certificates.json').write_text(json.dumps(certificate,indent=2)+'\n')
    result={'scope':'Leading invariant core, weighted one-step access and exact bulk exit segments. No global entry, terminal clock, first-capture or full nine-digit law is established.',
            'rational_identities':guards['rational_identities'],'corruption_controls':guards['corruption_controls'],
            'direct_mass':str(direct),'direct_mass_decimal':float(direct),
            'one_step_mass':str(predecessor),'one_step_mass_decimal':float(predecessor),
            'vertices':len(vertices),'simplices':len(simplices),'inverse_checks':algebra,
            'bulk_repetition_error_bound':guards['bulk_repetition_error_bound'],
            'cubature_moments_checked':126,
            'exit_masses':[{k:p[k] for k in ('face','direct_mass','one_step_mass')} for p in pieces],
            'one_step_type_A_mass':str(sum(Q(p['one_step_mass']) for p in pieces if p['face'][0]=='A')),
            'one_step_type_B_mass':str(sum(Q(p['one_step_mass']) for p in pieces if p['face'][0]=='B')),**numeric}
    (root/'nine_entry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    report=(f"PASS: {result['rational_identities']} rational guard identities; three corruptions rejected.\n"
            f"Core original-input mass {float(direct):.12g}; one-step preimage mass {float(predecessor):.12g}.\n"
            f"PASS: {algebra} inverse checks; {len(numeric['independent_segments'])} exact timed segments, "
            f"{numeric['independent_digit_steps']} independent subtraction steps.\n")
    (root/'nine_entry_verification.txt').write_text(report);print(report,end='')


if __name__=='__main__':main()
