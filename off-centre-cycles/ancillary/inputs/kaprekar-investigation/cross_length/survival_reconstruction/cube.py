"""Exact phase residues and rational cube volumes for the order-8 cohorts."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import factorial
from pathlib import Path
import json
import sys

from common import (add, complement_period, determinant, insertions, lift, lift_state, literal,
                    loss, periodic_tail, physical_period, root_kick, root_step)

REFERENCES = [
    ((6353,5202,1589,4459,8536,12967,9032,11185),
     [(0,8),(1414,1446),(7214,7478),(9454,9670),(10294,11294),(17091,17139),
      (21493,24197),(29639,31351)]),
    ((40763,7002,61702,54124,81085,1790,53786,36842),
     [(0,8),(1602,1618),(9422,9574),(71560,71688),(93079,93391),(128476,130108),
      (142481,144665),(171781,172149)]),
]


def project(state, active):
    return tuple(0 if j==active else v for j,v in enumerate(state))


def shift(state, slot, value):
    return tuple(value if i==slot else v for i,v in enumerate(state))


def projected_cycle(section, slot):
    x = shift(section,slot,0)
    initial = x,slot
    phases = []
    while True:
        phases.append((x,slot))
        x = root_kick(x,1,0) if slot==len(x)-1 else root_step(x)
        slot = (slot+1)%len(x)
        x = project(x,slot)
        if (x,slot)==initial:
            return phases
        assert len(phases)<10000


def phase_marks(section, slot):
    phases = projected_cycle(section,slot)
    canonical = min(x for x,i in phases if i==len(x)-1)
    length = len(phases)
    x0 = canonical
    for j in range(length):
        active = (len(x0)-1+j)%len(x0)
        if j==0:
            offsets = []
            x = x0
        if active==len(x0)-1:
            offsets.append((x[-1],shift(x,active,0)))
        x = root_kick(x,1,0) if active==len(x0)-1 else root_step(x)
    assert project(x,len(x)-1)==x0
    d = -x[-1]
    lower = 1-min(v for v,_ in offsets)
    assert d>0
    good = []
    terminal_marks = {}

    def good_tail(start):
        chain = {}
        x = start
        accumulated = 0
        while x not in chain and x not in terminal_marks:
            chain[x] = accumulated
            accumulated += loss(x[-1])
            x = root_step(x)
        answer = terminal_marks[x] if x in terminal_marks else accumulated == chain[x]
        terminal_marks.update((state,answer) for state in chain)
        return answer

    # Heights in this interval finish one complete projected return.
    # The first nonpositive read of the next period fixes the unprojected state.
    for height in range(lower,lower+d):
        for drift,state in offsets:
            final = height-d+drift
            if final<=0:
                break
        else:
            raise AssertionError('height escaped the residue window')
        if good_tail(shift(state,len(state)-1,final)):
            good.append(height%d)
    return length,d,canonical,sorted(good)


def solve_rows(rows):
    matrix = [[Q(v) for v in row]+[Q(1)] for row in rows]
    order = len(matrix)
    for j in range(order):
        k = next((k for k in range(j,order) if matrix[k][j]),None)
        if k is None:
            return None
        matrix[j],matrix[k] = matrix[k],matrix[j]
        row = [v/matrix[j][j] for v in matrix[j]]
        matrix[j] = row
        for k in range(order):
            if k!=j and matrix[k][j]:
                matrix[k] = [v-matrix[k][j]*w for v,w in zip(matrix[k],row)]
    return [row[-1] for row in matrix]


def dot(a,b):
    return sum(v*w for v,w in zip(a,b))


def volume(rays):
    """Choose the support of each vertex, then use a pulling triangulation."""
    m = len(rays)
    rows = set(zip(*rays))
    rows = sorted(row for row in rows
                  if not any(row!=other and all(a<=b for a,b in zip(row,other)) for other in rows))
    vertices = {(Q(0),)*m}
    for support_size in range(1,len(rows)+1):
        for support in combinations(range(m),support_size):
            for walls in combinations(rows,support_size):
                result = solve_rows([[row[j] for j in support] for row in walls])
                if result is None or any(v<=0 for v in result):
                    continue
                x = [Q(0)]*m
                for j,v in zip(support,result):
                    x[j] = v
                if all(dot(x,row)<=1 for row in rows):
                    vertices.add(tuple(x))
    vertices = sorted(vertices)
    faces = [{j for j,x in enumerate(vertices) if x[i]==0} for i in range(m)]
    faces += [{j for j,x in enumerate(vertices) if dot(x,row)==1} for row in rows]

    def rank(vectors):
        work = [list(row) for row in vectors]
        basis = 0
        for j in range(m):
            pivot = next((k for k in range(basis,len(work)) if work[k][j]),None)
            if pivot is None:
                continue
            work[basis],work[pivot] = work[pivot],work[basis]
            for k in range(basis+1,len(work)):
                if work[k][j]:
                    ratio = work[k][j]/work[basis][j]
                    work[k] = [v-ratio*w for v,w in zip(work[k],work[basis])]
            basis += 1
        return basis

    @lru_cache(maxsize=None)
    def dim(indices):
        first,*rest = indices
        return rank(add(vertices[j],vertices[first],-1) for j in rest)

    def triangulate(indices,dimension):
        first,*_ = sorted(indices)
        if dimension==0:
            yield (first,)
            return
        candidates = {tuple(sorted(face&indices)) for face in faces if first not in face}
        for facet in sorted(candidates):
            if len(facet)<dimension or dim(facet)!=dimension-1:
                continue
            for simplex in triangulate(set(facet),dimension-1):
                yield (first,*simplex)

    simplices = list(triangulate(set(range(len(vertices))),m))
    total = sum(abs(determinant([add(vertices[j],vertices[s[0]],-1) for j in s[1:]]))
                for s in simplices)/factorial(m)
    return total,len(vertices),len(simplices)


def verify(repo):
    details = []
    retirement_words = set()
    for source,intervals in REFERENCES:
        start,cycle = periodic_tail(source,root_step)
        end = start+len(cycle)
        route = insertions(source,intervals,end,start)
        retired = sorted(range(8),key=lambda j:route.first_retire[j])
        assert len(set(retired))==8
        for i,ray in enumerate(route.rays):
            assert all(v>=0 for v in ray)
            assert not any(ray[j] for j in retired[:i]) and ray[retired[i]]>0
        assert min(route.first_retire)>=8
        for rotation in range(8):
            word = tuple((j-rotation)%8 for j in retired)
            assert word not in retirement_words
            retirement_words.add(word)
        x = source
        for _ in range(intervals[-1][0]):
            x = root_step(x)
        slot = (retired[-1]+intervals[-1][0])%8
        period,descent,canonical,good = phase_marks(x,slot)
        assert route.rays[-1] == tuple(descent*(j==retired[-1]) for j in range(8))
        V,vertices,simplices = volume(route.rays)
        assert V>0
        D_data = lift(51,source,route,{3:2,5:1,9:3,17:0})
        labels,complement,pi,D,delta = D_data
        B = D*1000+11
        g = lift_state(B,source,labels,complement,pi)
        literal(g,50)
        terminal = lift_state(B,route.terminal,labels,complement,pi,
                              tuple(50*j for j in range(sum(map(len,complement)))))
        from math import lcm
        p = lcm(route.period,complement_period(complement))
        literal(terminal,p)
        assert physical_period(terminal,p)==p
        print(f'Cohort {period}: descent={descent}; good={len(good)}; V={V}; '
              f'{vertices=} {simplices=} covered={8*descent*V}; {50+p} literal steps at B={B}',flush=True)
        details.append((period,descent,canonical,good,V,retired,D,delta,route))

    lower = 8*sum(len(good)*V for _,_,_,good,V,*_ in details)
    upper = 1-8*sum((descent-len(good))*V for _,descent,_,good,V,*_ in details)
    print(f'Whole-cube bounds [{lower}, {upper}]',flush=True)

    # Compare accepted numerical output only after both reconstructions.
    output = json.loads((repo/'kaprekar-investigation/cross_length/z1_components/weights.json').read_text())
    for (period,descent,canonical,good,V,retired,D,delta,route), old in zip(details,output['components']):
        assert (period,descent,list(canonical),good,V,retired)==(
            old['period'],old['descent'],old['canonical_section'],old['good_residues'],
            Q(old['unrotated_weight_per_residue']),old['retirement_order'])
        assert abs(determinant(route.rays))==int(old['determinant'])
        assert D in (7997661360,50745451536) and delta==Q(2,4563)
    assert lower==Q(output['zero_loss_lim_inf_lower']) and upper==Q(output['zero_loss_lim_sup_upper'])
    print('Exact archived residues, weights, and bounds agree.',flush=True)


if __name__=='__main__':
    verify(Path(sys.argv[1]))
