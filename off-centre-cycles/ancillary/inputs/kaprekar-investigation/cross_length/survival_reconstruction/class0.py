"""Reconstruct both independent class-0 period-177 source cells."""
from fractions import Fraction as Q
from itertools import product
from math import factorial, prod
import json
import re
from pathlib import Path
import sys

from common import add, determinant, periodic_tail, sorted_word, subtract, word_gaps

ZERO = (0,)*6
B0 = 140000
SEED = (90594,64043,49868,17220,10073,13923)
OUTER = (-1988,1484,1491,-1974,315,322,350)
WALLS = [(13,6,1370,0,15),(8233,12,8,2,14),(8354,30,42,7,15),
         (9616,18,56,11,18),(10724,36,89,12,15)]
RAYS = [(42,35,28,21,14,7),(84,70,7,-7,28,14),(210,28,91,-91,-77,-112),
        (126,105,84,-84,42,21),(147,84,28,-28,-84,63),(84,21,7,-7,-21,63)]
RATES = tuple(map(Q, ('9584/7','169/21','902/21','1192/21','826/9','2063/63')))
WIDTHS = tuple(r/12 for r in RATES)


def gap_offsets(x):
    return add((0,*x),(*x,0),-1)


def slots(x, mask, unit):
    digits = [-unit]
    for i, val in enumerate(x):
        pair = (val-unit*(i==5), -x[5-i]-unit*(i<5))
        digits += pair[::-1] if mask & (1<<i) else pair
    return digits


def mask_at(x, unit):
    raw = slots(x,0,unit)
    return sum((raw[i]<raw[i+1])<<(i//2) for i in range(1,13,2))


def forced(x, mask, unit):
    digits = slots(x,mask,unit)
    return tuple(digits[i]-digits[12-i] for i in range(6))


def ordering(x, mask, unit, B=None):
    digits = slots(x,mask,unit)
    rows = ()
    if B is not None:
        rows = tuple(B+g-unit*(i in (0,6)) for i,g in enumerate(gap_offsets(x)))
        rows += tuple(B+digits[i]-digits[i+1] for i in range(0,12,2))
    return rows+tuple(digits[i]-digits[i+1]-unit*((mask>>(i//2))&1)
                      for i in range(1,13,2))


def follow(x, word, unit):
    for mask in word:
        x = forced(x,mask,unit)
    return x


def analytic_tail(x, unit):
    return periodic_tail(x, lambda z: forced(z,mask_at(z,unit),unit))


def literal_path(x, B, scale, count, start=0):
    g = tuple(Q(B,7)+Q(v,scale) for v in gap_offsets(x))
    assert all(v.denominator == 1 for v in g)
    word = sorted_word(tuple(map(int,g)))
    path, masks = {}, {}
    for i in range(start,start+count+1):
        path[i] = x
        if i == start+count:
            break
        masks[i] = mask_at(x,scale)
        assert min(ordering(x,masks[i],scale,scale*B//7)) >= 0
        word = subtract(word,B)
        offsets = tuple(scale*(v-B//7) for v in word_gaps(word,B))
        x = tuple(-sum(offsets[:j]) for j in range(1,7))
        assert x == forced(path[i],masks[i],scale)
    return path,masks


def tl3(evidence):
    data = json.loads((evidence/'tl3_survival/route.json').read_text())
    x = tuple(data['start'])
    rays = list(map(tuple,data['source_rays']))
    assert abs(determinant(rays)) == 84015792
    C = S = 0

    def bounds(base, directions):
        nonlocal C, S
        for const, *gradient in zip(gap_offsets(base),*(gap_offsets(v) for v in directions)):
            C = max(C,-7*const)
            S = max(S,sum(max(0,-7*d) for d in gradient))

    def affine_word(base, directions, word):
        for mask in word:
            assert min(ordering(base,mask,1)) >= 0
            assert all(min(ordering(ray,mask,0))>=0 for ray in directions)
            bounds(base,directions)
            base = forced(base,mask,1)
            directions = [forced(v,mask,0) for v in directions]
        bounds(base,directions)
        return base,directions

    total = index = 0
    masks = []
    for stage in data['stages']:
        word = stage['masks']
        if 'copies' not in stage:
            x,rays = affine_word(x,rays,word)
            masks += word
            continue
        copies = stage['copies']
        growth = tuple(stage['growth'])
        counts = [row[index] for row in data['count_rays']]
        assert copies > 0 and min(counts) >= 0
        end = add(x,growth,-copies)
        end_rays = [add(v,growth,-c) for v,c in zip(rays,counts)]
        last = add(end,growth)
        assert affine_word(x,rays,word) == (add(x,growth,-1),rays)
        assert affine_word(last,end_rays,word) == (end,end_rays)
        x,rays = end,end_rays
        index += 1
        total += copies*len(word)
        masks += word*copies
    assert total < 3335 and len(masks) == 3335 and x == (0,2,0,1,2,2) and all(v==ZERO for v in rays)
    pre, cycle = analytic_tail(x,1)
    assert pre == 0 and len(cycle) == 177
    for state in cycle:
        bounds(state, [])
    assert max(-7*min(gap_offsets(state)) for state in cycle) == 28
    assert (C,S)==(22834,5663)
    D = 16*S
    cutoff = 16*(C+1)
    assert (D,cutoff)==(90608,365360)
    assert Q(factorial(13),8*9**6*D**6) == Q(3575,1350641342287205909485672071168)
    B = 365365
    path, actual = literal_path(tuple(data['start']), B, 1, len(masks)+177)
    assert list(actual.values())[:len(masks)] == masks and path[len(masks)] == path[len(masks)+177]
    print(f'Task 158: det=84015792 C={C} S={S} D={D} period=177 radius numerator=28 '
          f'{len(actual)} literal updates on class 0', flush=True)


def robust(values):
    """Lower slope on the entire normalized six-direction box."""
    return values[0]-sum(w*abs(v) for w,v in zip(WIDTHS,values[1:]))


def guards_over_box(state, rays, mask, require_depth=False):
    assert min(ordering(state,mask,7,B0)) >= 0
    rows = zip(ordering(rays[0],mask,0,B0), *(ordering(ray,mask,0,0) for ray in rays[1:]))
    assert all(robust(values) >= 0 for values in rows)
    if require_depth:
        outer_witness([(state,rays)])


def outer_witness(endpoints):
    for i in range(7):
        if all(gap_offsets(x)[i] <= -2001 and
               robust(tuple(-gap_offsets(ray)[i] for ray in directions)) >= 0
               for x,directions in endpoints):
            return
    raise AssertionError('no common outside witness')


def push(rays, mask):
    return [forced(ray,mask,0) for ray in rays]


def integral():
    """Expand only terms that survive symmetric integration."""
    center = gap_offsets(SEED)
    derivatives = list(zip(*(gap_offsets(v) for v in RAYS)))
    mass = Q(0)
    for repeated in product(range(-1,6), repeat=7):
        powers = [repeated.count(i) for i in range(6)]
        if any(k%2 for k in powers):
            continue
        numerator = prod(B0+center[j] if i==-1 else derivatives[j][i]
                         for j,i in enumerate(repeated))
        mass += numerator*prod(w**k/(k+1) for w,k in zip(WIDTHS,powers))
    return mass*factorial(13)*prod(2*w for w in WIDTHS)/(7**7*B0**13)


def marked(evidence):
    end = 14108
    path,masks = literal_path(SEED,B0,7,end-13,start=13)
    assert min(min(gap_offsets(x)) for x in path.values()) < -2000
    assert all(min(gap_offsets(path[t])) < -2000 for t in range(13,end))
    assert gap_offsets(path[end]) == OUTER
    assert abs(determinant(RAYS)) > 0
    assert all(63*r == int(63*r) for r in RATES)
    assert all(v%7 == 0 for ray in (SEED,*RAYS) for v in ray)
    assert 8*robust((-gap_offsets(SEED)[0], *[-gap_offsets(v)[0] for v in RAYS])) > B0
    # g0 < B/8 for every base in the lane. The floor remainder is
    # in [0,63 B0), so the tightest scale is n=63.
    g0_slope = robust((-B0-8*gap_offsets(SEED)[0], *[-8*gap_offsets(v)[0] for v in RAYS]))
    assert 64*g0_slope > 63*B0
    for const,*coefficients in zip(gap_offsets(SEED),*(gap_offsets(v) for v in RAYS)):
        assert 18*B0+25*const > 0
        assert robust((18*B0+25*const,*(25*c for c in coefficients))) > 0
    rays = [SEED,*RAYS]
    clock = 13
    all_stages = [*WALLS,(13984,18,6,16,None)]
    for i,(start,length,copies,phase,critical) in enumerate(all_stages):
        for t in range(clock,start):
            guards_over_box(path[t],rays,masks[t],True)
            rays = push(rays,masks[t])
        word = [masks[t] for t in range(start,start+length)]
        drift = add(path[start+length],path[start],-1)
        assert follow(drift,word,0) == drift
        assert all(follow(v,word,0)==v for v in rays)
        next_rays = [add(v,drift,r) for v,r in zip(rays,(RATES[i],*(int(j==i) for j in range(6))))]
        first,last = rays[:],next_rays[:]
        for j,mask in enumerate(word):
            t0,t1 = start+j,start+(copies-1)*length+j
            guards_over_box(path[t0],first,mask)
            guards_over_box(path[t1],last,mask)
            if i<5:
                outer_witness([(path[t0],first),(path[t1],last)])
            first,last = push(first,mask),push(last,mask)
        if critical is not None:
            # This seed row fails at exactly copy "copies", phase "phase".
            # The checked last-copy row has zero slope in the residue lane.
            at_phase = [follow(v,word[:phase],0) for v in next_rays]
            coefficients = [ordering(v,word[phase],0,B0 if j==0 else 0)[critical]
                            for j,v in enumerate(at_phase)]
            assert not any(coefficients)
            advanced = add(path[start+phase],follow(drift,word[:phase],0),copies)
            row = ordering(advanced,word[phase],7,B0)[critical]
            assert row < 0
        else:
            assert all(v==ZERO for v in next_rays)
            first = rays[:]
            for j,mask in enumerate(word):
                d = gap_offsets(follow(drift,word[:j],0))
                last_copy = copies-int(j>=phase)
                first_state,last_state = path[start+j],add(path[start+j],follow(drift,word[:j],0),last_copy)
                witnesses = [k for k,v in enumerate(d) if v>=0 and gap_offsets(last_state)[k] < -2000]
                assert witnesses and any(robust(tuple(-gap_offsets(v)[k] for v in first)) >= 0
                                         for k in witnesses)
                if j<phase:
                    guards_over_box(last_state,next_rays,mask)
                first = push(first,mask)
        if i==0:
            assert set(word)=={0}
            first = rays[:]
            for j,mask in enumerate(word):
                assert min(ordering(first[0],mask,0,B0)) > 0
                assert all(robust(values)>0 for values in zip(
                    ordering(first[0],mask,0,B0),*(ordering(v,mask,0,0) for v in first[1:])))
                first = push(first,mask)
            gate_slack = robust((64*next_rays[0][0]-B0, *(64*v[0] for v in next_rays[1:])))
            assert gate_slack > 0
            assert 64*path[start+copies*length][0]-64*B0+63*gate_slack>0
        rays = next_rays
        clock = start+copies*length
    assert clock+16 == end
    pre,cycle = analytic_tail(path[end],7)
    assert pre==437 and len(cycle)==177
    x = path[end]
    suffix = [x]
    for _ in range(pre+len(cycle)):
        x = forced(x,mask_at(x,7),7)
        suffix.append(x)
    assert max(-min(gap_offsets(x)) for x in suffix)==1988
    assert B0 >= 8*1988+7
    tail,_ = literal_path(path[end],B0,7,pre+len(cycle))
    assert tail[pre] == tail[pre+len(cycle)]
    mu = integral()
    reference = (evidence/'e13b_passage/README.md').read_text()
    numerator, denominator = re.search(r'\\frac\{(377[0-9]+)\}\{([0-9]+)\}', reference).groups()
    assert mu == Q(int(numerator),int(denominator))
    # Repeat the literal route on a fresh eligible progression point, using
    # the scaled center and independent signed lattice variations.
    n = 63
    B = (1+n)*B0+7
    x = tuple((1+n)*v for v in SEED)
    for j,ray in enumerate(RAYS):
        x = add(x,ray,(-1)**j*int(n*WIDTHS[j]))
    count = end-13+sum(length*(int(n*r)+(-1)**j*int(n*WIDTHS[j]))
                       for j,((_,length,*_),r) in enumerate(zip(all_stages,RATES)))
    fresh,_ = literal_path(x,B,7,count)
    assert gap_offsets(fresh[count]) == OUTER and all(min(gap_offsets(fresh[t]))<-2000 for t in range(count))
    print(f'Task 167: det={abs(determinant(RAYS))} first-outer key={OUTER} preperiod={pre} period={len(cycle)} '
          f'literal size={count} {B=} mass={mu}', flush=True)


if __name__=='__main__':
    evidence = Path(sys.argv[1])/'kaprekar-investigation/cross_length'
    tl3(evidence)
    marked(evidence)
