"""Attempt whole-box coverage in intrinsic integer checkpoint coordinates."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import cache
from itertools import islice
import json
from math import ceil, lcm
import sys

import cohort1 as orbit
import verify as route
from structural import sections


@dataclass(frozen=True)
class Checkpoint:
    label: tuple
    entry: tuple
    pivots: tuple
    source_intercept: tuple


def intrinsic(checkpoint, guard, label):
    entry = tuple(route.centered(row,v) for row,v in zip(checkpoint,orbit.V))
    first = entry[0]
    _,b,c,_ = first.v
    for j,s in enumerate(entry[1:],1):
        _,p,r,_ = s.v
        determinant = b*r-c*p
        if determinant:
            second = s
            break
    else:
        raise ValueError(('checkpoint lacks a second direction',label))
    a,_,_,k = first.v
    q,_,_,l = second.v
    normal = []
    for row in entry:
        u,v,w,d = row.v
        alpha = Q(v*r-w*p,determinant)
        beta = Q(b*w-c*v,determinant)
        remaining_base = Q(u-alpha*a-beta*q,route.SCALE)
        constant = Q(d-alpha*k-beta*l,route.SCALE)
        normal.append(route.Form.create((remaining_base,alpha,beta),constant))
    return j, tuple(normal), Checkpoint(label,entry,(first,second),(guard.p0,guard.r0))


class CheckpointGuards:
    def __init__(self, key, checkpoints):
        self.key = key
        self.minimum_base = route.START
        self.current = checkpoints[0].label
        self.lattice = None
        self.leading = tuple(row.leading() for row in checkpoints[0].pivots)
        assert all(tuple(row.leading() for row in s.pivots) == self.leading for s in checkpoints)
        offsets = tuple(sorted({
            tuple(row.at((0,*s.source_intercept)) for row in s.pivots) for s in checkpoints
        }))
        self.origin = offsets[0]
        self.offset_denominator = lcm(*(q.denominator for pair in offsets for q in pair))
        self.offsets = hull(tuple(tuple(int(q*self.offset_denominator) for q in pair)
                                 for pair in offsets))

    @cache
    def coordinates(self, chart, base):
        _,p,r = route.leading_coordinates(chart)
        return tuple(row.at((base,base*p,base*r)) + bias
                     for row,bias in zip(self.leading,self.origin))

    @cache
    def leading_guard(self, linear):
        h,b,c = (Q(x,route.SCALE) for x in linear)
        leading = h*route.BASE + b*self.leading[0] + c*self.leading[1]
        return route.coefficient_bounds(leading.linear)

    def parts(self, guard):
        lo,hi = self.leading_guard(guard.linear)
        _,b,c,k = guard.v
        d = self.offset_denominator
        intercept = Q(d*k+min(b*u+c*v for u,v in self.offsets),d*route.SCALE)
        return lo,hi,intercept

    def nonnegative(self, guard, location):
        lo,hi,offset = self.parts(guard)
        if lo > 0:
            self.minimum_base = max(self.minimum_base,ceil(-offset/lo))
            return
        if lo == 0 and offset >= 0:
            return
        raise route.GuardFailure({'guard': guard.record(), 'location': location,
                                  'leading_range': list(map(str,(lo,hi))), 'intercept_min': str(offset),
                                  'lattice': self.lattice.record()})


def hull(points):
    def side(seq):
        turns = []
        for x,y in seq:
            while len(turns) >= 2:
                (u,v),(a,b) = turns[-2:]
                if (a-u)*(y-b) > (b-v)*(x-a):
                    break
                turns.pop()
            turns.append((x,y))
        return turns

    return tuple(sorted(set(side(points)+side(reversed(points)))))


def physical_lattice(index, entry):
    lattice = route.Lattice((5,5,5*(index+1)), ((6,0,0),(0,6,0),(0,0,6)))
    for i,row in enumerate(entry,1):
        congruence = (row-5*i)/6
        children = list(lattice.split(congruence))
        viable = [child for child in children if child.residue(congruence) == 0]
        assert len(viable) == 1 and 1 <= len(children) <= 600, (index,i,lattice,congruence)
        lattice = viable[0]
    return lattice


def weight_units(checkpoint):
    # Each residue modulo 36 covers two residues modulo 72, and cohort 1
    # includes its same-base cohort-0 copy.
    return 4 if checkpoint.label[0] == 1 else 2


def groups(limit=None):
    families = defaultdict(list)
    for label, guard, state, word, drift, _ in islice(sections(),limit):
        index, entry, given = intrinsic(state,guard,label)
        given_word = tuple(route.paired_order(w) for w in word)
        lin = tuple(row.linear for row in given.pivots)
        key = index,entry,given_word,tuple(6*d for d in drift),lin
        families[key].append(given)
    if limit is None:
        assert sum(weight_units(s) for seq in families.values() for s in seq) == 4752
    return families


def certify(key, sources):
    index,entry,word,drift,_ = key
    guard = CheckpointGuards(key,sources)
    lattice = physical_lattice(index,entry)
    clock = -(entry[0]+5*route.BASE)/6
    profiles = set()
    leaves = 0
    greatest_depth = 0
    for end,jumps,_,_ in route.walk_macro(entry,clock,word,drift,lattice,guard,()):
        a,b,c,_ = (Q(x,route.SCALE) for x in end.v)
        clock_at_input = a*route.BASE + b*guard.leading[0] + c*guard.leading[1]
        profiles.add(tuple(clock_at_input.record()[:3]))
        greatest_depth = max(greatest_depth,len(jumps))
        leaves += 1
    return {'leaves':leaves,'profiles':sorted(profiles),'depth':greatest_depth,
            'cutoff':guard.minimum_base}


def main():
    count_only = sys.argv[1:] == ['--count']
    families = groups(int(sys.argv[1]) if len(sys.argv) > 1 and not count_only else None)
    print(f'{len(families)} checkpoint families',file=sys.stderr)
    if count_only:
        print(json.dumps({'families':len(families),
                          'by_multiplicity':dict(sorted(Counter(len(rows) for rows in families.values()).items()))}))
        return
    all_profiles = set()
    misses = []
    units = leaves = successful = cutoff = max_depth = 0
    for i,(key,sources) in enumerate(families.items(),1):
        size = sum(map(weight_units,sources))
        try:
            val = certify(key,sources)
        except route.GuardFailure as error:
            misses.append({'family':i,'labels':[s.label for s in sources],
                           'units':size,'reason':error.args[0]})
        except RuntimeError as error:
            misses.append({'family':i,'labels':[s.label for s in sources],
                           'units':size,'reason':str(error)})
        else:
            all_profiles.update(tuple(p) for p in val['profiles'])
            leaves += val['leaves']
            cutoff = max(cutoff,val['cutoff'])
            max_depth = max(max_depth,val['depth'])
            successful += 1
            units += size
        if i % 10 == 0:
            print(f'{i} families seen, {successful} closed, {units}/4752 mass units, {leaves} leaves',
                  file=sys.stderr)
    print(json.dumps({'families':len(families),'closed_families':successful,'covered_units':units,
                      'leaves':leaves,'cutoff':cutoff,'max_depth':max_depth,
                      'leading_clock_forms':sorted(all_profiles),'failed_families':misses},indent=2))


if __name__ == '__main__':
    main()
