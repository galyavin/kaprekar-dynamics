"""Certify a constant list of all difference predecessors for B=5m, m>=20."""

import argparse
from collections import Counter
from functools import total_ordering
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path

from preimages import MACROS, RESIDUALS, basin, step, whole_step, preimages


MIN_M = 20


@total_ordering
class Affine:
    """The integer function slope*m + offset on the complete claimed range."""

    def __init__(self, slope, offset):
        self.slope = slope
        self.offset = offset

    def at(self, m):
        return self.slope*m + self.offset

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Affine) else Affine(0, other)

    def __add__(self, other):
        other = self.coerce(other)
        return Affine(self.slope+other.slope, self.offset+other.offset)

    __radd__ = __add__

    def __neg__(self):
        return Affine(-self.slope, -self.offset)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def sign(self):
        v = self.at(MIN_M)
        if (self.slope > 0 and v <= 0) or (self.slope < 0 and v >= 0):
            raise ArithmeticError(f"Sign can change on certificate range: {self}")
        return (v > 0) - (v < 0)

    def __eq__(self, other):
        return (self-other).sign() == 0

    def __lt__(self, other):
        return (self-other).sign() < 0

    def __hash__(self):
        # Dict and set tests skip equality checks at mismatched hashes.
        # Hashing at the left endpoint forces equality checks on any
        # two distinct functions which meet there. Every comparison
        # that affects membership also explicitly uses affine equality.
        return hash(self.at(MIN_M))

    def __repr__(self):
        return f"Affine({self.slope}, {self.offset})"


def reverse_certificate():
    base = Affine(5, 0)
    phases = [tuple(Affine(*row) for row in zip(p,e)) for p,e in zip(MACROS,RESIDUALS)]
    for index, q in enumerate(phases):
        assert q not in phases[:index]
        assert step(q, base) == phases[(index+1)%4]

    seen = []
    layers = [phases]
    while layers[-1]:
        # Explicit pairwise equality guarantees the validity of membership
        # comparisons at every m, including equalities beyond MIN_M.
        seen.extend(layers[-1])
        layer = []
        for q in layers[-1]:
            for p in preimages(q, base):
                if p not in seen and p not in layer:
                    layer.append(p)
        layers.append(layer)
    assert [len(layer) for layer in layers] == [4, 27, 31, 38, 57, 36, 25, 17, 0]
    small_gaps = Counter()
    weighted_box_factor = 0
    for q in seen:
        gaps = (base-q[0], *(a-b for a,b in zip(q,q[1:])), q[-1])
        gap = Affine.coerce(min(gaps))
        assert gap.slope == 0 and 0 <= gap.offset <= 3
        small_gaps[gap.offset] += 1
        factor = gap.offset if gaps[0] == gap else gap.offset + 1
        weighted_box_factor += factor
    assert small_gaps == {0: 141, 1: 48, 2: 30, 3: 16}
    # Keep the entire inverse certificate as canonical data.
    rows = []
    for layer in layers[:-1]:
        aligned = [tuple(Affine.coerce(x) for x in q) for q in layer]
        rows.append(sorted(tuple((x.slope, x.offset) for x in q) for q in aligned))
    capsule = json.dumps(rows, separators=(',', ':')).encode()
    assert weighted_box_factor == sum((a+1)*b for a,b in small_gaps.items()) == 391
    return layers, rows, weighted_box_factor, hashlib.sha256(capsule).hexdigest()


def independent_controls(layers):
    for m in (20, 21, 25, 40, 101, 1000, 10003):
        base = 5*m
        previous = set()
        expected, sizes = basin(m)
        seen = set()
        for index, layer in enumerate(layers[:-1]):
            points = {tuple(Affine.coerce(x).at(m) for x in q) for q in layer}
            assert len(points) == len(layer)
            for q in points:
                actual = whole_step(q, base)
                assert actual == step(q, base)
                assert actual in (previous if index else points)
            previous = points
            assert not (seen & points)
            seen.update(points)
        assert expected == seen
        assert sizes == [len(layer) for layer in layers]
    for base in (2,3,4,6,9):
        successors = {}
        for ascending in combinations_with_replacement(range(base),5):
            q = tuple(reversed(ascending))
            target = whole_step(q, base)
            assert step(q, base) == target
            if q[0] > 0:
                successors.setdefault(target, set()).add(q)
        for ascending in combinations_with_replacement(range(base),5):
            q = tuple(reversed(ascending))
            assert preimages(q,base) == successors.get(q,set())


def main():
    if not __debug__:
        raise RuntimeError('Assertions are required for this certificate.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    layers, rows, weight_factor, digest = reverse_certificate()
    independent_controls(layers)
    if args.output:
        args.output.write_text(json.dumps(rows, separators=(',', ':'))+'\n')
    print('Layer sizes:', [len(layer) for layer in layers])
    print('Original-word box factor before 11!:', weight_factor)
    print('Canonical affine basin data SHA-256:', digest)
    print('Independent whole-word and exhaustive small-base inverse controls passed.')


if __name__ == '__main__':
    main()
