"""Exact homogeneous sorting language of eleven-digit leading returns."""

import argparse
import ctypes as ct
from ctypes.util import find_library
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations, product
import json
from math import gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "e1_basins"))
from catalogue_checks import (I, N, difference, dot, family, jacobian, multiply,
                              normal_basis, null_basis, row_times, signed_groups,
                              orbit_subspace, integer_row, verify_orbit)

EARLIER_NULL = {147, 176, 177, 179}


def signed_primitive(row):
    divisor = gcd(*row)
    assert divisor
    return tuple(a // divisor for a in row)


def pull(guards, matrix):
    return tuple(signed_primitive(row_times(g, matrix)) for g in guards)


def negate(row):
    return tuple(-a for a in row)


def library(path=None):
    path = path or find_library("yices")
    candidates = ("/opt/homebrew/lib/libyices.dylib", "/usr/local/lib/libyices.so")
    return ct.CDLL(path or next((p for p in candidates if Path(p).is_file()), "libyices.so"))


class ExactCones:
    """Exact Yices QF_LRA queries, all strict homogeneous guards scaled to >= 1."""

    def __init__(self, lib=None):
        self.lib = library(lib)
        integer, pointer = ct.c_int32, ct.c_void_p
        methods = [
            ("init", None, []), ("new_config", pointer, []),
            ("default_config_for_logic", integer, [pointer, ct.c_char_p]),
            ("new_context", pointer, [pointer]), ("free_context", None, [pointer]),
            ("real_type", integer, []), ("new_uninterpreted_term", integer, [integer]),
            ("int32", integer, [integer]), ("parse_rational", integer, [ct.c_char_p]),
            ("mul", integer, [integer, integer]),
            ("sum", integer, [ct.c_uint32, ct.POINTER(integer)]),
            ("arith_geq_atom", integer, [integer, integer]),
            ("assert_formulas", integer, [pointer, ct.c_uint32, ct.POINTER(integer)]),
            ("check_context", integer, [pointer, pointer]),
        ]
        for name, result, args in methods:
            method = getattr(self.lib, f"yices_{name}")
            method.restype, method.argtypes = result, args
            setattr(self, name, method)
        self.init()
        self.config = self.new_config()
        assert self.default_config_for_logic(self.config, b"QF_LRA") == 0
        self.variables = tuple(self.new_uninterpreted_term(self.real_type()) for _ in range(N))
        self.count = 0
        self.one = self.int32(1)
        self.zero = self.int32(0)

    @lru_cache(maxsize=100_000)
    def linear(self, row):
        terms = []
        for a, x in zip(row, self.variables):
            if a:
                coefficient = self.parse_rational(str(a).encode())
                assert coefficient >= 0
                terms.append(self.mul(coefficient, x))
        return self.sum(len(terms), (ct.c_int32 * len(terms))(*terms)) if terms else self.zero

    @lru_cache(maxsize=100_000)
    def guard_atom(self, row):
        return self.arith_geq_atom(self.linear(row), self.one)

    def check(self, inequalities, equations=(), closed=()):
        terms = [self.guard_atom(g) for g in inequalities]
        terms.extend(self.arith_geq_atom(self.linear(g), self.zero) for g in closed)
        for row in equations:
            terms.extend(self.arith_geq_atom(self.linear(r), self.zero)
                         for r in (row, negate(row)))
        if not terms:
            return True
        assert min(terms) >= 0
        self.count += 1
        context = self.new_context(self.config)
        try:
            assert self.assert_formulas(context, len(terms), (ct.c_int32 * len(terms))(*terms)) == 0
            result = self.check_context(context, None)
            assert result in (3, 4), result
            return result == 3
        finally:
            self.free_context(context)

    def feasible(self, guards):
        return self._feasible(tuple(sorted(set(guards))))

    @lru_cache(maxsize=100_000)
    def _feasible(self, guards):
        signs = set(guards)
        if any(negate(g) in signs for g in guards):
            return False
        return self.check(guards)

    def contains(self, wider, narrower):
        """Both cones are nonempty and open, with nonzero guard rows."""
        narrow_rows = set(narrower)
        return all(g in narrow_rows or not self.feasible((*narrower, negate(g)))
                   for g in wider)

    def facets(self, guards):
        result = list(sorted(set(guards)))
        for g in tuple(result):
            other = tuple(r for r in result if r != g)
            if not self.feasible((*other, negate(g))):
                result.remove(g)
        return tuple(result)


def phases_for(row):
    """Return the minimal cyclic phase period and its exact branch tables."""
    all_matrices, tangent, _, points = family(row)
    period = row["period"]
    branches = []
    axis = tuple(I)
    for j, q in enumerate(points[:period]):
        assert 1 > q[0] > q[1] > q[2] > q[3] > q[4] > 0
        by_phase = []
        groups = signed_groups(q)
        sign = tuple((x > Q(1, 2)) - (x < Q(1, 2)) for x in q)
        for orders in product(*(permutations(group) for group in groups)):
            order = sum(orders, ())
            for signs in product(*((sign[i],) if sign[i] else (1, -1) for i in order)):
                rows = {i: (axis[i] if s > 0 else negate(axis[i])) for i, s in zip(order, signs)}
                guards = [difference(rows[a], rows[b]) for group in orders
                          for a, b in zip(group, group[1:])]
                if sign[order[-1]] == 0:
                    guards.append(rows[order[-1]])
                by_phase.append((jacobian(tuple(zip(order, signs))), tuple(guards)))
        assert tuple(a for a, _ in by_phase) == all_matrices[j]
        branches.append(tuple(by_phase))
    assert all(a == all_matrices[i % period] for i, a in enumerate(all_matrices))
    phase_tangents = []
    for choices in branches:
        phase_tangents.append(tuple(tangent))
        assert all(dot(g, t) == 0 for _, guards in choices for g in guards for t in tangent)
        next_tangent = tuple(tuple(dot(a, v) for a in choices[0][0]) for v in tangent)
        assert all(tuple(dot(a, v) for a in A) == w for A, _ in choices
                   for v, w in zip(tangent, next_tangent))
        tangent = next_tangent
    assert tuple(tangent) == phase_tangents[0]
    return tuple(branches), tuple(phase_tangents)


def words(cones, phases, start, length):
    """Enumerate exactly the jointly strict words of the given finite length."""
    stack = [((), I, ())]
    while stack:
        labels, matrix, guards = stack.pop()
        if len(labels) == length:
            yield labels, matrix, guards
            continue
        phase = (start + len(labels)) % len(phases)
        for label, (branch, local) in reversed(tuple(enumerate(phases[phase]))):
            extended = tuple(sorted(set((*guards, *pull(local, matrix)))))
            if cones.feasible(extended):
                stack.append(((*labels, label), multiply(branch, matrix), extended))


def branch_word(phases, start, labels):
    """Reconstruct the exact matrix and source guard for a labelled word."""
    guards = set()
    matrix = I
    for i, j in enumerate(labels):
        branch, local = phases[(start+i) % len(phases)][j]
        guards.update(pull(local, matrix))
        matrix = multiply(branch, matrix)
    return matrix, tuple(sorted(guards))


def open_rows(records):
    """Recheck the accepted subspace obstruction and 71 positive certificates."""
    zero, closed = [], []
    for i, row in enumerate(records):
        if row["paired"]:
            continue
        phases, tangent, _, _ = family(row)
        candidates = orbit_subspace(phases, tangent)
        if candidates:
            assert any(verify_orbit(phases, tangent, integer_row(v), 64, 4096)
                       for v in candidates), i
            closed.append(i)
        else:
            zero.append(i)
    assert len(zero) == 90 and len(closed) == 71 and EARLIER_NULL <= set(zero)
    return sorted(set(zero)-EARLIER_NULL)


def load_rows(catalogue=HERE.parent / "entry_pi/catalogue.json"):
    return json.loads(catalogue.read_text())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=HERE.parent / "entry_pi/catalogue.json")
    parser.add_argument("--rows", nargs="*", type=int)
    parser.add_argument("--returns", type=int, default=1)
    args = parser.parse_args()
    records = load_rows(args.catalogue)
    choices = args.rows if args.rows is not None else open_rows(records)
    cones = ExactCones()
    for i in choices:
        phases, tangents = phases_for(records[i])
        result = tuple(words(cones, phases, 0, args.returns * len(phases)))
        print(i, len(tangents[0]), [len(choices) for choices in phases], len(result), flush=True)
    print("Exact queries:", cones.count)
