"""Rebuild the P4 and P8 instance kernels from their original route descriptions."""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import lcm
from sys import argv

from physical import literal
from transport import (bounded_prefix, branch, check_guards, exit_children, initial,
                       linear, minus, real_rank, sum_rows, support, times, whole_word)


def profile_form(form):
    return (*form[:-1], 0)


def terminal_kind(sheet):
    """Determine the strip return length and displacement by literal subtraction."""
    B = sheet.base
    point = sheet.point
    seen = {}
    while len(seen) < 100:
        u, v, w, z = point
        offset = (v - 2 * u + B, w - 3 * u + 2 * B, z + u - B)
        if offset in seen:
            t, first = seen[offset]
            return len(seen) - t, u - first
        seen[offset] = (len(seen), u)
        point = literal(point, B)
    raise AssertionError(("strip did not close", sheet.point))


def terminal_route(sheet, reference):
    nine = (-2, -4, -6, -3)
    forty_two = (-10, -20, -30, -20)
    for _ in range(49):
        compiled = branch(sheet, 9, nine)
        if compiled:
            check_guards(sheet, compiled[1], reference)
            return sheet, 9, compiled
        compiled = branch(sheet, 42, forty_two)
        if compiled:
            return sheet, 42, compiled
        word, guards, path = whole_word(sheet.point, sheet.base, 1)
        assert path[-1] == literal(sheet.point, sheet.base)
        check_guards(sheet, guards, reference)
        sheet = sheet.push(word)
    raise AssertionError(("terminal return word not found", sheet.point))


def terminal_q(u, B):
    return (u, minus(times(2, u), B), minus(times(3, u), times(2, B)), minus(B, u))


def p4():
    # The P4 interior witness radius is (18,14,6,3)/25.
    reference = (Q(43, 50), Q(16, 25), Q(2, 5), Q(9, 50))
    Bform = (0, 0, 0, 0, 1, 0)
    L = (Q(7, 4), -1, 0, Q(-1, 4), Q(-3, 4), 0)
    a = (Q(-1, 4), Q(1, 2), 0, Q(-1, 4), Q(3, 4), 0)
    fixed = times(Q(4, 5), Bform)
    span = minus(a, fixed)
    stationary = times(Q(1, 6), sum_rows(a, times(4, Bform)))
    targets = {}
    for c in (Q(7, 2), Q(9, 2), Q(50, 9), Q(37, 6)):
        end = terminal_q(stationary if c == Q(7, 2) else fixed, Bform)
        targets[(sum_rows(L, times(c, span)), end)] = c
    routes = ((3, (-3, -2, -1, -1), 0),
              (12, (-6, -12, -6, -4), 30),
              (24, (-12, -24, -12, -18), 60))
    expected = (("11/48", "5/48", "137/288", "55/288"),
                ("29/128", "5/48", "47/96", "23/128"),
                ("29/144", "5/48", "607/1152", "193/1152"),
                ("107/576", "5/48", "69/128", "197/1152"),
                ("233/1152", "5/48", "73/144", "215/1152"))
    for beta in range(5):
        sheets = [initial(2500000 + beta, 5, reference)]
        sizes = []
        for length, drift, bound in routes:
            children = []
            for sheet in sheets:
                start, compiled, _ = bounded_prefix(sheet, length, drift, bound, reference)
                children.extend(child for _, _, child in exit_children(start, compiled, length, reference))
            sheets = children
            sizes.append(len(sheets))
        result = Counter()
        for sheet in sheets:
            start, length, compiled = terminal_route(sheet, reference)
            if length == 9:
                start_u = profile_form(start.macroscopic[0])
                assert start_u == a
                clock = sum_rows(start.clock, times(Q(9, 2), minus(start_u, fixed)))
                result[targets[clock, terminal_q(fixed, Bform)]] += start.weight
                continue
            for _, _, child in exit_children(start, compiled, length, reference):
                child_leading = tuple(map(profile_form, child.macroscopic[:-2]))
                actual_stationary = child_leading[0]
                assert child_leading == terminal_q(actual_stationary, Bform)
                period, drift = terminal_kind(child)
                if drift:
                    assert (period, drift) in ((32, -2), (37, -3))
                    clock = sum_rows(child.clock, times(Q(period, -drift), minus(actual_stationary, fixed)))
                    terminal = terminal_q(fixed, Bform)
                else:
                    assert period in (1, 2, 6)
                    clock, terminal = child.clock, child_leading
                result[targets[clock, terminal]] += child.weight
        weights = tuple(result[c] for c in sorted(result))
        assert weights == tuple(map(Q, expected[beta])) and sum(weights) == 1
        print(f"[G4] base class {beta} children after three exits {sizes} weights",
              *map(str, weights), flush=True)


def p8():
    # P8's accepted class-six control point is given in relative gap coordinates.
    gap = tuple(map(Q, ("12379624856777/184832467793424", "82001661370921/539094697730820",
                       "2173374886826093/12938272745539680", "9688666835593/40180971259440",
                       "389387410061627/2156378790923280", "494841812881277/2587654549107936")))
    assert sum(gap) == 1
    reference = tuple(sum(gap[i:]) for i in range(1, 6))
    # The eight denominator factors divide this integer, so base 6 is a
    # compatible representative for every guard discovered along the control.
    period = 6 * lcm(*range(1, 31)) ** 2
    sheet = initial(period + 6, period, reference)
    records = (
        (5, (-5, -4, -3, -2, -1), 0, 5, 3, 4),
        (10, (-10, -8, -12, -10, -8), 0, 1, 2, 3),
        (20, (-20, -18, -24, -20, -16), 0, 5, 12, 2),
        (10, (-10, -8, 0, -4, -2), 8, 9, 18, 2),
        (20, (-20, -16, -6, -20, -34), 18, 16, 4, 2),
        (10, (-10, -14, -12, -16, -8), 8, 5, 6, 1),
        (30, (-16, -10, 0, -16, -10), 25, 6, 5, 1),
        (20, (-16, -32, -32, -32, -16), 82, None, 4, 0),
    )
    for step, (length, drift, expected_prefix, wanted, index, rank) in enumerate(records, 1):
        entrant, compiled, prefix = bounded_prefix(sheet, length, drift, expected_prefix, reference)
        assert prefix == expected_prefix
        kids = list(exit_children(entrant, compiled, length, reference))
        assert all(i == index and real_rank(child.columns) == rank for _, i, child in kids)
        print(f"[G8] exit {step} prefix {prefix} rank {rank} index {index} residues",
              sorted(r for r, _, _ in kids), flush=True)
        if wanted is None:
            assert sorted(r for r, _, _ in kids) == [1, 5, 9, 13]
            assert sheet.weight == Q(1, 155520)
            outcomes = {tuple(map(profile_form, c.macroscopic[:-2])) for _, _, c in kids}
            (outcome,) = outcomes
            assert outcome == tuple((0, 0, 0, 0, 0, Q(i, 6), 0) for i in (5, 4, 3, 2, 1))
            centered = [(r, tuple(6 * q - (5 - i) * c.base for i, q in enumerate(c.point)))
                        for r, _, c in kids]
            print("[G8] exact centered exits", sorted(centered))
            assert len({c.clock for _, _, c in kids}) == 1
        else:
            (sheet,) = [c for r, _, c in kids if r == wanted]


def counterexample():
    pairs = {(u % 2, (2 * u + 1) % 4) for u in range(4)}
    first = {a for a, _ in pairs}
    second = {b for _, b in pairs}
    assert pairs == {(0, 1), (1, 3)} and len(first) * len(second) == 4
    assert len(support((1,), 2)) == 2 and len(support((2,), 4)) == 2
    print("tied-guard support", sorted(pairs), "joint weight", Q(1, len(pairs)),
          "product marginal weight", Q(1, len(first) * len(second)))


def borrow_controls():
    tested = 0
    for k in range(1, 5):
        for base in range(2, 9):
            for increasing in combinations_with_replacement(range(base), k):
                q = increasing[::-1]
                word, _, path = whole_word(q, base, 1)
                assert linear(word, (*q, base, 1))[:k] == path[-1] == literal(q, base)
                tested += 1
    print("whole-word digit compiler controls across all heights", tested)


if __name__ == "__main__":
    counterexample()
    borrow_controls()
    if not argv[1:] or "p4" in argv[1:]:
        p4()
    if not argv[1:] or "p8" in argv[1:]:
        p8()
