"""Generate exact real-linear certificates for the five high-energy cells."""

import argparse
from fractions import Fraction
import gzip
import json
from pathlib import Path


PREAMBLE = """(set-logic QF_LRA)
(define-fun abs_r ((x Real)) Real (ite (>= x 0) x (- x)))
(define-fun max ((x Real) (y Real)) Real (ite (>= x y) x y))
(define-fun min ((x Real) (y Real)) Real (ite (<= x y) x y))
(define-fun pos ((x Real)) Real (max x 0))
"""

SINGULAR_WORD = (((1,), (-5,), (2, -4), (-3,)),)
TRIANGLE_FIRST = ((1, -5), (-4,), (2,), (-3,))
TRIANGLE_SECOND = ((1,), (2, -5), (3,), (-4,))
TRIANGLE_THIRD = SINGULAR_WORD[0]
TRIANGLE_WORDS = {
    "triangle-interior": (TRIANGLE_FIRST, TRIANGLE_SECOND, TRIANGLE_THIRD),
    "triangle-zero-radius-edge": (
        TRIANGLE_FIRST, (*TRIANGLE_SECOND[:-1], (4, -4)), TRIANGLE_THIRD
    ),
    "triangle-equal-radius-edge": (
        ((1, -5), (-4,), (2, -3)), TRIANGLE_SECOND, TRIANGLE_THIRD
    ),
}
SPHERE_WORD = (
    ((1,), (2, -5), (3,), (4, -4)),
    ((1, -5), (-4,), (2, -3)),
    ((1,), (2, -5), (-4,), (3,)),
)


def declare_state(dim):
    return [f"(declare-fun x0_{i} () Real)" for i in range(dim)]


def append_step(lines, step, terms):
    names = [f"x{step}_{i}" for i in range(len(terms))]
    lines.extend(f"(define-fun {name} () Real {term})" for name, term in zip(names, terms))
    return names


def singular(state):
    v, c, d, z = state
    other = f"(- (+ {c} {d}))"
    return f"(- (+ {d} {z}))", f"(max {v} {other})", f"(min {v} {other})", f"(- (+ {v} {c}))"


def triangle(state):
    u, v, w = state
    shift = f"(- (+ {u} {w}) {v})"
    return f"(abs_r {shift})", f"(- (+ (abs_r {v}) (pos {w})) (pos {shift}))", f"(- {shift} (abs_r {w}))"


def singular_normal(state):
    a, b, c, d, e = state
    return f"(+ {a} {b})", c, d, f"(- {e} {a})"


def triangle_normal(state):
    a, b, c, d, e = state
    return f"(- {c} {d})", f"(+ {a} {d} {e})", f"(- (* 2 {a}) (+ {c} {d} {e}))"


def offset_step(state, groups):
    # Different groups have strict separation. Sort only the offsets within
    # a tied group. A group (i, -i) denotes a single zero signed radius.
    raw = [state[0], *(f"(+ {a} {b})" for a, b in zip(state, state[1:]))]
    terms = []
    for group in groups:
        if len(group) == 2 and group[0] == -group[1]:
            terms.append(f"(abs_r {raw[group[0] - 1]})")
            continue
        signed = [raw[i - 1] if i > 0 else f"(- {raw[-i - 1]})" for i in group]
        if len(signed) == 1:
            terms.extend(signed)
        else:
            left, right = signed
            terms.extend((f"(max {left} {right})", f"(min {left} {right})"))
    assert len(terms) == len(state)
    return terms


def full_return(lines, start, word):
    end = start
    for step, groups in enumerate(word, 1):
        end = append_step(lines, step, offset_step(end, groups))
    return end


def check_normal_return(word, normal, update):
    start = [f"x0_{i}" for i in range(5)]
    lines = declare_state(5)
    end = full_return(lines, start, word)
    differences = " ".join(
        f"(distinct {actual} {claimed})"
        for actual, claimed in zip(normal(end), update(normal(start)))
    )
    lines.append(f"(assert (or {differences}))")
    return lines


def normalize(lines, state):
    # Positive homogeneity allows normalization to the boundary of the unit box.
    lines.extend(f"(assert (<= (abs_r {x}) 1))" for x in state)
    unit = " ".join(f"(= (abs_r {x}) 1)" for x in state)
    lines.append(f"(assert (or {unit}))")


def expansion(dim, guard, update, count):
    start = [f"x0_{i}" for i in range(dim)]
    lines = declare_state(dim)
    normalize(lines, start)
    lines.append(f"(assert {guard})")
    end = start
    for step in range(1, count + 1):
        end = append_step(lines, step, update(end))
    lines.extend(f"(assert (<= (abs_r {x}) 2))" for x in end)
    return lines


def sphere():
    dim = 5
    start = [f"x0_{i}" for i in range(dim)]
    lines = declare_state(dim)
    normalize(lines, start)
    lines.append("(assert (>= x0_1 x0_2))")
    end = full_return(lines, start, SPHERE_WORD)
    # With anchor at most its limiting radius, compare the first energy
    # variation at both sections. See the README for the quadratic remainder.
    for state in (start, end):
        a, b, c, d, e = state
        lines.append(f"(assert (<= {a} 0))")
        lines.append(f"(assert (>= (+ (- {b}) {c} (* 2 {d}) {e}) 0))")
    return lines


def radius(q):
    return tuple(sorted((abs(2 * x - 1) for x in q), reverse=True))


def radius_step(r):
    return (r[0], *(a + b - 1 for a, b in zip(r, r[1:])))


def pattern(raw):
    groups = {}
    for i, x in enumerate(raw, 1):
        indices = (i,) if x > 0 else (-i,) if x < 0 else (i, -i)
        groups.setdefault(abs(x), []).extend(indices)
    return tuple(tuple(group) for _, group in sorted(groups.items(), reverse=True))


def check_catalogue(path):
    with gzip.open(path, "rt") as stream:
        rows = json.load(stream)
    high = [row for row in rows if row["location"] == "interior" and row["meets_high_energy"]]
    expected = {
        SINGULAR_WORD: "Singular ray",
        SPHERE_WORD: "isolated sphere orbit",
        **{word: label for label, word in TRIANGLE_WORDS.items()},
    }
    assert len(high) == len(expected)
    for row in high:
        r = radius(map(Fraction, row["q_witness"]))
        start = r
        word = []
        for _ in range(row["witness_period"]):
            raw = radius_step(r)
            word.append(pattern(raw))
            r = tuple(sorted(map(abs, raw), reverse=True))
        assert r == start
        label = expected.pop(tuple(word))
        print(f"{label}: {word}")
    assert not expected


def borrowed_step(q, base):
    k = len(q)
    assert base > q[0] >= q[-1] >= 1
    forward = tuple(value - int(i == k - 1) for i, value in enumerate(q))
    reflected = tuple(base - value - int(i < k - 1) for i, value in enumerate(reversed(q)))
    paired_guard = tuple(
        min(forward[i], reflected[i]) - max(forward[i + 1], reflected[i + 1])
        for i in range(k - 1)
    )
    digits = sorted((*forward, base - 1, *reflected), reverse=True)
    return tuple(digits[i] - digits[-1 - i] for i in range(k)), paired_guard


def check_coincidence_controls():
    q = (420, 179)
    margins = []
    five_guards = []
    for _ in range(5):
        margins.append(sum(q) - 600)
        q, guards = borrowed_step(q, 600)
        five_guards.append(guards[0])
    assert margins == [-1, 61, 2, 58, -1]
    assert five_guards == [240, 120, 236, 122, 236]
    print(f"five-digit signed leading margins {margins}, borrowed guards {five_guards}")

    q = (94085, 73810, 57978, 41773, 23784)
    for _ in range(1525):
        q, guards = borrowed_step(q, 100001)
        assert min(guards) >= 0
    _, guards = borrowed_step(q, 100001)
    assert guards == (3043, -2, -2, 3044)
    print(f"eleven-digit first borrowed failure at step 1525, guards {guards}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("catalogue", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    check_catalogue(args.catalogue)
    check_coincidence_controls()
    cases = {
        "singular-expansion": expansion(4, "(>= x0_1 x0_2)", singular, 10),
        "singular-normal-return": check_normal_return(SINGULAR_WORD, singular_normal, singular),
        "sphere-anchor-and-energy": sphere(),
    }
    for label, word in TRIANGLE_WORDS.items():
        cases[f"{label}-normal-return"] = check_normal_return(word, triangle_normal, triangle)
    args.output.mkdir(parents=True, exist_ok=True)
    for label, lines in cases.items():
        (args.output / f"{label}.smt2").write_text(
            PREAMBLE + "\n".join(lines) + "\n(check-sat)\n(exit)\n"
        )


if __name__ == "__main__":
    main()
