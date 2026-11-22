"""Find integer witnesses in a disposable copy, then replay with verify.py."""

import argparse
import json
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from verify import (EXAMPLES, LIMIT, SHORT_LIMIT, Length, cone_certificate,
                    distinct_normals, dot, echelon, kernel, phase_zero_rows,
                    primitive, relation_value, verify_faces)


def strict_point(rows, width):
    if not rows:
        return (0,) * width
    if not width or any(not any(row) for row in rows):
        return None
    ones = (1,) * width
    if all(dot(row, ones) > 0 for row in rows):
        return ones
    if width == 1:
        return (-1,) if all(row[0] < 0 for row in rows) else None
    result = linprog(np.zeros(width), A_ub=-np.array(rows, dtype=float),
                     b_ub=-np.ones(len(rows)), bounds=(None, None), method="highs")
    if result.success:
        for scale in (1, 100, 10_000, 1_000_000):
            point = tuple(round(float(v) * scale) for v in result.x)
            if all(dot(row, point) > 0 for row in rows):
                return primitive(point)
    return None


def exact_solve(equations, rhs, guess):
    """Use the proposed support. Solve pivot values exactly and check."""
    width = len(guess)
    basis = echelon([(*row, value) for row, value in zip(equations, rhs)])
    assert width not in basis, "inconsistent proposed support"
    point = [Fraction(float(v)).limit_denominator(100_000) for v in guess]
    for index, row in sorted(basis.items(), reverse=True):
        point[index] = (row[-1] - dot(row[:-1], point) + row[index] * point[index]) / row[index]
    assert all(dot(row, point) == value for row, value in zip(equations, rhs))
    return point


def separating_relation(rows):
    """Find a nonzero nonnegative dependence, then verify it over Q."""
    for i, row in enumerate(rows):
        if not any(row):
            return [(i, 1)]
    if len(rows[0]) == 1:
        i = next(i for i, row in enumerate(rows) if row[0] > 0)
        j = next(j for j, row in enumerate(rows) if row[0] < 0)
        return [(i, -rows[j][0]), (j, rows[i][0])]
    width = len(rows[0])
    equations = (*zip(*rows), (1,) * len(rows))
    rhs = (0,) * width + (1,)
    result = linprog(np.zeros(len(rows)), A_eq=np.array(equations, dtype=float),
                     b_eq=rhs, bounds=(0, None), method="highs")
    assert result.success, "no exact certificate proposed"
    support = [i for i, value in enumerate(result.x) if value > 1e-9]
    submatrix = [tuple(row[i] for i in support) for row in equations]
    weights = primitive(exact_solve(submatrix, rhs, result.x[support]))
    assert weights and all(weight > 0 for weight in weights)
    assert not any(relation_value([rows[i] for i in support], weights))
    return list(zip(support, weights))


def positive_relation(rows):
    """Find a homogeneous relation with positive weight on every row."""
    if not rows:
        return ()
    if not any(map(any, rows)):
        return (1,) * len(rows)
    equations = list(zip(*rows))
    # Solve for the extra coefficients after assigning one to every row.
    rhs = tuple(-sum(row) for row in equations)
    if not any(rhs):
        return (1,) * len(rows)
    result = linprog(np.zeros(len(rows)), A_eq=np.array(equations, dtype=float),
                     b_eq=rhs, bounds=(0, None), method="highs")
    assert result.success, "claimed implicit rows have no positive relation proposed"
    support = [i for i, value in enumerate(result.x) if value > 1e-9]
    submatrix = [tuple(row[i] for i in support) for row in equations]
    extra = exact_solve(submatrix, rhs, result.x[support])
    weights = [1] * len(rows)
    for j, w in zip(support, extra):
        weights[j] += w
    weights = primitive(weights)
    assert all(w > 0 for w in weights) and not any(relation_value(rows, weights))
    return weights


def restrict(rows, columns):
    return [primitive(dot(row, column) for column in columns) for row in rows]


def relative_interior(rows, width):
    """Identify active rows, then certify exact point in their common kernel."""
    direct = strict_point(rows, width)
    if direct is not None:
        return direct, []
    if not width:
        return (), list(range(len(rows)))
    count = len(rows)
    result = linprog(np.concatenate([np.zeros(width), -np.ones(count)]),
                     A_ub=np.hstack([-np.array(rows, dtype=float), np.eye(count)]),
                     b_ub=np.zeros(count), bounds=[(None, None)] * width + [(0, 1)] * count,
                     method="highs")
    assert result.success
    active = [i for i, slack in enumerate(result.x[width:]) if slack < 0.5]
    columns = kernel([rows[i] for i in active], width)
    remaining = [row for i, row in enumerate(rows) if i not in active]
    in_kernel = strict_point(restrict(remaining, columns), len(columns))
    assert in_kernel is not None, "no exact relative interior point proposed"
    point = primitive(dot(in_kernel, row) for row in zip(*columns)) if columns else (0,) * width
    assert all(dot(row, point) == 0 if i in active else dot(row, point) > 0
               for i, row in enumerate(rows))
    return point, active


def find_faces(length, implicit, span):
    normals = distinct_normals(length, implicit, span)
    certificates = {}
    for i, normal in phase_zero_rows(length, implicit, span).items():
        others = [row for row in normals if row != normal]
        face_columns = kernel([normal], len(span))
        # Keep the exact scaling to recover the relation in the original dual coordinates.
        rational_rows = [tuple(dot(row, col) for col in face_columns) for row in others]
        denominator = lcm(*(value.denominator for row in rational_rows for value in row))
        face_rows = [tuple(int(denominator * a) for a in row) for row in rational_rows]
        projected, active = relative_interior(face_rows, len(face_columns))
        point = (primitive(dot(projected, row) for row in zip(*face_columns))
                 if face_columns else (0,) * len(span))
        weights = positive_relation([face_rows[j] for j in active])
        record = {"point": point, "weights": weights}
        if weights:
            summed = relation_value([others[j] for j in active], weights)
            pivot = next(j for j, v in enumerate(normal) if v)
            scale = Fraction(summed[pivot], normal[pivot])
            record["scale"] = scale.numerator
            record["weights"] = tuple(w * scale.denominator for w in weights)
        certificates[str(i)] = record
    verify_faces(length, implicit, span, certificates)
    return certificates


def short_certificates():
    for n in sorted(set(range(7, SHORT_LIMIT + 1, 2)) | set(EXAMPLES)):
        length = Length(n)
        rows = [row for _, row in length.invariant_rows().values()]
        point, active = relative_interior(rows, length.s)
        record = {"interior": point, "relation": positive_relation([rows[j] for j in active])}
        implicit, span = cone_certificate(length, record)
        if n <= SHORT_LIMIT and len(span) < length.q:
            record["faces"] = find_faces(length, implicit, span)
        yield n, record


def census(limit):
    for n in range(7, limit, 2):
        length = Length(n)
        if not length.q or length.labels[3] != 0:
            continue
        retained = {i: row for i, (keep, row) in length.invariant_rows().items() if keep}
        rows = list(retained.values())
        point = strict_point(rows, length.s)
        if point is not None:
            yield [n, "+", point]
        else:
            indices = list(retained)
            yield [n, "-", [(indices[i], w) for i, w in separating_relation(rows)]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("short", "census"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--limit", type=int, default=LIMIT)
    args = parser.parse_args()
    if args.mode == "short":
        args.output.write_text(json.dumps(dict(short_certificates()), separators=(",", ":")) + "\n")
    else:
        with args.output.open("w") as out:
            for record in census(args.limit):
                out.write(json.dumps(record, separators=(",", ":")) + "\n")
