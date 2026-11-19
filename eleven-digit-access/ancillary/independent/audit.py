"""Audit the independently generated arithmetic catalogue without an LP solver."""
import argparse
from collections import Counter
from fractions import Fraction
from math import lcm
from pathlib import Path
import gzip
import json

from exhaust import (N, ROWS, Failure, checked_certificate, dot, energy, gaps_from_q,
                     leading, q_from_gaps, sign, sorting_faces, write_json)
from work import maximum_phase_work, residual_groups


def records(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            yield json.loads(line)


def point(entries):
    return tuple(map(Fraction, entries))


def ratio(nums, denominator):
    return tuple(Fraction(x, denominator) for x in nums)


def cyclic(g):
    phases = [g]
    while True:
        g = leading(g)
        if g == phases[0]:
            return phases
        assert g not in phases and len(phases) < 720
        phases.append(g)


def geometry(word, vertices, faces_by_id):
    center = ratio((1,) * N, N)
    ray_tip = ratio((0, 1, 1, 0, 0, 1), 3)
    triangle = {
        ratio((2, 4, 2, 0, 3, 3), 14),
        ratio((1, 2, 1, 1, 1, 2), 8),
        ratio((1, 2, 1, 1, 2, 1), 8),
    }
    apex = ratio((2, 4, 2, 0, 3, 3), 14)
    isolated = ratio((1, 1, 2, 3, 1, 2), 10)
    for face_id in word:
        vertices = set(vertices)
        if vertices == {center, ray_tip}:
            return "singular ray"
        if vertices == {isolated}:
            return "isolated orbit"
        if apex in vertices and len(vertices) >= 2 and vertices <= triangle:
            return "triangle" if vertices == triangle else "triangle edge"
        perm = faces_by_id[face_id].permutation
        vertices = [tuple(v[i] for i in perm) for v in vertices]
    raise AssertionError((word, vertices))


def actual_next(q, base):
    descending = (*q, 0, 0, 0, 0, 0, 0)
    asc = descending[::-1]
    lhs = rhs = 0
    for a, b in zip(descending, asc):
        lhs = lhs * base + a
        rhs = rhs * base + b
    value = lhs - rhs
    digits = []
    for _ in range(11):
        value, digit = divmod(value, base)
        digits.append(digit)
    assert value == 0
    digits.sort(reverse=True)
    return tuple(digits[i] - digits[-1 - i] for i in range(5))


def realize_offsets(q, offsets, base):
    realized = [base * p + d for p, d in zip(q, offsets)]
    for radius, positions in residual_groups(q):
        if radius == 0:
            continue
        first = positions[0]
        first_sign = sign(q[first] - Fraction(1, 2))
        translation = first_sign * (int(realized[first]) - realized[first])
        for i in positions:
            realized[i] += translation * sign(q[i] - Fraction(1, 2))
    assert all(d.denominator == 1 for d in realized)
    result = tuple(map(int, realized))
    assert base > result[0] and result[-1] > 0
    assert all(a >= b for a, b in zip(result, result[1:]))
    return result


def physical_phase_check(gaps, w, doubled_offsets):
    q = q_from_gaps(gaps)
    step = 2 * lcm(*(p.denominator for p in q))
    first = 4096
    if Fraction(1, 2) in q:
        first += doubled_offsets[q.index(Fraction(1, 2))] % 2
    offsets = tuple(Fraction(i, 2) for i in doubled_offsets)

    def check_base(base):
        src = realize_offsets(q, offsets, base)
        dst = actual_next(src, base)
        source_gaps = gaps_from_q(src, total=base)
        target_gaps = gaps_from_q(dst, total=base)
        return dot(target_gaps, target_gaps) - dot(source_gaps, source_gaps) - base * w
    assert check_base(first) == check_base(first + step)


def audit(output):
    faces_by_id, _ = sorting_faces()
    faces_by_sign = {face.signs: face for face in faces_by_id.values()}
    summary = json.loads((output / "summary.json").read_text())
    ncertificates = ninstances = 0
    for certificate in records(output / "certificates.jsonl.gz"):
        rows = tuple(tuple(r) for r in certificate["rows"])
        failure = Failure(rows, tuple(certificate["strict"]))
        terms = [(tuple(row), Fraction(weight)) for row, weight in certificate["terms"]]
        checked_certificate(failure, terms)
        ncertificates += 1
        ninstances += certificate["count"]
    assert summary["certificates"] == {"instances": ninstances, "distinct": ncertificates}

    high = Counter()
    geometry_rows = []
    work_dimensions = Counter()
    work_vertices = set()
    total_tests = 0
    cutoff = Fraction(1, 30)
    for item in records(output / "catalogue.jsonl.gz"):
        word = tuple((tuple(map(tuple, groups)), center, edge)
                     for groups, center, edge in item["word"])
        if item["category"] == "positive" and item["at_or_above_sphere"]:
            label = geometry(word, map(point, item["vertices"]), faces_by_id)
            high[label] += 1
            erange = item["energy"]
            geometry_rows.append({"word": word, "geometry": label, "energy": erange})
        if item["ordinary"]:
            assert all(Fraction(test["constant"]) >= Fraction(2, 7)
                       for test in item["work"]["tests"] if test["constant"] is not None)
            work_dimensions[item["energy"]["dimension"]] += 1
            total_tests += len(item["work"]["tests"])
            work_vertices.update(point(row["point"]) for row in item["work"]["tests"])

    closure_cases = Counter()
    for gaps in sorted(work_vertices):
        if not 0 < gaps[0] < Fraction(1, 6):
            closure_cases["nonordinary"] += 1
            continue
        phases = cyclic(gaps)
        work = sum(maximum_phase_work(
            faces_by_sign[tuple(sign(dot(row, phase)) for row in ROWS)], phase)[0]
                   for phase in phases)
        radius = max(abs(gap - Fraction(1, N)) for gap in gaps)
        factor = Fraction(4, 7) if energy(gaps) < cutoff else Fraction(2, 7)
        assert work + factor * len(phases) * radius <= 0, (gaps, work)
        closure_cases["ordinary"] += 1
        if energy(gaps) < cutoff:
            closure_cases["ordinary_in_ball"] += 1

    target = gaps_from_q(ratio((12, 8, 6, 6, 3), 14))
    peak = []
    for gaps in cyclic(target):
        face = faces_by_sign[tuple(sign(dot(row, gaps)) for row in ROWS)]
        w, offsets = maximum_phase_work(face, gaps)
        physical_phase_check(gaps, w, offsets)
        peak.append({"q": q_from_gaps(gaps), "work": w, "doubled_offsets": offsets})
    assert [p["work"] for p in peak] == [0, Fraction(2, 7), Fraction(-3, 7)]
    write_json(output / "audit.json",
               {"certificate_checks": ncertificates, "certificate_instances": ninstances,
                "high_energy_geometry": dict(high), "high_energy_cells": geometry_rows,
                "work_dimensions": dict(work_dimensions), "work_vertex_tests": total_tests,
                "distinct_work_vertices": len(work_vertices), "closure_strata": dict(closure_cases),
                "global_sharp_orbit": peak})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    audit(parser.parse_args().output)
