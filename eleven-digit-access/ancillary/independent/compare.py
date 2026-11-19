"""Match the frozen independent catalogue to the original rational records."""
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path
import gzip
import json

from audit import point, records
from exhaust import ROWS, dot, gaps_from_q, leading, rotations, sign, sorting_faces, write_json
from work import analyze_work


def own_word(word):
    return tuple((tuple(map(tuple, groups)), center, edge) for groups, center, edge in word)


def old_word(word):
    return tuple((tuple(map(tuple, groups)), center, edge) for center, groups, edge in word)


def archived(evidence, stem):
    with gzip.open(evidence / f"{stem}.json.gz", "rt") as stream:
        return json.load(stream)


def check_witness(g, word, faces):
    start = g
    for face_id in word:
        actual = tuple(sign(dot(row, g)) for row in ROWS)
        assert actual == faces[face_id].signs, (g, face_id, actual)
        g = leading(g)
    assert g == start


def change_phase(g, word, shift, faces):
    for face_id in word[:shift]:
        g = tuple(g[i] for i in faces[face_id].permutation)
    return g


def compare(frozen, evidence, output):
    faces, _ = sorting_faces()
    results = {"checks": {}, "differences": [],
               "vertex_counts": {"recorded": 0, "independent": 0}}
    counts = Counter()

    def check(label, identifier, observed, recorded):
        counts[label] += 1
        if observed != recorded:
            results["differences"].append({"field": label, "word": identifier,
                                           "independent": observed, "recorded": recorded})

    independent_oriented = {own_word(row["word"]): point(row["witness"])
                            for row in records(frozen / "oriented.jsonl.gz")}
    for stem in ("positive", "zero"):
        old = {old_word(r["faces"]): gaps_from_q(point(r["q"])) for r in archived(evidence, stem)}
        own = {w: g for w, g in independent_oriented.items() if bool(g[0]) == (stem == "positive")}
        check(stem + "_words_missing", None, sorted(old.keys() - own.keys()), [])
        check(stem + "_words_extra", None, sorted(own.keys() - old.keys()), [])
        for word, g in old.items():
            check_witness(g, word, faces)
            counts[stem + "_witness_replays"] += 1
    independent = {own_word(row["word"]): row for row in records(frozen / "catalogue.jsonl.gz")}
    location = {"interior": "positive", "zero-anchor": "zero_anchor", "boundary": "zero_other"}
    old_energy = archived(evidence, "energy")
    check("energy_family_keys", None, sorted(independent),
          sorted(min(rotations(old_word(r["faces"]))) for r in old_energy))
    for r in old_energy:
        word = old_word(r["faces"])
        canonical = min(rotations(word))
        if canonical not in independent:
            continue
        item = independent[canonical]
        erange = item["energy"]
        shift = rotations(word).index(canonical)
        recorded_max = sorted(change_phase(point(g), word, shift, faces)
                              for g in r["maximizing_gaps"])
        comparisons = {
            "location": (item["category"], location[r["location"]]),
            "energy_infimum": (erange["min"], r["energy_inf"]),
            "energy_supremum": (erange["max"], r["energy_sup"]),
            "infimum_attainment": (erange["min_attained"], r["inf_attained"]),
            "supremum_attainment": (erange["max_attained"], r["sup_attained"]),
            "ball_membership": (item["in_ball"], r["in_ball"]),
            "sphere_or_exterior": (item["at_or_above_sphere"], r["meets_high_energy"]),
            "energy_maximizing_points": (sorted(point(g) for g in erange["max_points"]), recorded_max)
        }
        for label, (a, b) in comparisons.items():
            check(label, word, a, b)
    old_work = archived(evidence, "work")
    old_work_keys = {min(rotations(old_word(r["faces"]))) for r in old_work}
    own_work = {w: row for w, row in independent.items() if row["ordinary"]}
    check("work_family_keys", None, sorted(own_work), sorted(old_work_keys))
    for r in old_work:
        word = old_word(r["faces"])
        canonical = min(rotations(word))
        if canonical not in own_work:
            continue
        item = own_work[canonical]
        g = point(r["point"])
        from_recorded_point = analyze_work(word, faces, (g,))["tests"][0]
        check("work_ratio", word, Fraction(item["work"]["constant"]), -Fraction(r["ratio"]))
        check("archive_point_work", word, from_recorded_point["sum"], Fraction(r["work"]))
        check("work_vertex_count", word, len(item["work"]["tests"]), r["vertices"])
        results["vertex_counts"]["recorded"] += r["vertices"]
        results["vertex_counts"]["independent"] += len(item["work"]["tests"])
    results["checks"] = dict(counts)
    write_json(output, results)
    return not results["differences"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("frozen", type=Path)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(not compare(args.frozen, args.evidence, args.output))
