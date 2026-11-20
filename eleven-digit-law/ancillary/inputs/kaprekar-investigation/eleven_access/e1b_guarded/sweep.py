"""Search or replay finite guard-cone and limiting-basin certificates."""

import argparse
from collections import Counter
from gzip import compress, open as gzip_open
import json
from pathlib import Path

from anchor import constrained_counts, gap_order_profiles, majorization_guards
from covectors import add_fixed_constraints, candidate_rows, guarded_orbit, monoid_space
from covectors import verify_cover
from guarded import ExactCones, branch_word, dot, family, load_rows, normal_basis
from guarded import null_basis, open_rows, phases_for, words
from noncontraction import survivors, verify_norm
from repeated import invariant_repeat, repeatable, replay_ray


TRANSVERSE = {
    66: (6, ((1, 0, -1, 2, -1), (2, -1, 0, 0, -1), (0, 1, -2, 1, 0),
             (1, 0, 0, 1, -2), (2, -2, 1, 0, 0))),
    143: (3, ((0, 0, 1, -2, 0), (0, 2, 1, 2, -4), (3, -3, -1, -1, 3),
              (3, -2, -1, 1, 1), (3, -1, 1, -1, -1), (3, 0, 1, 1, -3))),
}


def try_cover(cones, phases, tangents, candidates, caps):
    misses = Counter()
    for u in candidates:
        cover, diagnosis = guarded_orbit(cones, phases, tangents, 0, u, *caps)
        if cover is not None:
            final_rows = verify_cover(cones, phases, tangents, 0, u, cover)
            return {"method": "covector", "seed": u, "cover": cover,
                    "return_rows": len(final_rows), "phase_rows": diagnosis[1], "cones": len(cover)}
        misses[diagnosis[0]] += 1
    return {"caps_hit": dict(misses)}


def majorization_proof(cones, phases, tangents, points, depth):
    """Cover an isolated or one-dimensional cell by its exact gap order profiles."""
    if len(tangents[0]) == 0:
        profiles = ((0, points),)
    elif len(tangents[0]) == 1:
        profiles = gap_order_profiles(tangents, points)
    else:
        return None
    witnesses = []
    for parameter, profile_points in profiles:
        constraints = []
        for q, T in zip(profile_points, tangents):
            strict, weak = [], []
            for row in majorization_guards(q):
                (strict if all(dot(row, t) == 0 for t in T) else weak).append(row)
            constraints.append((tuple(strict), tuple(weak)))
        counts, exhausted = constrained_counts(cones, phases, 0, constraints, depth)
        if not exhausted or counts[-1]:
            return None
        last = next(i for i, c in enumerate(counts) if c == 0)
        witnesses.append({"parameter": str(parameter), "counts": counts[:last+1]})
    return witnesses


def exact_transverse(cones, phases, tangents, rounds, norm):
    verify_norm(norm, tangents[0])
    counts = survivors(cones, phases, 0, norm, rounds)
    assert counts[-1] == 0
    return {"method": "noncontraction", "norm": norm, "counts": counts}


def repetition_test(cones, returns, tangents):
    columns = normal_basis(tangents[0])
    exhibits = []
    for labels, W, guards in returns:
        rays = repeatable(cones, W, guards)
        cone = None if rays else invariant_repeat(cones, W, guards, cap=3)
        if not rays and cone is None:
            continue
        for ray in rays or (cone,):
            replay_ray(cones, W, guards, ray)
            exhibits.append({"labels": labels, **ray})
        columns = add_fixed_constraints(columns, W)
        if len(columns) == 5:
            break
    return exhibits, null_basis(columns)


def search(cones, index, row, caps):
    phases, tangents = phases_for(row)
    period = len(phases)
    _, _, equations, points = family(row)
    returns = tuple(words(cones, phases, 0, period))
    dimensions = {"dimension": row["dimension"], "period": period,
                  "one_return_words": len(returns),
                  "one_return_monoid_upper_dimension": len(monoid_space(
                      [W for _, W, _ in returns], tangents[0]))}
    short_cover_first = len(returns) > 64 or period != 5
    basis = null_basis(tangents[0])
    candidates = candidate_rows(basis, tangents[0], equations)
    result = try_cover(cones, phases, tangents, candidates[:4] if short_cover_first else candidates,
                       caps)
    if "cover" in result:
        return {**dimensions, **result}

    exhibits, candidate_basis = repetition_test(cones, returns, tangents)
    metadata = {**dimensions, **result, "repeated_words": exhibits,
                "repetition_upper_dimension": len(candidate_basis)}
    if candidate_basis:
        retry = candidate_rows(candidate_basis, tangents[0], equations)
        retry = tuple(u for u in retry if not short_cover_first or u not in candidates[:4])
        if short_cover_first and retry:
            result = try_cover(cones, phases, tangents, retry, caps)
            if "cover" in result:
                return {**metadata, **result}
    profiles = majorization_proof(cones, phases, tangents, points[:period], min(2*period+2, 14))
    if profiles is not None:
        return {**metadata, "method": "majorization", "profiles": profiles}
    if index in TRANSVERSE:
        rounds, norm = TRANSVERSE[index]
        return {**metadata, **exact_transverse(cones, phases, tangents, rounds, norm)}
    return metadata


def replay(cones, index, row, evidence):
    phases, tangents = phases_for(row)
    period = len(phases)
    assert (row["dimension"], period) == (evidence["dimension"], evidence["period"])
    returns = tuple(words(cones, phases, 0, period))
    assert len(returns) == evidence["one_return_words"]
    assert len(monoid_space([W for _, W, _ in returns], tangents[0])) == (
        evidence["one_return_monoid_upper_dimension"])
    if "repeated_words" in evidence:
        columns = normal_basis(tangents[0])
        for ray in evidence["repeated_words"]:
            assert len(ray["labels"]) % period == 0
            W, guards = branch_word(phases, 0, ray["labels"])
            assert cones.feasible(guards)
            replay_ray(cones, W, guards, ray)
            columns = add_fixed_constraints(columns, W)
        assert len(null_basis(columns)) == evidence["repetition_upper_dimension"]
    if evidence.get("method") == "covector":
        final_rows = verify_cover(cones, phases, tangents, 0, evidence["seed"], evidence["cover"])
        assert len(final_rows) == evidence["return_rows"]
        assert len(evidence["cover"]) == evidence["cones"]
        assert len({(p, tuple(u)) for p, u, _ in evidence["cover"]}) == evidence["phase_rows"]
    elif evidence.get("method") == "majorization":
        depth = max(len(p["counts"])-1 for p in evidence["profiles"])
        _, _, _, points = family(row)
        assert majorization_proof(cones, phases, tangents, points[:period], depth) == evidence["profiles"]
    elif evidence.get("method") == "noncontraction":
        rounds = len(evidence["counts"]) - 1
        norm = tuple(map(tuple, evidence["norm"]))
        assert exact_transverse(cones, phases, tangents, rounds, norm)["counts"] == evidence["counts"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rows", type=int, nargs="*")
    parser.add_argument("--row-cap", type=int, default=256)
    parser.add_argument("--cone-cap", type=int, default=2048)
    parser.add_argument("--replay", type=Path)
    args = parser.parse_args()
    records = load_rows()
    ids = open_rows(records) if args.rows is None else args.rows
    cones = ExactCones()
    result = {}
    if args.replay:
        with gzip_open(args.replay, "rt") as source:
            result = json.load(source)
        assert sorted(map(int, result)) == ids
        for i in ids:
            replay(cones, i, records[i], result[str(i)])
    else:
        for i in ids:
            result[str(i)] = search(cones, i, records[i], (args.row_cap, args.cone_cap))
            row = result[str(i)]
            print(i, {k: v for k, v in row.items() if k not in ("cover", "repeated_words")},
                  f"exhibits={len(row.get('repeated_words', ()))}", flush=True)
        if args.output:
            data = (json.dumps(result, separators=(",", ":")) + "\n").encode()
            args.output.write_bytes(compress(data, mtime=0))
    if args.rows is None:
        assert all(r["one_return_monoid_upper_dimension"] == 0 for r in result.values())
        assert Counter(r.get("method") for r in result.values()) == {
            "covector": 74, "noncontraction": 2, "majorization": 9, None: 1,
        }
        assert [int(i) for i, r in result.items() if "method" not in r] == [181]
    for method in ("covector", "noncontraction", "majorization"):
        print(method, [int(i) for i, r in result.items() if r.get("method") == method])
    print("Open rows:", [int(i) for i, r in result.items() if "method" not in r])
    print("Exact queries:", cones.count)


if __name__ == "__main__":
    main()
