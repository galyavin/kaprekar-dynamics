#!/usr/bin/env python3
"""Enumerate eleven-digit leading periodic cells by rational cone clipping."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd
from pathlib import Path
import gzip
import io
import json

N = 6
ONE = (1,) * N
ZERO = (0,) * N
UNIT = tuple(tuple(int(i == j) for i in range(N)) for j in range(N))
INITIAL = tuple(range(N))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def lin(a, b, m=1, n=-1):
    return tuple(m * x + n * y for x, y in zip(a, b))


def minus(a):
    return tuple(-x for x in a)


def primitive(row):
    scale = gcd(*row)
    return tuple(x // scale for x in row) if scale else row


def pull(row, labels):
    result = [0] * N
    for i, j in enumerate(labels):
        result[j] = row[i]
    return tuple(result)


def q_from_gaps(g):
    return tuple(sum(g[i:]) for i in range(1, N))


def gaps_from_q(q, total=1):
    extended = (total, *q, 0)
    return tuple(a - b for a, b in zip(extended, extended[1:]))


def leading(g):
    q = q_from_gaps(g)
    digits = sorted((1, *q, *(1 - x for x in q)), reverse=True)
    return gaps_from_q(tuple(digits[i] - digits[-1 - i] for i in range(5)))


def representative(groups, central, edge):
    radii = tuple(Fraction(len(groups) - i, 2 * (len(groups) + (not edge)))
                  for i in range(len(groups)))
    return (tuple(Fraction(1, 2) + r for r, (hi, _) in zip(radii, groups) for _ in range(hi))
            + (Fraction(1, 2),) * central
            + tuple(Fraction(1, 2) - r for r, (_, lo) in reversed(tuple(zip(radii, groups)))
                    for _ in range(lo)))


def group_lists(size):
    if size == 0:
        yield ()
    for first in range(1, size + 1):
        for hi in range(first + 1):
            for rest in group_lists(size - first):
                yield ((hi, first - hi), *rest)


def circle_zero_loss(q):
    multiplicities = defaultdict(lambda: [0, 0])
    for coord in (0, *q):
        multiplicities[coord % 1][0] += 1
        multiplicities[-coord % 1][1] += 1
    sites = [multiplicities[coord] for coord in sorted(multiplicities)]
    colors = [None if a and b else int(bool(b)) for a, b in sites]
    if any(a + b > 1 and color is not None for (a, b), color in zip(sites, colors)):
        return False
    return all(x is None or x != y for x, y in zip(colors, sites_wrap(colors)))


def sites_wrap(sequence):
    return sequence[1:] + sequence[:1]


def symbolic_permutation(groups, central, edge, q):
    # Each form is twice the indicated normalized coordinate; index 0 is constant.
    radii = [UNIT[0] if edge and i == 0 else UNIT[i + (not edge)]
             for i in range(len(groups))]
    forms = (tuple(lin(UNIT[0], r, n=1) for r, (hi, _) in zip(radii, groups)
                   for _ in range(hi))
             + (UNIT[0],) * central
             + tuple(lin(UNIT[0], r) for r, (_, lo) in reversed(tuple(zip(radii, groups)))
                     for _ in range(lo)))
    total = lin(UNIT[0], UNIT[0], n=1)
    order = sorted(zip((1, *q, *(1 - x for x in q)),
                       (total, *forms, *(lin(total, x) for x in forms))),
                   key=lambda digit: digit[0], reverse=True)
    digits = [x[1] for x in order]
    output_q = [lin(digits[i], digits[-1 - i]) for i in range(5)]
    source_coords = (total, *forms, ZERO)
    target_coords = (total, *output_q, ZERO)
    source = [lin(a, b) for a, b in zip(source_coords, source_coords[1:])]
    target = [lin(a, b) for a, b in zip(target_coords, target_coords[1:])]
    remaining = list(enumerate(source))
    permutation = []
    for form in target:
        matches = [i for i, (_, other) in enumerate(remaining) if form == other]
        if not matches:
            return None
        label, _ = remaining.pop(matches[0])
        permutation.append(label)
    return tuple(permutation)


def all_guard_rows():
    q = tuple(tuple(int(j >= i) for j in range(N)) for i in range(1, N))
    return tuple(dict.fromkeys([
        UNIT[0], UNIT[-1],
        *(primitive(lin(a, b)) for a, b in combinations(q, 2)),
        *(primitive(lin(lin(a, b, n=1), ONE)) for i, a in enumerate(q) for b in q[i:]),
    ]))


ROWS = all_guard_rows()


def sign(x):
    return (x > 0) - (x < 0)


@dataclass(frozen=True)
class SortingFace:
    radii: tuple
    central: int
    edge: bool
    signs: tuple
    permutation: tuple
    exemplar: tuple

    @property
    def ID(self):
        return self.radii, self.central, self.edge


def sorting_faces():
    counts = Counter()
    faces = {}
    signatures = set()
    for central in range(6):
        for groups in group_lists(5 - central):
            for edge in ([False, True] if groups else [False]):
                q = representative(groups, central, edge)
                g = gaps_from_q(q)
                signs = tuple(sign(dot(row, g)) for row in ROWS)
                assert signs not in signatures
                signatures.add(signs)
                counts["all"] += 1
                permutation = symbolic_permutation(groups, central, edge, q)
                lossless = circle_zero_loss(q)
                assert lossless == (permutation is not None)
                if lossless:
                    counts["zero_loss"] += 1
                    assert leading(g) == tuple(g[i] for i in permutation)
                    if permutation[0] == 0:
                        counts["anchor_stationary"] += 1
                        face = SortingFace(groups, central, edge, signs, permutation, g)
                        faces[face.ID] = face
    return faces, dict(counts)


@dataclass(frozen=True)
class Trie:
    rules: tuple
    split: int = None
    children: tuple = ()
    face: object = None


def build_trie(faces, available=tuple(range(len(ROWS)))):
    constant = tuple((i, faces[0].signs[i]) for i in available
                     if all(face.signs[i] == faces[0].signs[i] for face in faces))
    if len(faces) == 1:
        return Trie(constant, face=faces[0])
    remaining = tuple(i for i in available if all(i != row for row, _ in constant))
    split = min(remaining, key=lambda i: max(Counter(face.signs[i] for face in faces).values()))
    children = []
    for branch in (0, 1, -1):
        selected = tuple(face for face in faces if face.signs[split] == branch)
        if selected:
            children.append((branch, build_trie(selected, tuple(i for i in remaining if i != split))))
    return Trie(constant, split, tuple(children))


@dataclass(frozen=True)
class Cone:
    rows: tuple
    strict: frozenset
    rays: tuple

    def interior_sample(self):
        ray_sum = tuple(map(sum, zip(*self.rays)))
        return tuple(Fraction(x, sum(ray_sum)) for x in ray_sum)

    def vertices(self):
        return tuple(tuple(Fraction(x, sum(v)) for x in v) for v in self.rays)


SIMPLEX = Cone(UNIT, frozenset(), UNIT)


@lru_cache(maxsize=30000)
def dot_values(row, rays):
    return tuple(dot(row, ray) for ray in rays)


@lru_cache(maxsize=30000)
def active_masks(rows, rays):
    masks = [0] * len(rays)
    for i, row in enumerate(rows):
        for j, value in enumerate(dot_values(row, rays)):
            if value == 0:
                masks[j] |= 1 << i
    return tuple(masks)


def ray_clipping(cone, row):
    heights = dot_values(row, cone.rays)
    negative = [i for i, h in enumerate(heights) if h < 0]
    if not negative:
        return cone.rays
    keep = [r for r, h in zip(cone.rays, heights) if h >= 0]
    positive = [i for i, h in enumerate(heights) if h > 0]
    if not positive:
        return tuple(keep)
    masks = active_masks(cone.rows, cone.rays)
    for i in positive:
        for j in negative:
            common = masks[i] & masks[j]
            if any(k != i and k != j and mask & common == common
                   for k, mask in enumerate(masks)):
                continue
            keep.append(primitive(lin(cone.rays[i], cone.rays[j], -heights[j], heights[i])))
    return tuple(sorted(set(keep)))


@dataclass(frozen=True)
class Failure:
    rows: tuple
    strict_row: tuple


def impose(cone, row, strict=False, equality=False):
    row = primitive(row)
    rays = ray_clipping(cone, row) if row not in cone.rows else cone.rays
    new_rows = tuple(sorted(set((*cone.rows, row))))
    new_strict = cone.strict | {row} if strict else cone.strict
    if equality:
        reverse = minus(row)
        rays = tuple(ray for ray in rays if dot(row, ray) == 0)
        new_rows = tuple(sorted(set((*new_rows, reverse))))
    if not rays:
        return Failure(new_rows, ONE)
    for required in sorted(new_strict):
        if not any(dot_values(required, rays)):
            # Any positive dual weight must be on a row vanishing on the cone.
            active_rows = tuple(r for r in new_rows if not any(dot_values(r, rays)))
            return Failure(active_rows, required)
    return Cone(new_rows, new_strict, rays)


def restrict(cone, row, side):
    if side == 0:
        return impose(cone, row, equality=True)
    return impose(cone, row if side > 0 else minus(row), strict=True)


def rotations(word):
    return [word[i:] + word[:i] for i in range(len(word))]


def align(cone, labels):
    for i, j in enumerate(labels):
        if i != j:
            cone = impose(cone, lin(UNIT[i], UNIT[j]), equality=True)
            if isinstance(cone, Failure):
                return cone
    return cone


class Search:
    def __init__(self):
        self.prunings = Counter()
        self.counts = Counter()
        self.words = {}

    def check(self, result, operation):
        self.counts[operation] += 1
        if isinstance(result, Failure):
            self.prunings[result] += 1
            self.counts[operation + "_rejected"] += 1
            return None
        return result

    def run(self, cone, trie):
        self.visit(cone, trie, (INITIAL,), ())

    def visit(self, cone, trie, labels, word):
        current = labels[-1]
        for index, side in trie.rules:
            cone = self.check(restrict(cone, pull(ROWS[index], current), side), "guard")
            if cone is None:
                return
        if trie.face is None:
            row = pull(ROWS[trie.split], current)
            for side, child in trie.children:
                branch = self.check(restrict(cone, row, side), "guard")
                if branch is not None:
                    self.visit(branch, child, labels, word)
            return
        face = trie.face
        next_labels = tuple(current[i] for i in face.permutation)
        next_word = (*word, face.ID)
        if next_labels not in labels:
            self.visit(cone, self.root, (*labels, next_labels), next_word)
            return
        self.counts["repeated_labels"] += 1
        cycle_start = labels.index(next_labels)
        period = len(next_word) - cycle_start
        self.counts["max_path"] = max(self.counts["max_path"], len(next_word))
        periodic = self.check(align(cone, labels[period] if period < len(labels) else next_labels),
                              "return")
        if periodic is None:
            return
        self.counts["periodic_leaves"] += 1
        self.words.setdefault(next_word[:period], periodic.interior_sample())


def cell_for_word(word, faces, found_sample, search):
    assert all(x >= 0 for x in found_sample)
    cone = restrict(SIMPLEX, UNIT[0], sign(found_sample[0]))
    labels = INITIAL
    for face_id in word:
        face = faces[face_id]
        for index, side in enumerate(face.signs):
            cone = search.check(restrict(cone, pull(ROWS[index], labels), side), "cell")
            assert cone is not None
        labels = tuple(labels[i] for i in face.permutation)
    cone = search.check(align(cone, labels), "cell_return")
    assert cone is not None
    assert all(dot(row, found_sample) >= 0 for row in cone.rows)
    assert all(dot(row, found_sample) > 0 for row in cone.strict)
    g = found_sample
    for face_id in word:
        assert tuple(sign(dot(row, g)) for row in ROWS) == faces[face_id].signs
        g = leading(g)
    assert g == found_sample
    return cone


def linear_system(columns, target):
    """Solve a full-column-rank linear system, or return None."""
    n = len(columns)
    rows = [[Fraction(x) for x in (*entries, value)]
            for entries, value in zip(zip(*columns), target)]
    if not columns:
        return () if not any(target) else None
    for j in range(n):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x / divisor for x in rows[j]]
        for i, row in enumerate(rows):
            if i != j and row[j]:
                multiplier = row[j]
                rows[i] = [x - multiplier * y for x, y in zip(row, rows[j])]
    if any(row[-1] for row in rows[n:]):
        return None
    return tuple(row[-1] for row in rows[:n])


def checked_certificate(failure, terms):
    assert all(row in failure.rows and weight >= 0 for row, weight in terms)
    assert all(sum(weight * row[i] for row, weight in terms) == -x
               for i, x in enumerate(failure.strict_row))


def certificate(failure, optimizer):
    target = minus(failure.strict_row)
    rows = tuple(row for row in failure.rows if row not in (failure.strict_row, ZERO))
    if target in rows:
        terms = ((target, Fraction(1)),)
    elif target == ZERO:
        terms = ()
    else:
        terms = optimizer(rows, target)
    checked_certificate(failure, terms)
    return terms


def scipy_optimizer(rows, target):
    from scipy.optimize import linprog
    result = linprog([1] * len(rows), A_eq=list(zip(*rows)), b_eq=target, method="highs-ds")
    assert result.success, (rows, target, result.message)
    support = tuple(row for row, value in zip(rows, result.x) if value > 1e-9)
    values = linear_system(support, target)
    if values is not None and all(value >= 0 for value in values):
        return tuple(zip(support, values))
    # A floating point optimizer only proposes supports. Try exact supports if needed.
    for size in range(1, N + 1):
        for support in combinations(rows, size):
            values = linear_system(support, target)
            if values is not None and all(value >= 0 for value in values):
                return tuple(zip(support, values))
    raise ArithmeticError((rows, target))


def energy(g):
    return dot(g, g) - Fraction(1, N)


def project_affine(vertices):
    origin = vertices[0]
    basis = []
    for vertex in vertices[1:]:
        difference = lin(vertex, origin)
        if linear_system((*basis, difference), ZERO) is not None:
            basis.append(difference)
    gram = tuple(tuple(dot(u, v) for v in basis) for u in basis)
    coordinates = linear_system(gram, tuple(-dot(u, origin) for u in basis))
    return tuple(o + sum(c * d[j] for c, d in zip(coordinates, basis))
                 for j, o in enumerate(origin))


def energy_range(cone):
    vertices = cone.vertices()
    whole = frozenset(range(len(vertices)))
    hyperfaces = {frozenset(i for i, vertex in enumerate(vertices) if dot(row, vertex) == 0)
                  for row in cone.rows} - {frozenset(), whole}
    faces = {whole}
    for hyperface in hyperfaces:
        faces |= {face & hyperface for face in faces if face & hyperface}
    projections = [project_affine(tuple(vertices[i] for i in sorted(face))) for face in faces]
    feasible = [p for p in projections if all(dot(row, p) >= 0 for row in cone.rows)]
    minimum = min(feasible, key=energy)
    maximum = max(map(energy, vertices))
    high = tuple(v for v in vertices if energy(v) == maximum)
    attained = lambda p: all(dot(r, p) > 0 for r in cone.strict)
    dimension = max(len(face) for face in faces)
    # Affine dimension is rank of the vertex differences, not vertex count.
    if dimension > 1:
        basis = []
        for vertex in vertices[1:]:
            difference = lin(vertex, vertices[0])
            if linear_system((*basis, difference), ZERO) is not None:
                basis.append(difference)
        dimension = len(basis)
    else:
        dimension = 0
    return {
        "min": energy(minimum),
        "min_attained": attained(minimum),
        "min_point": minimum,
        "max": maximum,
        "max_attained": any(map(attained, high)),
        "max_points": high,
        "dimension": dimension,
    }


def encode(obj):
    if isinstance(obj, Fraction):
        return str(obj)
    raise TypeError(type(obj))


def write_json(path, obj):
    with path.open("w") as stream:
        print(json.dumps(obj, default=encode, separators=(",", ":"), sort_keys=True), file=stream)


def gzip_lines(path, records):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped) as stream:
                for record in records:
                    print(json.dumps(record, default=encode, separators=(",", ":"),
                                     sort_keys=True), file=stream)


def certify_prunings(searches, output, frozen):
    combined = Counter()
    for search in searches:
        combined.update(search.prunings)
    print(f"Certifying {sum(combined.values())} prunings using {len(combined)} distinct systems.",
          flush=True)

    def records():
        reference = gzip.open(frozen, "rt") if frozen else None
        try:
            yield from checked_records(reference)
        finally:
            if reference is not None:
                reference.close()

    def checked_records(reference):
        for i, failure in enumerate(sorted(combined, key=lambda f: (f.strict_row, f.rows)), 1):
            if reference is None:
                terms = certificate(failure, scipy_optimizer)
            else:
                expected = json.loads(reference.readline())
                assert expected["rows"] == [list(r) for r in failure.rows]
                assert expected["strict"] == list(failure.strict_row)
                assert expected["count"] == combined[failure]
                terms = tuple((tuple(row), Fraction(weight)) for row, weight in expected["terms"])
                checked_certificate(failure, terms)
            if i % 20000 == 0:
                print(f"Checked {i} dual systems.", flush=True)
            yield {"strict": failure.strict_row, "rows": failure.rows,
                   "terms": terms, "count": combined[failure]}
        if reference is not None:
            assert not reference.readline()
    gzip_lines(output / "certificates.jsonl.gz", records())
    return {"instances": sum(combined.values()), "distinct": len(combined)}


def exhaust(output, frozen):
    from work import analyze_work

    faces, counts = sorting_faces()
    print(f"Sorting faces: {counts}", flush=True)
    searches = []
    for anchor_sign in (1, 0):
        search = Search()
        cone = restrict(SIMPLEX, UNIT[0], anchor_sign)
        search.root = build_trie(tuple(f for f in faces.values() if f.signs[0] == anchor_sign))
        search.run(cone, search.root)
        searches.append(search)
        print(f"Anchor sign {anchor_sign}: {len(search.words)} oriented words, "
              f"search counts {dict(search.counts)}", flush=True)
    summary = {"faces": counts, "searches": [dict(s.counts) for s in searches],
               "oriented": [len(s.words) for s in searches]}
    catalogue = []
    for search in searches:
        for word in sorted(search.words):
            if word != min(rotations(word)):
                continue
            witness = search.words[word]
            cone = cell_for_word(word, faces, witness, search)
            erange = energy_range(cone)
            cutoff = Fraction(1, 30)
            anchor_cut = lin(ONE, UNIT[0], n=-6)
            inside = search.check(impose(cone, anchor_cut, strict=True), "anchor_cut")
            category = ("zero_anchor" if witness[0] == 0 else
                        "positive" if all(witness) else "zero_other")
            item = {"word": word, "witness": witness, "category": category, "energy": erange,
                    "vertices": cone.vertices(), "in_ball": erange["min"] < cutoff,
                    "at_or_above_sphere": (erange["max"] > cutoff or
                                           erange["max"] == cutoff and erange["max_attained"]),
                    "ordinary": inside is not None and witness[0] > 0}
            if inside is not None and witness[0] > 0:
                below = energy_range(inside)
                item["ordinary_in_ball"] = below["min"] < cutoff
                item["work"] = analyze_work(word, faces, inside.vertices())
            else:
                item["ordinary_in_ball"] = False
            catalogue.append(item)
    gzip_lines(output / "catalogue.jsonl.gz", catalogue)
    gzip_lines(output / "oriented.jsonl.gz",
               ({"word": w, "witness": s.words[w]} for s in searches for w in sorted(s.words)))
    summary["categories"] = dict(Counter(item["category"] for item in catalogue))
    for key in ("in_ball", "ordinary", "ordinary_in_ball", "at_or_above_sphere"):
        summary[key] = dict(Counter(item["category"] for item in catalogue if item[key]))
    summary["interior_sphere_suprema"] = dict(Counter(
        "attained" if item["energy"]["max_attained"] else "open"
        for item in catalogue
        if item["category"] == "positive" and item["energy"]["max"] == Fraction(1, 30)))
    constants = {k: [item["work"]["constant"] for item in catalogue
                     if item.get(k) and item.get("work")
                     and item["work"]["constant"] is not None]
                 for k in ("ordinary", "ordinary_in_ball", "in_ball")}
    summary["constants"] = {k: min(values) for k, values in constants.items() if values}
    summary["certificates"] = certify_prunings(searches, output, frozen)
    summary["searches"] = [dict(search.counts) for search in searches]
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, default=encode, sort_keys=True), flush=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--frozen-certificates", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    exhaust(args.output, args.frozen_certificates)
