"""Exact full-dimensional cuts with edge incidence on the five-simplex."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import lcm

from pi_geometry import Poly as ArchivedPoly, dot


@lru_cache(maxsize=512)
def integer_vertices(vertices):
    result = []
    for vertex in vertices:
        scale = lcm(*(v.denominator for v in vertex))
        result.append(tuple(int(scale * v) for v in vertex))
    return tuple(result)


@lru_cache(maxsize=256)
def incidences(facets, count):
    labels = tuple(frozenset(k for k, facet in enumerate(facets) if i in facet)
                   for i in range(count))
    faces = defaultdict(list)
    for i, active in enumerate(labels):
        for shared in combinations(sorted(active), 4):
            faces[shared].append(i)
    # In dimension five an edge is cut out by four independent facet
    # equalities. A common facet set with exactly two vertices gives an edge,
    # because every nonempty face is the convex hull of its vertices.
    edges = tuple(sorted({tuple(on_face) for on_face in faces.values() if len(on_face) == 2}))
    return labels, edges


class Poly(ArchivedPoly):
    def cut(self, row):
        if not any(row):
            return self
        scale = lcm(*(getattr(a, "denominator", 1) for a in row))
        normal = tuple(int(scale * a) for a in row)
        values = tuple(dot(normal, vertex) for vertex in integer_vertices(self.vertices))
        if min(values) >= 0:
            return self
        if max(values) <= 0:
            return None

        labels, edges = incidences(self.facets, len(values))
        points = [v for v, value in zip(self.vertices, values) if value >= 0]
        new_face = len(self.facets)
        marks = [f | ({new_face} if value == 0 else set())
                 for f, value in zip(labels, values) if value >= 0]
        for i, j in edges:
            if values[i] * values[j] >= 0:
                continue
            v, w = self.vertices[i], self.vertices[j]
            a, b = dot(normal, v), dot(normal, w)
            intersection = tuple((a * y - b * x) / (a - b) for x, y in zip(v, w))
            points.append(intersection)
            marks.append((labels[i] & labels[j]) | {new_face})
        assert len(points) == len(set(points))
        faces = set()
        for k in range(len(self.facets) + 1):
            restricted = frozenset(i for i, label in enumerate(marks) if k in label)
            if restricted:
                faces.add(restricted)
        maximal = tuple(sorted((f for f in faces if not any(f < other for other in faces)),
                               key=lambda f: tuple(sorted(f))))
        return Poly(tuple(points), maximal)
