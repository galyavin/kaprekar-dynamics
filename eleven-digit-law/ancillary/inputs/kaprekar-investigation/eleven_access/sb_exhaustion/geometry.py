"""Exact geometry and original limiting mass in the central gap simplex."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import accumulate, combinations, product
from math import factorial, gcd, lcm

D = 6
U_CENTER = tuple(Q(a, 180000) for a in (26293, 27557, 28481, 40447, 29865, 27357))
U_HALF = Q(1, 10**6)


def dot(row, x):
    return sum(a * b for a, b in zip(row, x))


def integer_row(row):
    row = tuple(row)
    scale = lcm(*(Q(v).denominator for v in row))
    result = [int(scale * v) for v in row]
    common = gcd(*result)
    return tuple(v // common for v in result) if common else tuple(result)


def basis(i):
    return tuple(int(i == j) for j in range(D))


def xchart(gaps):
    return tuple(i - 6 * a for i, a in enumerate(accumulate(gaps[:5]), 1))


def conversion(normal):
    # h>=0, sum h=1, y=(1+h)/7, and x=7X/B on each vertex h=e_i.
    return integer_row(dot(normal, (i - 6 * (j < i) for i in range(1, D)))
                       for j in range(D))


U_CUTS = tuple(integer_row((1 - 7 * (U_CENTER[i] - U_HALF)) + (j == i)
                           if side == 1 else (7 * (U_CENTER[i] + U_HALF) - 1) - (j == i)
                           for j in range(D))
               for i, side in product(range(5), (1, -1)))


def det(matrix):
    a = [list(row) for row in matrix]
    result = Q(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Q(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        p = a[i][i]
        result *= p
        for j in range(i + 1, len(a)):
            ratio = a[j][i] / p
            for k in range(i + 1, len(a)):
                a[j][k] -= ratio * a[i][k]
    return result


@dataclass(frozen=True)
class Poly:
    vertices: tuple
    facets: tuple

    @classmethod
    def root(cls):
        return cls(tuple(tuple(Q(x) for x in basis(i)) for i in range(D)),
                   tuple(frozenset(j for j in range(D) if j != i) for i in range(D)))

    def cut(self, row):
        values = [dot(row, v) for v in self.vertices]
        if min(values) >= 0:
            return self
        if max(values) <= 0:
            return None
        vertex_facets = [frozenset(k for k, face in enumerate(self.facets) if i in face)
                         for i in range(len(values))]
        points = [v for v, value in zip(self.vertices, values) if value >= 0]
        markings = [f | ({len(self.facets)} if value == 0 else set())
                    for f, value in zip(vertex_facets, values) if value >= 0]
        for i, j in combinations(range(len(values)), 2):
            a, b = values[i], values[j]
            if a * b >= 0:
                continue
            common = vertex_facets[i] & vertex_facets[j]
            # The intersection of the common active facets defines the
            # minimal face containing these vertices. It is an edge iff
            # there are exactly two vertices on it.
            if sum(common <= other for other in vertex_facets) != 2:
                continue
            points.append(tuple((a * w - b * v) / (a - b)
                                for v, w in zip(self.vertices[i], self.vertices[j])))
            markings.append(common | {len(self.facets)})
        assert len(points) == len(set(points))
        faces = {frozenset(i for i, labels in enumerate(markings) if k in labels)
                 for k in range(len(self.facets) + 1)}
        # The maximal nonempty restricted old faces and the new cut face
        # are exactly the current facets. Any smaller old face is redundant.
        maximal = tuple(sorted((face for face in faces if face and not any(face < f for f in faces)),
                               key=lambda face: tuple(sorted(face))))
        return Poly(tuple(points), maximal)

    def restrict(self, rows):
        poly = self
        for row in rows:
            poly = poly.cut(row)
            if poly is None:
                break
        return poly

    def difference(self, rows):
        poly = self
        for row in rows:
            outside = poly.cut(tuple(-x for x in row))
            inside = poly.cut(row)
            if outside is poly and inside is poly:
                continue
            if outside is not None:
                yield outside
            if inside is None:
                return
            poly = inside

    def simplices(self):
        def pull(vertices, faces, dimension):
            if dimension == 0:
                assert len(vertices) == 1
                yield (next(iter(vertices)),)
                return
            if len(vertices) == dimension + 1:
                yield tuple(sorted(vertices))
                return
            anchor = min(vertices)
            for face in faces:
                if anchor in face:
                    continue
                subfaces = {face & other for other in faces if other != face}
                maximal = [g for g in subfaces if g and not any(g < f for f in subfaces)]
                for simplex in pull(face, maximal, dimension - 1):
                    yield (anchor, *simplex)
        yield from pull(frozenset(range(len(self.vertices))), self.facets, D - 1)

    def mass(self):
        return sum(simplex_mass(tuple(self.vertices[i] for i in s)) for s in self.simplices())

    def center(self):
        return tuple(sum(v[i] for v in self.vertices) / len(self.vertices) for i in range(D))


def simplex_mass(vertices):
    # Expand prod_i (1+h_i) in each simplex's barycentric coordinates.
    # Integer numerator coefficients keep the polynomial integration exact.
    scale = lcm(*(v.denominator for row in vertices for v in row))
    bernstein = {(0,) * D: 1}
    for j in range(D):
        values = [int(scale * (1 + row[j])) for row in vertices]
        expanded = {}
        for power, coeff in bernstein.items():
            for i, val in enumerate(values):
                key = power[:i] + (power[i] + 1,) + power[i + 1:]
                expanded[key] = expanded.get(key, 0) + val * coeff
        bernstein = expanded
    numerator = sum(coeff * factorial_product(power) for power, coeff in bernstein.items())
    jac = abs(det(vertices))
    assert jac > 0
    return jac * Q(numerator, scale**D * 7**11)


def factorial_product(power):
    result = 1
    for exponent in power:
        result *= factorial(exponent)
    return result


def u_mass():
    result = factorial(11) * (2 * U_HALF)**5
    for center in U_CENTER[:5]:
        result *= center
    return result * (U_CENTER[-1] - U_HALF**2 * sum(1 / v for v in U_CENTER[:5]) / 3)


INITIAL = ((1, 0, 0, 0, 1), (0, 1, 0, 1, 0), (0, 0, 1, 0, 0))


def entrance_cells():
    root = Poly.root()
    result = []
    for signs in product((1, -1), repeat=3):
        cuts = tuple(conversion(tuple(s * v for v in row)) for s, row in zip(signs, INITIAL))
        poly = root.restrict(cuts)
        result.append((poly.mass(), signs, poly, cuts))
    return sorted(result, reverse=True, key=lambda row: row[0])


if __name__ == '__main__':
    total = Poly.root().mass()
    print('S', total, float(total))
    assert total == Q(805597, 1977326743)
    cells = entrance_cells()
    assert sum(row[0] for row in cells) == total
    for mass, signs, poly, _ in cells:
        print('cell', *signs, mass / total, float(mass / total),
              'center', *integer_row(xchart(tuple((1 + h) / 7 for h in poly.center()))))
    uc = Poly.root().restrict(U_CUTS)
    assert uc is not None and uc.mass() == u_mass(), (uc.mass(), u_mass())
    print('U', u_mass(), float(u_mass()))
