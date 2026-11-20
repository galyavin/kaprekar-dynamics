"""Uncompress one period-66 ingress to obtain its original source mass."""
from argparse import Namespace
from fractions import Fraction as Q
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory

from stride_inventory import verify
from geometry import det, xchart


ROOT = (0, 6)
PREFIX = tuple((ROOT, ROOT, remainder, index, descriptor) for remainder, index, descriptor in (
    (6, 3, (0, 5, 5, 0, 5, 36, -30, -24, -18, -12, -6, 0, 0, 2, 0, 0, 0, 6)),
    (0, 2, (0, 10, 10, 0, 4, 12, -60, -48, 0, -24, -12, -2, 2, 2, 2, -1, 1, 6)),
    (36, 11, (0, 10, 10, 1, 3, 66, -60, -48, 0, -18, -12, 2, -1, -2, -1, 1, -1, -42)),
    (18, 12, (5, 120, 10, 0, 2, 144, -60, -12, 0, 12, -12, 2, 0, -2, 0, 2, 0, -102)),
    (18, 7, (30, 120, 20, 0, 2, 126, -66, -24, 0, 24, 66, -2, 3, -5, 6, -5, 0, -120)),
))


def small_representative(poly, chart, rng, retry):
    vertices = tuple(poly.vertices[i] for i in next(poly.simplices()))
    orientation = det(vertices)
    assert orientation
    weights = [rng.randrange(1, 13) if retry else 1 for _ in vertices]
    centre = tuple(sum(w * v[i] for w, v in zip(weights, vertices)) / sum(weights) for i in range(6))
    for scale in (2**20, 2**27, 2**34, 2**41):
        numerators = [round(value * scale) for value in centre[:5]]
        gaps = (*[Q(value, scale) for value in numerators], Q(scale - sum(numerators), scale))
        if any(det((*vertices[:i], gaps, *vertices[i+1:])) * orientation <= 0 for i in range(6)):
            continue
        centered = xchart(gaps)
        leading = [(sum(chart[i + j] * centered[j] for j in range(5)) + chart[i + 5]) / chart[i + 6]
                   for i in range(0, 35, 7)]
        denominator = lcm(*(v.denominator for v in leading))
        return (denominator, *(int(v * denominator) for v in leading))
    raise AssertionError('no dyadic representative in the selected interior simplex')


if __name__ == '__main__':
    verify.cover.exact_point = small_representative
    with TemporaryDirectory(prefix='task144-witness-') as temp:
        reader = verify.build(Path(temp) / 'chart', False)
        config = Namespace(offsets=[0], rank=1, compress_planes=False, saturate_planes=False,
                           compress_lines=False, shorts=500, limit=2500, keep_limit=2,
                           prefix=PREFIX, stop_signature='543221')
        assert verify.cover.run_reader(reader, config) == 2
