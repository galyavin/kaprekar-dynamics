"""Recover P4, P8 and isolated-wall weights from their physical certificates."""

from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import json
import sys


INVESTIGATION = Path(__file__).resolve().parents[2]


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def linear(a, x, b, y):
    return tuple(a * u + b * v for u, v in zip(x, y))


def control_p4():
    sys.path.insert(0, str(INVESTIGATION))
    import nine_kernel as compiler

    def split(state, row, speed, threshold, drift, stride, base):
        slope, point, basis = state
        h, base_term = row[:4], row[4]
        # Solve s*t + h*Lz = e - h*point - base_term*B, with t = -N.
        # Maintain the exact pushforward of the generators through (t,z) -> Lz - t*d.
        divisor = speed
        generator = tuple(-x for x in drift)
        kernel = []
        for col in zip(*basis):
            value = dot(h, col)
            common = gcd(divisor, value)
            modulus = divisor // common
            v = pow(value // common, -1, modulus) if modulus > 1 else 0
            u = (common - v * value) // divisor
            assert u * divisor + v * value == common
            kernel.append(linear(-value // common, generator, modulus, col))
            generator = linear(u, generator, v, col)
            divisor = common
        assert divisor == gcd(speed, *(dot(h, col) for col in zip(*basis)))
        outgoing_basis = tuple(tuple(col[i] for col in kernel) for i in range(4))
        h_slope = dot(h, slope) + base_term * stride
        h_point = dot(h, point) + base_term * base
        assert h_slope % divisor == 0
        outgoing_slope = linear(1, slope, -h_slope // divisor, generator)
        for e in range(threshold - speed, threshold):
            if (e - h_point) % divisor == 0:
                outgoing_point = linear(1, point, (e - h_point) // divisor, generator)
                # Reuse only the compiler's formatting of the exact outgoing lattice.
                child = compiler.canonical(outgoing_slope, outgoing_point, outgoing_basis)
                assert compiler.constant_form(child, row, stride, base) == e
                yield e, Fraction(divisor, speed), child

    compiler.return_split = split
    route = ('A0', 'A:B1', 'AB')
    centers = json.loads((INVESTIGATION / 'nine_core_test_points.json').read_text())['centers']
    radius, = (rec['radius'] for rec in centers if tuple(rec['route']) == route)
    center = compiler.core.C(tuple(map(Fraction, radius)))
    saved = json.loads((INVESTIGATION / 'nine_kernel_results.json').read_text())['kernels']
    checked = 0
    for record in saved:
        name = tuple(record['route'])
        if name != route:
            continue
        base = record['base_residue']
        expected = {(tuple(map(Fraction, atom['clock'])),
                     tuple(tuple(map(Fraction, row)) for row in atom['terminal'])):
                    Fraction(atom['weight']) for atom in record['atoms']}
        got = compiler.compile_route(name, center, 5, base)
        assert got == expected, (name, base)
        checked += 1
        moving = [weight for (_, terminal), weight in got.items()
                  if any(coefficient for row in terminal for coefficient in row[:4])]
        expected_weight = ('11/48', '29/128', '29/144', '107/576', '233/1152')[base]
        assert moving == [Fraction(expected_weight)]
        print('[G4] base mod 5:', base, 'moving-terminal weight', moving[0], flush=True)
    assert checked == 5
    print('[G4] matched exact kernel tables:', checked)


def control_p8():
    sys.path.insert(0, str(INVESTIGATION / 'eleven_access' / 'restricted_theorem_check'))
    import check_controls
    from exact_model import Sheet

    split = Sheet.exits
    indices = []

    def checked_exits(self, guard, speed, drift):
        expected_phase = dot(guard, (*self.point, self.base, 1))
        assert (6 - self.base) % self.stride == 0
        expected_phase += (6 - self.base) // self.stride * dot(guard, (*self.slope, self.stride, 0))
        coefficients = tuple(dot(guard, col) for col in self.basis)
        divisor = gcd(speed, *coefficients)
        index = speed // divisor
        children = list(split(self, guard, speed, drift))
        supported = {r for r, i, c in children if (6 - c.base) % c.stride == 0 and i == index}
        assert supported == set(range(expected_phase % divisor, speed, divisor))
        assert len([1 for _, _, c in children if (6 - c.base) % c.stride == 0]) == index
        indices.append(index)
        print('[G8] guard, speed, coefficients, index, residues:', guard, speed, coefficients,
              index, tuple(sorted(supported)))
        yield from children

    Sheet.exits = checked_exits
    check_controls.class_six()
    print('[G8] selected path indices:', indices)


def control_isolated_wall():
    sys.path.insert(0, str(INVESTIGATION / 'cross_length' / 's17_share'))
    from wall_kernel import DRIFT_ROOT, drift, kernel, profile, root, step

    cases = ((15, {6: 5, 8: 16}, 9), (17, {8: 20, 12: 5, None: 2}, 13))
    for n, expected, expected_signature_count in cases:
        orbit, _, threshold, speed, _ = drift(n)
        section = range(threshold - speed, threshold)
        periods = {}
        signature_count = 0
        k = n // 2
        complement = tuple(i for i in range(k + 1) if i not in (0, 1, 2, 4, k))
        for c in section:
            pre, period, _, _ = kernel(n, c, orbit)
            periods[c] = period
            if period is None:
                continue
            state = profile(n, 0, DRIFT_ROOT, [c if i == k - 1 else 100000 + i
                                                for i in complement])
            cycle = []
            for t in range(pre + period):
                if t >= pre:
                    cycle.append(root(state))
                state, _, _ = step(state, t % len(orbit) == 0)
            signature_count += min(cycle) == (-1, 2, -1, 0) and period == 8
        counts = Counter(periods.values())
        assert counts == expected
        assert signature_count == expected_signature_count
        # The source proof gives a direction increasing c by one while
        # decreasing an independent high gap by one, at fixed anchor and base.
        admitted = sum(Fraction(v, speed) for p, v in counts.items() if p is not None)
        print('isolated-wall section:', n, 'section range:', (section.start, section.stop - 1),
              'index:', speed, 'admitted:', admitted, 'period weights:',
              {p: str(Fraction(v, speed)) for p, v in counts.items()},
              'period-eight signature weight:', Fraction(signature_count, speed),
              'unresolved:', tuple(c for c, period in periods.items() if period is None))


if __name__ == '__main__':
    control_p4()
    control_p8()
    control_isolated_wall()
