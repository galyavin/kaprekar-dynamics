"""Reproduce the marked box projection and selected canonical normal strides."""
from ast import literal_eval
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
sys.path.insert(0, str(ACCESS / 'z1b_rank_one'))
import verify


def box_table():
    directory = ACCESS / 'r2_routes'
    fresh = subprocess.check_output([sys.executable, '-B', str(directory / 'control_weights.py')])
    assert fresh == (directory / 'control-weights.json').read_bytes()
    data = json.loads(fresh)
    table = {}
    pools = {}
    for base, labels in data['marked_profiles'].items():
        base = int(base)
        pool = {}
        names = []
        for rows in labels.values():
            names.append({'-11/324': 'S', '29/54': 'M2', '139/144': 'L'}[rows[0]['coefficients'][0]])
            for row in rows:
                key = tuple(row['coefficients'])
                pool[key] = pool.get(key, Q(0)) + Q(row['weight']) / 2
        assert sum(pool.values()) == 1
        assert pools.setdefault(base % 12, pool) == pool
        assert table.setdefault(base % 24, names) == names
    return dict(marked_profiles=table, box_original_mass=data['box_original_mass'],
                label_original_mass=data['each_marked_label_original_mass'])


def encode_line(base, period, slope, point, tangent):
    return (base, period, 1, *slope, *point, *tangent)


def original_rank_one_sections():
    witness = ACCESS / 'z1_stages/third-keep-route.out'
    for line in witness.read_text().splitlines():
        if not line.startswith('stage') or line.endswith('None'):
            continue
        sheet = literal_eval(line[line.rfind('('):])
        if sheet[2] == 1:
            yield sheet


def direct_normal_period(sheet):
    _, period, _, *rows = sheet
    common = gcd(*rows[10:])
    tangent = [v // common for v in rows[10:]]
    # A shift by 6m is feasible exactly when an integer vector w solves
    # period*w - m*slope = c*tangent for some integer c.  Since the
    # tangent is primitive, it suffices to search c modulo period.
    slope = rows[:5]
    for m in range(1, period // 6 + 1):
        if any(all((m * v + c * u) % period == 0 for u, v in zip(tangent, slope))
               for c in range(period)):
            return 6 * m
    raise AssertionError('the original period must be feasible')


def normal_table():
    sections = [('archived', s) for s in original_rank_one_sections()]
    # Known physical lifts with nonzero normal slope also exercise the
    # denominator test. Every base and point below lies in root zero.
    for period, slope in ((12, (6, 0, 0, 0, 0)), (18, (6, 0, 0, 0, 0)),
                          (72, (0, 12, 18, 0, 0))):
        tangent = (6, 12, 0, -6, 0)
        sections.append(('constructed', encode_line(0, period, slope, (0,) * 5, tangent)))
    output = []
    for source, sheet in sections:
        section, _, _, _ = verify.canonical(sheet)
        period = direct_normal_period(sheet)
        assert section[1] == period
        output.append(dict(source=source, sheet=sheet, canonical=section, direct_period=period))
    assert len(output) == 5
    assert [row['direct_period'] for row in output] == [6, 6, 12, 18, 72]
    return output


if __name__ == '__main__':
    print(json.dumps(dict(box=box_table(), normal=normal_table()), indent=2))
