"""Read the inherited base periods before rank one and the full line torsion."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from math import gcd, lcm
from pathlib import Path
import fcntl
import json
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = ACCESS.parents[1]
sys.path.insert(0, str(ACCESS / 'z1b_rank_one'))
import verify


def tangent_period(sheet):
    _, period, rank, *rows = sheet
    assert rank == 1
    stride = gcd(*rows[10:])
    primitive = tuple(v // stride for v in rows[10:])
    k = next(i for i, v in enumerate(primitive) if v)
    slope = Q(rows[k], primitive[k])
    assert slope.denominator == 1
    assert tuple(slope * v for v in primitive) == tuple(rows[:5])
    return period * stride // gcd(stride, int(slope))


def run(exe, root):
    rank_periods = Counter()
    periods = Counter()
    strata = Counter()
    reference = verify.cover.reply

    def record_reply(reader, point, chart, sheet, max_short):
        rank_periods[sheet[2], sheet[1]] += 1
        return reference(reader, point, chart, sheet, max_short)

    class Inventory(verify.Sections):
        def arrival(self, poly, chart, sheet, ranks):
            super().arrival(poly, chart, sheet, ranks)
            period = tangent_period(sheet)
            periods[sheet[1], period] += 1
            strata[sheet[0] % period, period] += 1

    original_sections = verify.Sections
    verify.Sections = Inventory
    verify.cover.reply = record_reply
    try:
        with redirect_stdout(StringIO()):
            sections = verify.capture_ingress(exe, root, argparse.Namespace(shorts=500, ingress_budget=50000))
    finally:
        verify.cover.reply = reference
        verify.Sections = original_sections
    assert sum(periods.values()) == sections.stats['ingress', 'all']
    return dict(base_root=root, source_reads=[(*k, v) for k, v in sorted(rank_periods.items())],
                rank_one_periods=[(*k, v) for k, v in sorted(periods.items())],
                rank_one_strata=[(*k, v) for k, v in sorted(strata.items())],
                common_modulus=lcm(12, *(k[1] for k in periods)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roots', type=int, nargs='+', choices=range(6), default=list(range(6)))
    parser.add_argument('--jobs', type=int, default=6)
    args = parser.parse_args()
    with (ROOT / '.task-monitor/compute.lock').open() as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix='task144-') as tmp:
            executable = verify.build(Path(tmp) / 'charts', False)
            with ProcessPoolExecutor(max_workers=min(args.jobs, len(args.roots))) as pool:
                rows = list(pool.map(run, [executable] * len(args.roots), args.roots))
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
