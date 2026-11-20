"""Read every inherited arrival with the separate ordinary-difference reconstruction.

The Task 145 reader finds the physical words, clips source cells and transports
the ordinary q lattice. This audit computes centered coset orders from that
output and tests all exceptional saturated intervals directly.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
import fcntl
import json
from math import gcd, lcm
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

BALANCE = (5, 4, 3, 2, 1)


def private_root(root):
    import cover_exact as model

    class CountingReader(model.Reader):
        def __init__(self, executable):
            super().__init__(executable)
            self.rank_strides = set()

        def __call__(self, chart, point, sheet, budget=None):
            if len(sheet.basis) >= 2:
                self.rank_strides.add((len(sheet.basis), sheet.stride))
                assert sheet.stride in (6, 12)
            return super().__call__(chart, point, sheet, budget)

    class Ingress:
        def __init__(self, reader):
            self.reader = reader
            self.arrivals = Counter()
            self.exceptional = model.Lines(reader, root)
            self.decrements = Counter()

        def add_sheet(self, poly, chart, sheet):
            tangent, = sheet.basis
            step = tuple(6*x for x in tangent)
            translation = tuple(6*x-sheet.stride*v for x, v in zip(sheet.slope, BALANCE))
            period = sheet.stride * lcm(*(
                abs(t)//gcd(t, s) for t, s in zip(step, translation) if t))
            # Check homogeneity of the exact joint sheet and restoration, including
            # coordinates with no tangent variation.
            i = next(i for i, t in enumerate(step) if t)
            assert all(t*translation[i] == s*step[i] for t, s in zip(step, translation))
            assert (period // sheet.stride * translation[i]) % step[i] == 0
            self.arrivals[sheet.stride, period] += 1
            if period == 66:
                assert root % 2 == 0 and sheet.stride == 6
                self.exceptional.add_sheet(poly, chart, sheet)
            else:
                assert period in (6, 12)

        def run(self):
            # Every insertion on the same normal sheet uses the same affine chart.
            # Read its final interval union once. These cells cover all eleven
            # inherited tangent phases simultaneously, on the common physical line.
            charts = {key: chart for key, chart, _ in self.exceptional.work}
            for key, covered in self.exceptional.cover.items():
                sheet, _ = key
                chart = charts[key]
                for lo, hi in covered:
                    self.read_interval(sheet, chart, lo, hi)

        def read_interval(self, sheet, chart, low, high):
            work = [(low, high)]
            while work:
                lo, hi = work.pop()
                for n in range(2, 35):
                    point = (lo + (hi-lo)/n, 0, 0, 0, 0)
                    kind, read, guards = self.reader(chart, point, sheet, budget=250)
                    if kind != 'wall':
                        break
                assert kind == 'stage'
                for g in guards:
                    accepted = model.interval_test(lo, hi, g)
                    assert accepted is not None
                    if accepted != (lo, hi):
                        rejected = model.interval_test(lo, hi, tuple(-a for a in g))
                        assert rejected is not None
                        work.append(rejected)
                    lo, hi = accepted
                _, _, unpaired, prechart, _, drift, races = read
                assert unpaired == 0
                winners = 0
                for time, (guard, s) in races.items():
                    interval = lo, hi
                    for rival in races:
                        if rival != time:
                            interval = model.interval_test(*interval, model.combine(1, rival, -1, time))
                            if interval is None:
                                break
                    if interval is None:
                        continue
                    winners += 1
                    assert sum(ai*vi for ai, vi in zip(guard, BALANCE)) + 6*guard[5] == 0
                    assert model.central(model.exit_chart(prechart, drift, time))
                    self.decrements[6*s] += 1
                assert winners

    reader = CountingReader('./reader')
    audit = model.Ranks(reader, root, 2)
    inventory = Ingress(reader)
    audit.lines = inventory
    try:
        with redirect_stdout(sys.stderr):
            audit.run()
    finally:
        reader.close()
    return {
        'root': root, 'rank_strides': sorted(reader.rank_strides),
        'arrivals': sorted([*key, count] for key, count in inventory.arrivals.items()),
        'exceptional_normal_lines': len({sheet for sheet, _ in inventory.exceptional.cover}),
        'exceptional_sections': len(inventory.exceptional.cover),
        'terminal_decrements': sorted(inventory.decrements.items()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--private', action='store_true')
    args = parser.parse_args()
    if args.private:
        with ProcessPoolExecutor(max_workers=6) as pool:
            print(json.dumps(list(pool.map(private_root, range(6))), indent=2))
        return
    root = Path(__file__).resolve().parents[4]
    model = root / 'kaprekar-investigation/eleven_access/restricted_theorem_check'
    with tempfile.TemporaryDirectory(prefix='astra-148-period-') as tmp:
        folder = Path(tmp)
        for name in ('reader.cpp', 'cover_exact.py', 'exact_model.py'):
            shutil.copyfile(model/name, folder/name)
        shutil.copyfile(__file__, folder/'period_audit.py')
        with (root/'.task-monitor/compute.lock').open('rb') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            print('Exclusive lock acquired, running six roots in private snapshot', file=sys.stderr, flush=True)
            subprocess.run(['clang++', '-std=c++17', '-O3', '-ftrapv', 'reader.cpp', '-o', 'reader'],
                           cwd=folder, check=True)
            subprocess.run([sys.executable, '-B', 'period_audit.py', '--private'], cwd=folder, check=True)


if __name__ == '__main__':
    main()
