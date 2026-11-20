"""Trace the canonical continuation downstream of the period-66 ingress."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from pathlib import Path
import fcntl
import json
import tempfile

from stride_inventory import ROOT, tangent_period, verify


def run(ingress, continuation, root):
    recipient_histories = Counter()
    sample_histories = {}

    class Selected(verify.Sections):
        def arrival(self, poly, chart, sheet, ranks):
            self.stats['arrivals', ranks] += 1
            if tangent_period(sheet) != 66:
                return
            canonical, intercept, direction, axis = verify.canonical(sheet)
            assert canonical[1] == 6 and not any(intercept)
            verify.check_ingress_chart(chart, intercept, direction, axis)
            bounds = verify.ingress_interval(poly, chart, axis)
            list(self.add(sheet, bounds))
            self.stats['ingress', 'all'] += 1

    class SelectedQueue(verify.IngressQueue):
        def append(self, item):
            _, _, sheet, ranks, history, _ = item
            if sheet[2] == 1 and tangent_period(sheet) == 66:
                signature = ''.join(map(str, ranks))
                recipient_histories[signature] += 1
                sample_histories.setdefault(signature, history)
            super().append(item)

    original_sections = verify.Sections
    original_queue = verify.IngressQueue
    verify.Sections = Selected
    verify.IngressQueue = SelectedQueue
    try:
        with redirect_stdout(StringIO()):
            sections = verify.capture_ingress(ingress, root, argparse.Namespace(shorts=500, ingress_budget=50000))
    finally:
        verify.Sections = original_sections
        verify.IngressQueue = original_queue

    assert sections.stats['ingress', 'all'] == 3248
    initial_count = len(sections.coverage)
    decrements = Counter()
    readers_of_eleven = {}
    reference = verify.cover.reply

    def record_reply(reader, point, chart, sheet, max_short):
        answer = reference(reader, point, chart, sheet, max_short)
        kind, _, exits, guards = answer
        if kind != 'stage':
            return answer
        inside = (Q(0), Q(1))
        for row in sorted(guards):
            assert len(set(row[1:])) == 1
            inside = verify.clip(inside, (row[0] - row[1], row[1]))
            assert inside is not None
        distances = [verify.gap_row(exit[0]) for exit in exits]
        for (_, descriptor, _, children), (slope, offset) in zip(exits, distances):
            winning = verify.restrict(inside, [(a-slope, b-offset) for a, b in distances
                                               if (a, b) != (slope, offset)])
            if winning:
                decrement = descriptor[0]
                decrements[decrement, bool(children)] += 1
                if decrement % 11 == 0 and children:
                    _, _, _, axis = verify.canonical(sheet)
                    width = Q(-6 * chart[7 * axis], chart[7 * axis + 6])
                    # The synthetic gap parameter has width `width` along
                    # the first changing coordinate. Undo that change and
                    # the exit decrement to evaluate the guard on 6u.
                    tangent_increment = decrement * slope * sheet[13 + axis] / width
                    assert tangent_increment.denominator == 1
                    first_prime_free = decrement
                    while first_prime_free % 11 == 0:
                        first_prime_free //= 11
                    if int(tangent_increment) % (decrement // first_prime_free):
                        readers_of_eleven.setdefault(decrement, (sheet, descriptor, point, chart,
                                                                  int(tangent_increment)))
        return answer

    verify.cover.reply = record_reply
    try:
        with verify.reader_process(continuation) as reader, redirect_stdout(StringIO()):
            try:
                assert verify.walk(reader, sections, argparse.Namespace(shorts=500, line_budget=100000))
            finally:
                reader.stdin.close()
            assert reader.wait() == 0
    finally:
        verify.cover.reply = reference
    return dict(base_root=root, initial_sections=initial_count, final_sections=len(sections.coverage),
                period_66_history_counts=dict(recipient_histories), period_66_samples=sample_histories,
                decrement_counts=[(*key, n) for key, n in sorted(decrements.items())],
                nonterminal_readers_of_eleven=readers_of_eleven)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roots', type=int, nargs='+', choices=(0, 2, 4), default=(0, 2, 4))
    args = parser.parse_args()
    with (ROOT / '.task-monitor/compute.lock').open() as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix='task144-prime-') as temp:
            work = Path(temp)
            ingress = verify.build(work / 'ingress', False)
            continuation = verify.build(work / 'continuation', True)
            with ProcessPoolExecutor(max_workers=len(args.roots)) as pool:
                results = list(pool.map(run, [ingress]*len(args.roots), [continuation]*len(args.roots), args.roots))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
