#!/usr/bin/env python3
"""Regenerate the even exclusion tree and record every leaf.

Imports exclude.py read-only from a private snapshot.  explore() below is a
line-for-line copy of exclude.explore with one addition: each terminal class
(coset point + basis, preterminal time form, constant preterminal state,
route history) is appended to DUMP.  Standard output reproduces exclude.py's
__main__ output so its hash can be compared with exclusion.txt.

Usage: dump_tree.py SNAPSHOT_DIR DUMP.jsonl.gz > exclusion_rerun.txt
"""
import gzip
import json
import sys

sys.path.insert(0, sys.argv[1])
import exclude as ex                                  # noqa: E402
from affine_routes import VARIABLES, ZERO, normalized_interval, repeat   # noqa: E402
from primitive import step                            # noqa: E402

DUMP = gzip.open(sys.argv[2], "wt")


def explore(state, elapsed, x, basis, depth=0, history=()):
    ex.metrics['entry', depth] += 1
    for gap in range(ex.MAX_SHORT):
        word, delta, modulus, num = ex.macro(state)
        if modulus:
            break
        state, elapsed = step(state, word[0])[0], elapsed + 1
    else:
        raise ArithmeticError(f'no substantial return at depth {depth}')
    assert depth < ex.MAX_BLOCKS, history
    speed = (len(word), delta)
    ex.metrics['blocks', depth, speed] += 1
    ex.families[(depth, speed)].add((num.row, modulus))
    anchor, radial = state[:2]
    if (not any(state[2].row) and anchor.row == state[3].row
            and radial.row == state[4].row):
        fiber = ex.plain((state[2], state[3] - anchor, state[4] - radial))
        ex.sheet_fibers[speed].add(fiber)
    for remainder, point, next_basis in ex.remainders(num, modulus, x, basis):
        assert ex.positive((num - remainder) / modulus)
        next_state, time = repeat(state, elapsed, (word, delta, modulus, num), remainder)
        if all(not any(s.row) for s in next_state):
            assert ex.positive((num - remainder) / modulus + 1 - ex.PRETERMINAL_RETURNS)
            ex.record_terminal(next_state, time, word, delta)
            ex.routes.add(history + (speed,))
            endpoint = ex.plain(next_state)
            spot = tuple(e - ex.PRETERMINAL_RETURNS * d for e, d in zip(endpoint, delta))
            cut = time - ex.PRETERMINAL_RETURNS * len(word)
            DUMP.write(json.dumps(dict(point=list(point), basis=[list(b) for b in next_basis],
                                       cut=[str(c) for c in cut.row] + [str(cut.constant)],
                                       endpoint=list(endpoint), spot=list(spot),
                                       hist=[(h[0], list(h[1])) for h in history + (speed,)],
                                       word_len=len(word), delta=list(delta))) + "\n")
        else:
            explore(next_state, time, point, next_basis, depth + 1, history + (speed,))


if __name__ == '__main__':
    basis = tuple(tuple(6 * int(i == j) for j in range(5)) for i in range(5))
    explore(VARIABLES, ZERO, (0,) * 5, basis)
    DUMP.close()
    for metric, value in sorted(ex.metrics.items(), key=lambda p: repr(p[0])):
        print(metric, value)
    print('uniform base threshold', ex.base_threshold)
    for row, offsets in sorted(ex.bounds.items()):
        lo, hi = normalized_interval(row)
        print('leading time', row, 'range', (str(lo), str(hi)), 'constants', offsets)
    for key, values in sorted(ex.sheet_fibers.items()):
        print('sheet', key, 'fibers', sorted(values))
    for key, values in sorted(ex.terminal_fibers.items()):
        print('terminal', key, 'endpoints', sorted(values))
    for route in sorted(ex.routes):
        print('route', route)
    print(ex.macro.cache_info())
