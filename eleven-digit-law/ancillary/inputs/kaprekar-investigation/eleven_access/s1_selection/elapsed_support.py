"""Separate first elapsed-clock differences from first geometric differences."""
from collections import Counter
from itertools import combinations, zip_longest
import json
from pathlib import Path
import sys

from summarize import estimate, explore_parse, interval, read


def main():
    rows, hashes = read(Path(sys.argv[1]),'explore',4000,explore_parse)
    result = {}
    for r,s in combinations(range(24),2):
        first = Counter()
        per_stage = [[] for _ in range(15)]
        for entry,arms in rows:
            differences = [0.]*15
            if entry>=0 and all(not a[0] for a in arms):
                left,right = arms[r][2],arms[s][2]
                for k,(a,b) in enumerate(zip_longest(left,right)):
                    differences[k] = (a[0] if a else 0.)-(b[0] if b else 0.)
                for k,d in enumerate(differences):
                    if d:
                        first[(k+1,left[k][1][0] if k<len(left) else -1,
                               right[k][1][0] if k<len(right) else -1)] += 1
                        break
                else:
                    first['none'] += 1
            else:
                first['excluded'] += 1
            for col,x in zip(per_stage,differences):
                col.append(x)
        result[f'{r}-{s}'] = dict(first={str(k):v for k,v in first.items()},
                                  per_stage={i+1:interval(estimate(x)) for i,x in enumerate(per_stage)
                                             if any(x)})
    print(json.dumps(dict(shard_hashes=hashes,rows=len(rows),paired_support=result),indent=2))


if __name__=='__main__':
    main()
