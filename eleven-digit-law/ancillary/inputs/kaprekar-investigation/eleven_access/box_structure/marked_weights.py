"""Read the Task 071 affine leaves and check the marked route test."""
from collections import defaultdict
from fractions import Fraction as Q
import json
from math import gcd
from pathlib import Path


def determinant(a,b,c):
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))


def main():
    root = Path(__file__).resolve().parents[3]
    source = root / 'docs/drafts/research-direction/replay-task071/results/leaves.jsonl'
    weights = defaultdict(lambda: defaultdict(Q))
    with source.open() as lines:
        for line in lines:
            leaf = json.loads(line)
            if leaf['label'][:3] != [0,1,17] or leaf['label'][3] not in (603,639) or leaf['rho'] != 2:
                continue
            columns = leaf['columns']
            step = gcd(*(v[0] for v in columns))
            weight = Q(216*step, abs(determinant(*columns)))
            for base in range(5,72,6):
                if (base-leaf['point'][0]) % step == 0:
                    weights[base,leaf['label'][3]][leaf['lane'],tuple(leaf['profile'])] += weight
    out = []
    for (base,pin),table in sorted(weights.items()):
        assert sum(table.values()) == 1
        long_mass = sum(weight for (lane,_),weight in table.items() if lane == 'L')
        assert long_mass == (int(base % 24 in (11,17)) if pin == 603 else int(base % 24 in (5,23)))
        out.append({'base_mod_72':base,'r2':pin,'plateau':long_mass,
                    'profiles':[{'lane':lane,'coefficients':coeffs,'conditional_weight':str(weight)}
                                for (lane,coeffs),weight in sorted(table.items())]})
    assert len(out) == 24
    print(json.dumps(out,indent=2,default=int))


if __name__ == '__main__':
    main()
