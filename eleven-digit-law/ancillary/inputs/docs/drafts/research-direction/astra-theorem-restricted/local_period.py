"""Evaluate the original-position period witnesses at interior rational points."""
from fractions import Fraction as Q
from pathlib import Path
import json

ACCESS = Path(__file__).resolve().parents[4] / 'kaprekar-investigation/eleven_access'


def cyclic_chamber(q):
    start = q
    for _ in range(5):
        digits = (1, *(v for pair in zip(q, (1-u for u in reversed(q))) for v in pair))
        assert all(a>b for a,b in zip(digits,digits[1:]))
        q = tuple(digits[i]-digits[-1-i] for i in range(5))
    assert q == start


def mean_clipped(readers):
    assert sum(w for w,t in readers) == 1
    assert all(0<t<1 for w,t in readers)
    return sum(w*t for w,t in readers)


def cube():
    table = json.loads((ACCESS/'stage_graph/cell-lotteries.json').read_text())
    y = tuple(Q(n,table['center_grid']) for n in table['center_q'])
    cyclic_chamber(y)
    x = tuple(6*a-(5-i) for i,a in enumerate(y))
    means = {}
    for residue in (3, 5, 9):
        rows = table['tables'][table['mod12_to_table'][str(residue)]]['clock']
        costs = [(Q(row['weight']), sum(Q(a,row['denominator'])*b
                  for a,b in zip(row['numerator'],(*x,1)))) for row in rows]
        means[residue] = mean_clipped(costs)
    assert means[9] == means[3] > means[5]
    print('Cube centre, 9 minus 5:', means[9]-means[5])


def box():
    rows = json.loads((ACCESS/'box_structure/limit_weights.json').read_text())['graphs']
    a,b,c,h,g = (Q(left+right,200000) for left,right in (
        (94040,94130),(73770,73850),(57940,58015),(20,80),(200,300)))
    y = a,b,c,1-c-g,b-(1+h)/2
    cyclic_chamber(y)
    n = (a+c-Q(3,2)-Q(11,2)*h-5*g)/22
    p,r = c-h/2-2*n,b-Q(3,2)*h-g-4*n
    means = {}
    for residue in (5,11):
        costs = []
        for row in rows:
            A,P,R = map(Q,row['s'])
            costs.append((Q(row[f'arm_{residue}_units'])/28512, a+A+P*p+R*r))
        means[residue] = mean_clipped(costs)
    assert means[5] > means[11]
    print('Box midpoint, 5 minus 11:', means[5]-means[11])


if __name__ == '__main__':
    cube()
    box()
