"""Read the full borrowed digits with two formal parameters on the marked 66-period line."""
from fractions import Fraction as Q
from pathlib import Path
import ast
import re

ROOT = Path(__file__).resolve().parents[4]
ACCESS = ROOT / 'kaprekar-investigation/eleven_access'
V = (5, 4, 3, 2, 1)


def minus(a, b):
    return tuple(x-y for x,y in zip(a,b))


def raw_step(state, interval):
    # Rows denote a*B + b*h + c in six times physical digit units.
    digits = [(6, 0, -6)]
    for i, (base, slope, offset) in enumerate(state):
        digits.append((V[i]+base, slope, offset - 6*(i == 4)))
        base, slope, offset = state[4-i]
        digits.append((V[i]-base, -slope, -offset - 6*(i < 4)))
    lo, hi = interval
    mid = (lo+hi)/2
    ordered = sorted(digits, key=lambda row:(row[0]+mid*row[1],row[2]), reverse=True)
    guards = [minus(a,b) for a,b in zip(ordered, ordered[1:])]
    height = tuple(x+y for x,y in zip(state[-1],(1,0,-6)))
    guards += [height]
    for a,b,c in guards:
        assert a+b*lo >= 0 and a+b*hi >= 0
        assert a+b*mid > 0 or a == b == 0 and c >= 0
    return tuple(minus(minus(a,b),(v,0,0)) for a,b,v in zip(ordered, reversed(ordered), V)), guards


def verify():
    line = (ACCESS/'class_compat/period66-witness.out').read_text()
    vertices, rest = line.split(' vertices ')[1].split(' chart ')
    chart, rest = rest.split(' sheet ')
    sheet = ast.literal_eval(rest.split(' history ')[0])
    chart = ast.literal_eval(chart)
    vertices = ast.literal_eval(re.sub(r'Fraction\(([-0-9]+), (\d+)\)',r'"\1/\2"',vertices))
    base, stride, rank, *data = sheet
    slope, point, tangent = data[:5], data[5:10], data[10:]
    assert (base, stride, rank) == (0, 6, 1)
    assert all(11*s == t for s,t in zip(slope,tangent))
    u = tuple(t//66 for t in tangent)
    assert tangent == [66*x for x in u] and u[0] == 1
    parameters = []
    for raw in vertices:
        gaps = [Q(a) for a in raw]
        source = [6*sum(gaps[i:])-v for i,v in zip(range(1,6),V)] + [1]
        x = [sum(a*b for a,b in zip(chart[i:i+6],source))/chart[i+6] for i in range(0,35,7)]
        assert x == [x[0]*z for z in u]
        parameters.append(x[0])
    limits = min(parameters), max(parameters)
    assert limits == (0, Q(3,14))
    assert {6*k % 66 for k in range(11)} == {6*j for j in range(11)}
    # The inherited phases are h = 6*k mod 66. Work on the larger real h line.
    # Search for the same symbolic section with h replaced by h-m, m>0.
    state = tuple((0,s,p) for s,p in zip(u,point))
    states = [state]
    guards = []
    for end in range(1,400):
        state, rows = raw_step(state,limits)
        states.append(state)
        guards.append(rows)
        for start, old in enumerate(states[:-1]):
            if [v[:2] for v in old] != [v[:2] for v in state]:
                continue
            i = next(i for i,f in enumerate(old) if f[1])
            drop = Q(old[i][2]-state[i][2], old[i][1])
            if drop <= 0 or any(a[2]-b[2] != drop*a[1] for a,b in zip(old,state)):
                continue
            assert drop.denominator == 1
            falling = [g for step in guards[start:] for g in step if g[1]>0]
            # The maximum cutoff H(B) for h among declining guards is bounded:
            # all a>=0 at h/B=0, at least one a=0 has positive slope.
            assert falling and all(a>=0 for a,b,c in falling)
            cutoff = max(Q(-c,b) for a,b,c in falling if a==0)
            assert all(f[0]==0 for f in old)
            print('Uncompressed 66-period witness on 0 < h/B < 3/14')
            print('Subject to h = 6k mod 66 at B = 6k, each of the eleven tangent phases applies.')
            print('Short updates:',start,'word length:',end-start,'parameter decrement:',drop)
            print('Common leading clock slope in h/B:', Q(end-start,drop))
            print('Largest terminal lower guard threshold in h:',cutoff)
            return
    raise AssertionError('No translating return closed on the full interval')


if __name__ == '__main__':
    verify()
