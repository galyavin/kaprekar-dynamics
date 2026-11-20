"""Find the first incompatible S^5 phase on the archived enlarged-box ray."""
from fractions import Fraction as Q
import json
import subprocess

from cohort1 import digit_rows, indices

S = '5061728394a'
W = (0,-12,18,-6,-18)
V = (5,4,3,2,1)
START = (60005, (56516,44321,34724,25107,14328))


def walk(base, state, steps):
    ix = indices(S)
    for phase in range(steps):
        digits = digit_rows(base, state)
        bad = [(i, j, digits[i]-digits[j]) for i,j in zip(ix,ix[1:]) if digits[i]<digits[j]]
        if bad:
            return dict(phase=phase, state=state, violated=bad,
                        actual_order=''.join(format(i,'x') for i in sorted(range(11),key=lambda i:(-digits[i],i))))
        state = tuple(digits[ix[i]]-digits[ix[-1-i]] for i in range(5))
    raise ValueError('No failed guard found')


def main():
    subprocess.run(['cc','-O2','word2.c','-o','word-clock'],check=True)
    cases = []
    B, q = START
    for scale in (1, 1201):
        base = B*scale
        q0 = tuple(scale*x for x in q)
        q1 = tuple(a+b for a,b in zip(q0,W))
        for b, x in ((base,q0),(base+12,q1)):
            h = 2*(x[1]-x[4])-b+1
            cohort = (x[1]-x[4]) % 3
            n = (h-(-4,-2,-6)[cohort])//6
            cases.append(dict(base=b, source=x, H=h, n=n, first_bad=walk(b,x,5)))
    data = ''.join(' '.join(map(str,(r['base'],*r['source'])))+'\n' for r in cases)
    clocks = subprocess.check_output(['./word-clock'],input=data,text=True)
    for result, row in zip(cases,clocks.splitlines()):
        t,p = map(int,row.split())
        result['first_cyclic_word'],result['period'] = t,p
    assert [x['first_cyclic_word'] for x in cases[:2]] == [6980,6974]
    assert Q(2*(q[1]-q[4])-B,B) == -Q(19,B)
    print(json.dumps(cases,indent=2))


if __name__ == '__main__':
    main()
