"""Replay the frozen repeated-failure class with the paired checkpoint mark."""
from fractions import Fraction as Q
import json
import subprocess

import cohort1 as orbit
import quotient as box
import verify as route
from sections import prefixes


def midpoint_input(B,pin,kappa):
    a,b,c,h,g = (B*(low+high)/2 for low,high in orbit.BOX)
    z,r1 = 1,17
    p,r,N = orbit.prefix_coordinates(B,a,b,c,h,g,z,r1)
    n = round((h+3)/6)
    m = round((g-1+5*n+z)/2)
    N = round(N)
    r = 6*round((r-2)/6)+2
    residue = ((pin-kappa)//2-B+2*r)%36
    p = 36*round((p-residue)/36)+residue
    qR = tuple(int(row.eval(B,p,r)) for row in box.rows_at_switch(1,z,r1))
    state = orbit.translate(qR,-N,(-20,-4,-2,2,-4))
    state = orbit.translate(state,-m,(-10,-2,0,2,-2))
    state = orbit.translate(state,n,orbit.V)
    source = orbit.translate(state,-1,(-5,0,1,2,-1))
    a,b,c,d,e = source
    chart = tuple(Q(u,B) for u in (a,b,c,2*(b-e)-B,B-c-d))
    assert all(lo<=x<=hi for x,(lo,hi) in zip(chart,orbit.BOX))
    assert (b-e)%3 == 0 and (2*B+2*p-4*r+kappa-pin)%72 == 0
    h1 = 2*(state[1]-state[4])-B+1
    assert h1 == 6*n-2
    after_first = orbit.translate(state,-n,orbit.V)
    g1 = B-after_first[2]-after_first[3]-1
    assert g1 == 2*m-z
    after_prep = orbit.translate(after_first,m,(-10,-2,0,2,-2))
    assert after_prep[0]+after_prep[2]-3*(B-1)//2 == 22*N+r1
    return source, chart, N


def main():
    kappa = next(rows[-1] for rows in prefixes() if rows[0] == (1,1,17))
    observed = []
    for phase,f_expected in zip((5,11,17,23),((0,1),(1,0),(1,0),(0,1))):
        base = 90_000_000+phase
        for pin,answer in zip((603,639),f_expected):
            args,chart,N = midpoint_input(base,pin,kappa)
            observed.append((base,pin,answer,args,chart,N))
    subprocess.run(['cc','-O2','word2.c','-o','word-clock'],check=True)
    data = ''.join(' '.join(map(str,(B,*args)))+'\n' for B,_,_,args,_,_ in observed)
    clocks = subprocess.check_output(['./word-clock'],input=data,text=True).splitlines()
    assert len(clocks) == len(observed)
    output = []
    for (B,pin,answer,source,chart,N),row in zip(observed,clocks):
        t,period = map(int,row.split())
        f = min(1,max(0,500*(Q(t,B)-route.theta(chart)+Q(1,1000))))
        assert f == answer, (B,pin,t,f)
        output.append({'base': B, 'r2': pin, 'source': source, 'N': N,
                       'word_clock': t, 'period': period, 'F': int(f)})
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    main()
