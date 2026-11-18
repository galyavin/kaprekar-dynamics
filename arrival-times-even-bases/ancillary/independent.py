"""Independent evaluations of (5), (6), and the finite functional graphs."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, inf
import subprocess


def subtract(digits, base):
    """Return the four padded digits of a literal positional subtraction."""
    digits = sorted(digits)
    lo = hi = 0
    for a, b in zip(digits, reversed(digits)):
        lo, hi = lo*base+a, hi*base+b
    value = hi-lo
    out = []
    for _ in digits:
        value, digit = divmod(value, base)
        out.append(digit)
    return tuple(reversed(out))


def cycles(m):
    remaining = {(u,v) for v,u in combinations(range(1,m,2),2)}
    parts = []
    while remaining:
        p = min(remaining)
        cycle = set()
        while p in remaining:
            cycle.add(p)
            remaining.remove(p)
            a,b,c,d = sorted(subtract((*p,0,0),m))
            p = d-a,c-b
        assert p == min(cycle)
        parts.append(cycle)
    return parts


def brute(m, n, target):
    s, base = 2**n, m*2**n
    goal = {subtract((s*u,s*v,0,0),base) for u,v in target}
    distance = dict.fromkeys(goal,0)
    hist = Counter()
    for start in product(range(base), repeat=4):
        if len(set(start)) == 1:
            continue
        path = []
        p = start
        while p not in distance:
            path.append(p)
            distance[p] = inf
            p = subtract(p,base)
        t = distance[p]
        for point in reversed(path):
            t += 1
            distance[point] = t
        if distance[start] < inf:
            hist[distance[start]] += 1
    return dict(sorted(hist.items()))


class Coefficients:
    def __init__(self,m,n,target):
        self.m = m
        self.h = (m+1)//2
        self.residues = set()
        for u,v in target:
            for a,b in ((u,v),(v,u)):
                for i,j in product((-1,1),repeat=2):
                    self.residues.add((i*a % m,j*b % m))
        def is_seed(a,b):
            return (a%m,b%m) in self.residues
        self.seeds = {
            'Z': {r for r in range(m) if is_seed(r-1,-r)},
            'A': {r for r in range(m) if is_seed(2*r+1,2*r-1)}
        }
        self.vertices = {'S': ('Z',2**n), 'E': ('Z',2**(n+1)),
                         'U': ('A',self.h*2**n)}

    def coefficient(self, kind, u):
        seen = set()
        steps = 1
        while (kind,u%self.m) not in seen:
            r = u % self.m
            if r in self.seeds[kind]:
                return steps
            seen.add((kind,r))
            if kind == 'Z' and r in (1,self.h):
                name = 'S' if r == 1 else 'U'
            elif kind == 'A' and r in (self.h-1,self.h):
                name = 'E'
            else:
                break
            kind,u = self.vertices[name]
            steps += 1
        return inf

    def k(self,name):
        return self.coefficient(*self.vertices[name])


def progression(lower, upper, m, residue, valuation, offset):
    """u congruent to r mod m with exact valuation of u-offset."""
    half = 1 << valuation
    step = 2*half*m
    a = offset+half
    a += 2*half * (((residue-a)*pow(2*half,-1,m)) % m)
    first = a + ((lower-a+step-1)//step)*step
    return first, step, max(0,(upper-first)//step+1)


def quadratic_sum(fun, a, step, count):
    # Newton's finite-difference basis sums to successive binomial coefficients.
    y0,y1,y2 = (fun(a+i*step) for i in range(3))
    return count*y0 + comb(count,2)*(y1-y0) + comb(count,3)*(y2-2*y1+y0)


def root_weight(kind,base,u):
    if kind == 'Z':
        return (base-u)*(12*u-4)
    if kind == 'D':
        return 6*(base-u)
    return 24*(base-u)*(2*u-base)


def formula(m,n,target):
    """Evaluate (5) by residue and valuation, including the singleton roots."""
    s, base = 2**n, m*2**n
    data = Coefficients(m,n,target)
    rset = data.residues
    g4 = len(rset)*m*m
    g2 = 6*sum(m*abs(a-b)-(a-b)**2 for a,b in rset)-g4
    W = sum(24*(m-u)*(u-v) for u,v in target)
    hist = Counter({0:len(target), 1:s*s*W-len(target), 2:s*s*(g4+g2-W)})
    for t in range(2,n+2):
        hist[t+1] += 3*g4*4**(n+t-2)

    def clock(k):
        return k*(n+1)+1

    def arrival(kind,u,depth):
        r = u % m
        family = 'Z' if kind == 'Z' else 'A'
        if r in data.seeds[family]:
            return n+2
        if family == 'Z':
            if u == 1:
                return n+3 if (m-1) in data.seeds['Z'] else inf
            if r == data.h:
                return n+1+clock(data.k('U'))
            if r in (0,1):
                w = u-r
                t = n-((w & -w).bit_length()-1)
                return t+1+clock(data.coefficient('Z',2**t))
        else:
            if u == base//2:
                return n+3 if 1 in data.seeds['A'] else inf
            if r in (data.h-1,data.h):
                return n+1+clock(data.k('E'))
            if r == 0:
                t = n-depth-2
                return t+1+clock(data.coefficient('A',2**t))
        return inf

    def add(kind,a,step,count,depth):
        H = arrival(kind,a,depth)
        if H == inf:
            return
        hist[1+H] += quadratic_sum(lambda u:root_weight(kind,base,u),a,step,count)
        if kind == 'D':
            return
        for t in range(1,depth+1):
            def layer(u):
                if kind == 'Z':
                    return base**2*(4**t-4)+12*base*u-6*u*u
                return 2*base**2*(4**t-10)+48*base*u-24*u*u
            hist[1+H+t] += quadratic_sum(layer,a,step,count)

    add('Z',1,1,1,0)
    for kind in ('Z','D','A'):
        lower = {'Z':2, 'D':1, 'A':base//2+1}[kind]
        for r in range(m):
            # Split Z(1 mod m) by parity: for odd u the return uses v2(u-1).
            split = kind == 'Z' and r == 1
            for offset in ((0,1) if split else (0,)):
                for v in range(int(split),base.bit_length()):
                    a,step,count = progression(lower,base-1,m,r,v,offset)
                    if count:
                        add(kind,a,step,count,0 if offset else min(v,n))
    return dict(sorted((t,w) for t,w in hist.items() if w))


def limit5(r):
    m = 5
    c = Coefficients(m,r,{(3,1)})
    atoms = Counter()
    atoms[1] += Q(len(c.residues),25)
    atoms[2] += Q(2*sum(map(len,c.seeds.values())),25)
    for name,mass in (('U',2),('E',4),('S',2)):
        atoms[c.k(name)+2] += Q(mass,25)
    for t in range(1,5):
        atoms[c.coefficient('Z',2**t)+1] += Q(2*2**(4-t),3*25*15)
    for t in range(4):
        atoms[c.coefficient('A',2**t)+1] += Q(2**(4-t),6*25*15)
    values = [atoms[j] for j in range(1,6)]
    basin = sum(values)
    return values,basin,[sum(j**p*v for j,v in enumerate(values,1))/basin for p in (1,2)]


def graph(executable,m,n):
    result = subprocess.check_output([executable,str(m),str(n)],text=True)
    maxima, hist, points = {},{},{}
    for line in result.splitlines():
        label,*fields = line.split()
        a,b,*rest = map(int,fields)
        if label == 'cycle':
            maxima[a] = b
        elif label == 'weight':
            hist[a] = b
        else:
            points[label,a] = b,rest[0]
    assert sum(hist.values()) == (m<<n)**4-(m<<n)
    return maxima,hist,points


def tv(hist,m,n):
    total = sum(hist.values())
    model = {t:Q(3,4**(n+3-t)) for t in range(2,n+3)}
    model[1] = Q(1,4**(n+1))
    return sum(abs(Q(hist.get(t,0),total)-model.get(t,0)) for t in hist.keys() | model.keys())/2


def run(executable):
    cases = []
    for m,n,which in [(5,1,'all'),(5,2,'all'),(7,1,'empty'),(9,1,'first'),
                      (9,1,'all'),(13,1,'first two')]:
        parts = cycles(m)
        target = {'all':set.union(*parts),'empty':set(),'first':parts[0],
                  'first two':set.union(*parts[:2])}[which]
        if which == 'first two':
            assert len(parts)>2
        predicted, observed = formula(m,n,target),brute(m,n,target)
        assert predicted == observed,(m,n,which)
        cases.append(dict(m=m,n=n,target=which,formula=predicted,subtraction=observed))

    times, distances = [],[]
    for m in [*range(5,42,2),53]:
        parts = cycles(m)
        for n in range(7):
            check_clocks = m <= 41 and n >= 1
            check_tv = m in (27,29,31,53) and n <= 4
            if not (check_clocks or check_tv):
                continue
            maxima,hist,points = graph(executable,m,n)
            if check_tv:
                d = tv(hist,m,n)
                assert d*m < 26
                distances.append(dict(m=m,n=n,tv=str(d),scaled=str(m*d)))
            if not check_clocks:
                continue
            def label(part):
                u,v = min(part)
                return (u<<n)*((u<<n)+1)//2+(v<<n)-1
            targets = [([i],part) for i,part in enumerate(parts)]
            targets.append((list(range(len(parts))),set.union(*parts)))
            for ids,target in targets:
                labels = {label(parts[i]) for i in ids}
                largest = max(maxima[k] for k in labels)
                assert largest <= 5*n+7
                c = Coefficients(m,n,target)
                observed_times = []
                for (kind,u),(terminal,actual) in points.items():
                    k = c.coefficient(kind,u)
                    predicted = k*(n+1)+1
                    assert predicted == (actual if terminal in labels else inf),(m,n,kind,u,ids)
                    observed_times.append(actual if terminal in labels else None)
                times.append(dict(m=m,n=n,cycles=ids,max_time=largest,
                                  H=observed_times,critical_states=len(points)))
    rows,basins = [],[]
    for r in range(4):
        atoms,basin,moments = limit5(r)
        numerators = [1125*a for a in atoms]
        expected = [[360,394,185,94,0],[360,394,371,0,0],
                    [360,394,92,0,0],[360,394,94,92,185]][r]
        assert numerators == expected
        expected_moments = [(Q(2079,1033),Q(5105,1033)),(Q(2261,1125),Q(211,45)),
                            (Q(712,423),Q(1382,423)),(Q(2723,1125),Q(8879,1125))]
        assert tuple(moments) == expected_moments[r]
        rows.append(dict(phase=r,numerators=list(map(int,numerators)),
                         basin=str(basin),moments=list(map(str,moments))))
    for n in range(1,41):
        base = 5<<n
        count = sum(formula(5,n,{(3,1)}).values())
        p = Q(count,base**4)
        scaled = (p-Q(rows[n%4]['basin']))*2**n
        assert abs(scaled)<1
        basins.append(dict(n=n,count=count,raw=str(p),eligible=str(Q(count,base**4-base)),
                           scaled_error=str(scaled)))
    return dict(polynomials=cases,boundary_checks=times,table=rows,basins=basins,tv=distances)
