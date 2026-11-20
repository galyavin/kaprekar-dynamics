"""Check the common period-eight root chart and the two late wall feeders."""
from fractions import Fraction as Q
from math import factorial

from common import (add, cycles, determinant, difference_rows, gap_step, insertions,
                    literal, loss, permute, physical_period, root_kick, root_step, spread)

V = [(5,5,5,5), (0,12,12,12), (0,0,26,26), (0,0,27,0)]
CHART = [
    ('++++', (5,13,107,62), 0, (0,1,2,3)),
    ('+++-+++++++-++++', (0,20,114,69), 4, (1,2,3)),
    ('+++-++++++--+++0++-++++-++0+++-0+++-++-++++-++-++++0', (0,8,102,57), 12, (2,3)),
    ('++-++++-++++++-++++-++++++-++++-++++++-++++-++++', (0,0,71,28), 18, (2,3)),
    ('++-+0++-++++-++0++-++++-0+++-+-0++0-++-+0++-++-+-++0', (0,0,45,2), 13, (3,)),
]
MOVES = [tuple(-x for x in V[0]), tuple(-x for x in V[1]), (0,-8,-31,-29),
         tuple(-x for x in V[2]), tuple(-x for x in V[3])]
SIGN = {'+': 1, '0': 0, '-': -1}


def trace(x, count, step):
    path = [x]
    for _ in range(count):
        path.append(step(path[-1]))
    return path


def check_chart(shift):
    bounds = []
    assert abs(determinant(V)) == 42120
    for (word, apex, total, indices), move in zip(CHART, MOVES):
        base = add(apex, (0,0,shift,0))
        rays = [V[i] for i in indices]
        x = base
        time_loss = 0
        bound = 0
        for character in word:
            z = x[-1]
            e = SIGN[character]
            assert (z*e >= 1) if e else z == 0
            assert all(ray[-1]*e >= 0 if e else ray[-1] == 0 for ray in rays)
            if e < 0:
                assert all(ray[-1] == 0 for ray in rays)
            bound = max(bound, spread(x)+sum(map(spread, rays)))
            h = loss(z)
            rays = [add(root_kick(add(x, ray), e, loss(z+ray[-1]) if e < 0 else h),
                        root_kick(x,e,h), -1) for ray in rays]
            x = root_kick(x,e,h)
            time_loss += h
        assert x == add(base, move) and rays == [V[i] for i in indices]
        assert time_loss == total
        bounds.append(max(bound, spread(x)+sum(map(spread, rays))))
    suffix = trace((0,0,18+shift,2), 43 if shift == 0 else 40, root_step)
    assert max(map(spread, suffix)) <= max(bounds)
    used = sum(loss(x[-1]) for x in suffix[:-1])
    print(f'period-eight chart {shift=}: det=42120 spreads={bounds} suffix_loss={used}')
    if shift == 0:
        assert suffix[35] == suffix[-1] and len(set(suffix[:-1])) == 43
        assert min((0, *x) for x in suffix[35:-1]) == (0,-1,2,-1,0)
        assert all(loss(x[-1]) == 0 for x in suffix[35:-1])
        assert max(bounds) <= 179 and used == 10
    else:
        assert suffix[-1] == (0,0,1,-2) and max(bounds) <= 178 and used == 10
    return suffix[-1]


def prepare(n, B, pi, x, others):
    _, labels, *_ = cycles(n)
    g = [None]*len(pi)
    g[0] = (B*pi[0]-Q(sum(x), len(x)+1)).__floor__()
    for j, value in zip(labels, x):
        g[j] = g[0]+value
    for j, weight in others.items():
        g[j] = (B*weight).__floor__()
    return g


def close(g, B):
    missing = g.index(None)
    assert g.count(None) == 1
    g[missing] = B-sum(v for v in g if v is not None)
    return tuple(g)


def off_centre():
    for n, alpha, divisor in [(15,Q(47,400),100000), (17,Q(14,135),100000)]:
        _, labels, *comps = cycles(n)
        m = (n+1)//2
        beta = (1-5*alpha)/(m-5)
        pi = tuple(alpha if j in (0,*labels) else beta for j in range(m))
        assert permute(pi) == pi and sum(pi) == 1
        assert min(difference_rows(pi)[1:]) > Q(8*191,divisor)
        assert (m+1)*(min(pi)-Q(191,divisor)) > 1+Q(1,divisor)
        assert Q(1,m)-alpha-Q(191,divisor) > Q(1,200)
        assert min(pi)-Q(192,divisor) >= Q(1,11)
        assert all(8%len(part) == 0 for part in comps)
        B = 200003
        x = (0,8,49,31)
        for ray in V[2:]:
            x = add(x, ray)
        free = [j for part in comps for j in part]
        g = close(prepare(n,B,pi,x,{j: beta for j in free[:-1]}),B)
        count = len(CHART[2][0])+len(CHART[3][0])+len(CHART[4][0])+35
        end = literal(g,count)
        assert physical_period(end,8) == 8
        print(f'n={n} off-centre: delta=1/200 eta={n}!/(11^{m}*{divisor}^{m-1}); {count+8} literal transitions')


def direction(size,j):
    return tuple(i==j for i in range(size))


def coast_and_section(n, root, low, low_residue, highs, fixed, drift_length, section_length):
    _, labels, *_ = cycles(n)
    g = [0]*((n+1)//2)
    for i, value in zip(labels, root):
        g[i] = value
    g[low] = low_residue
    g = tuple(g)
    allowed = (2,) if n==17 else (2,3)

    def section_step(state):
        phase = len(section)-1
        active = len(allowed)+1 if locations[phase] == low else 1
        return gap_step(state, range(1, active+1))

    locations = [low]
    for _ in range(max(drift_length,section_length)+1):
        v = 2*locations[-1]
        locations.append(min(v,n-v))
    drift = trace(g, drift_length, lambda x: gap_step(x,(1,)))
    section = [g]
    for t in range(section_length):
        raw = section_step(section[-1])
        section.append(raw)
    # At low-facing phases, high or fixed directions must be absent from
    # every active secondary correction. Verify the actual tested rows.
    for t, location in enumerate(locations[:max(drift_length,section_length)]):
        if location == low:
            for row in allowed:
                assert all(difference_rows(direction(len(g),j))[row-1] == 0
                           for j in (*highs,*fixed))
        highs = tuple(min(2*j,n-2*j) for j in highs)
    d = add(drift[-1],g,-1)
    if n == 17:
        assert d == (12,12,12,-15,12,-15,-15,-15,12)
        assert min(difference_rows(x)[1] for t,x in enumerate(drift[:-1]) if locations[t]==low)==-6
    else:
        assert d == (42,42,42,-98,42,42,-98,0,42,-98,42)
        assert min(min(difference_rows(x)[1:3]) for t,x in enumerate(drift[:-1]) if locations[t]==low)==-3
    # Uniform physical margins. At each phase the finite trace is
    # affine in the high offsets and, on complete coast blocks, the
    # excess over the declared low offset. Extremizing each linear
    # guard over that surrounding box is conservative.
    cutoff = 10**6 if n==17 else 10**8
    alpha = Q(21,200) if n==17 else Q(87,1000)
    error = Q(191,10**6) if n==17 else Q(2500,10**8)
    high_min,high_max = (Q(16,1000),Q(19,1000)) if n==17 else (Q(5,1000),Q(10,1000))
    excess_bound = Q(1282,10**6) if n==17 else Q(105,100000)

    def prove_physical(trace_, coast):
        high_positions = tuple(j for j in range(len(g))
                               if j not in (0,*labels,low,*fixed))
        for phase,state in enumerate(trace_[:-1]):
            low_position = locations[phase]
            boxes = {j:(high_min,high_max) for j in high_positions}
            boxes.update({j:(Q(25,1000),Q(27,1000)) for j in fixed})
            boxes[low_position] = (Q(0),excess_bound) if coast else (Q(0),Q(0))
            for j,offset in enumerate(state):
                minimum = alpha-error+boxes.get(j,(0,0))[0]
                min_gap = Q(104 if n==17 else 86,1000)
                assert minimum>=min_gap and (minimum-min_gap)*cutoff+offset>0
            comparisons = difference_rows(state)
            for i,value in enumerate(comparisons,1):
                if i==1:
                    assert all(difference_rows(direction(len(g),j))[0]==0 for j in boxes)
                    continue
                if i in allowed and low_position==low:
                    if coast:
                        assert value + (27 if n==17 else 140) > 1
                    continue
                lead = sum(lo if j>=len(g)-i else -hi if j<i else 0
                           for j,(lo,hi) in boxes.items())
                assert lead>=0 and lead*cutoff+value>1, (n,phase,i,lead,value)
            high_positions = tuple(min(2*j,n-2*j) for j in high_positions)

    prove_physical(drift, True)
    prove_physical(section, False)
    return drift, section, locations


def check_section(n, section, high_indices, fixed_indices, start, expected_outside,
                  expected_contacts, length, peak_t, peak_offset):
    _, root_labels, *_ = cycles(n)
    all_outside, contacts = [], []
    B = 1000000 if n == 17 else 100000000
    m = len(section[0])
    A = 105350 if n == 17 else 8730000
    high_values = ([17500, 17600, 0] if n==17 else [690000, 0])
    fixed_values = [2550000] if n==21 else []
    g = add(section[0], (A,)*m)
    chosen = [(j,v) for j,v in zip(high_indices,high_values) if v] + list(zip(fixed_indices,fixed_values))
    for j,v in chosen:
        g = add(g,direction(m,j),v)
    final_high = high_indices[-1]
    g = add(g,direction(m,final_high),B-sum(g))
    assert sum(g) == B
    end_state = []

    def check(t, actual):
        d = difference_rows(actual)
        if min(d[1:]) < -1:
            all_outside.append((t, d[1:len(expected_contacts[0][1])+1]))
        if t >= start:
            end_state.append(actual)
            local = d[1:len(expected_contacts[0][1])+1]
            if max(local) <= 1:
                contacts.append((t,local))

    end = literal(g,start+length,check)
    assert section[start] == section[start+length]
    assert all(len(part)==1 or length%len(part)==0 for part in cycles(n)[2:])
    assert end == end_state[0] and physical_period(end,length) == length
    assert len(set(tuple(state[j]-state[0] for j in root_labels) for state in end_state))==length
    assert all_outside == expected_outside if isinstance(expected_outside,list) else len(all_outside)==expected_outside
    assert contacts == expected_contacts
    candidates = [(abs(Q(sum(state[:j]),B)-Q(j,m)), t, j)
                  for t,state in enumerate(end_state,start) for j in range(1,m)]
    peak, time, prefix = max(candidates)
    j = 5 if n==17 else 6
    assert (time,prefix,peak) == (peak_t,j,Q(j,m)-Q(j*A+peak_offset,B))
    # All-base comparison: each prefix/suffix has finite offsets plus a
    # collection of independent macroscopic high-direction terms.
    moving = list(high_indices)
    active_sets = {}
    for t in range(start+length):
        active_sets[t] = set(moving)|set(fixed_indices)
        moving = [min(2*i,n-2*i) for i in moving]
    lower_f, upper_f, high_floor, fixed_floor, cutoff = (
        (Q(1,200),Q(1,150),Q(16,1000),Q(0),10**6) if n==17 else
        (Q(35,10000),Q(4,1000),Q(5,1000),Q(25,1000),10**8))
    for t in range(start,start+length):
        for count in range(1,m):
            for indices in (range(count), range(m-count,m)):
                high = sum(i in active_sets[t] and i not in fixed_indices for i in indices)
                fixed = sum(i in fixed_indices for i in indices)
                edge = lower_f if j>=count else upper_f
                leading = (j-count)*edge+high*high_floor+fixed*fixed_floor
                constant = sum(section[t][i] for i in indices)-peak_offset
                assert leading >= 0 and leading*cutoff+constant >= 0, (n,t,count,list(indices),leading,constant)
    assert sum(section[peak_t][:j])==peak_offset
    assert not (active_sets[peak_t]&set(range(j)))
    print(f'n={n} section: {len(all_outside)} outside states; contacts={contacts}; '
          f'peak={j}/{m}-(%dA+%d)/B; {start+length} literal steps' % (j,peak_offset))


def prefix_profile(n, pi, errors, D, min_gap, phases):
    m = len(pi)
    assert sum(pi)==1
    for _ in range(phases):
        assert min(v-Q(e,D) for v,e in zip(pi,errors))>min_gap
        for i,margin in enumerate(difference_rows(pi),1):
            if i>1:
                assert margin>Q(sum(errors[:i])+sum(errors[-i:]),D)
        pi,errors = permute(pi),permute(errors)
    assert (m+1)*min_gap>1+Q(1,D)


def source_control(n, source, root_clock, drift, low_origin, low_end, mod, residue, terminal_clock, period):
    scale = 10**6 if n==17 else 10**8
    B = 2*scale+123
    if n==17:
        pi = [Q(21,200)]*9
        others = {3:Q(123,1000),5:Q(123,1000),7:Q(106,1000)}
    else:
        pi = [Q(87,1000)]*11
        others = {3:Q(94,1000),6:Q(88,1000),7:Q(113,1000)}
    partial = prepare(n,B,pi,source,others)
    initial = close(partial.copy(),B)
    formal = trace(initial,root_clock,lambda g:gap_step(g,(1,)))[-1]
    c0 = formal[low_end]-formal[0]
    partial[low_origin] += (residue-c0)%mod
    g = close(partial,B)
    capture = literal(g,root_clock)
    c = capture[low_end]-capture[0]
    blocks,remainder = divmod(c-residue,mod)
    assert remainder==0 and blocks>=2
    block = len(drift)-1
    change = add(drift[-1],drift[0],-1)
    assert literal(capture,block) == add(capture,change)
    last = add(capture,change,blocks-1)
    section = literal(last,block)
    assert section == add(capture,change,blocks)
    count_outside = []

    def book(t, state):
        if min(difference_rows(state)[1:]) < -1:
            count_outside.append(t)

    final = literal(section,terminal_clock+period,book)
    assert len(count_outside)==(1 if n==17 else 17)
    assert physical_period(final,period)==period
    print(f'n={n} source to capture plus first/last coast blocks and section: '
          f'{root_clock+2*block+terminal_clock+period} literal steps; {B=}; {blocks=}')


def verify():
    check_chart(0)
    off_centre()
    root17 = check_chart(-1)
    path = trace(root17,61,root_step)
    assert path[-1] == root17 and len(set(path[:-1]))==61
    assert sum(loss(x[-1]) for x in path[:-1])==15
    drift17, section17, _ = coast_and_section(17,root17,7,23,(3,5,6),(),244,210)
    check_section(17,section17,(3,5,6),(),202,[(200,(-3,))],[(208,(1,))],8,205,48)
    pi17 = tuple(Q(106,1000) if j==7 else Q(123,1000) if j in (3,5,6) else Q(21,200)
                 for j in range(9))
    err17 = tuple(91 if j in (3,5,6,7) else 191 for j in range(9))
    prefix_profile(17,pi17,err17,10**6,Q(104,1000),4)
    assert Q(1,1000)-Q(3*191+91,10**6) > 0
    assert Q(2,1000)-Q(4*191+2*91,10**6) > 0
    assert Q(19,1000)-Q(4*191+4*91,10**6) > 0
    assert Q(21,200)-Q(191,10**6)>Q(104,1000)
    # Bound F from the root at the start of the coast and its anchor gain.
    assert Q(1,1000)-Q(191+91,10**6) == Q(718,10**6)
    assert Q(1,1000)+Q(191+91,10**6) == Q(1282,10**6)
    assert Q(18,1000)-Q(191+91+1282+20,10**6) > Q(16,1000)
    assert Q(18,1000)+Q(191+91,10**6) < Q(19,1000)
    F17_min = Q(1,9)-Q(21,200)-Q(191,10**6)-Q(12,27)*Q(1282,10**6)
    F17_max = Q(1,9)-Q(21,200)+Q(191,10**6)
    assert F17_min > Q(1,200) and F17_max < Q(1,150)
    assert 5*F17_min-Q(48,10**6)>Q(1,40)
    assert Q(104,1000)-Q(1,10**6)>Q(1,11)

    source = (1961,320,1323,1639,770,1113)
    intervals = [(0,6),(280,304),(956,1034),(1633,1675),(2620,2764),(3609,3897)]
    route = insertions(source,intervals,4486)
    assert abs(determinant(route.rays)) == 224501760
    assert route.lc_budget==(5068,1729,2181)
    root21 = route.terminal
    path = trace(root21,326,root_step)
    assert path[-1] == root21 and len(set(path[:-1]))==326
    assert sum(loss(x[-1]) for x in path[:-1])==98
    drift21, section21, _ = coast_and_section(21,root21,9,140,(3,6),(7,),978,1213)
    check_section(21,section21,(3,6),(7,),1201,17,[(1203,(1,1)),(1206,(1,0))],12,1207,276)
    pi21 = tuple(Q({3:94,6:88,7:113,9:96}.get(j,87),1000) for j in range(11))
    prefix_profile(21,pi21,(2500,)*11,10**8,Q(86,1000),3)
    assert Q(1,1000) > Q(25000,10**8)
    assert Q(87,1000)-Q(2500,10**8)>Q(86,1000)
    assert Q(1,1000)-Q(5000,10**8) == Q(95,100000)
    assert Q(1,1000)+Q(5000,10**8) == Q(105,100000)
    assert Q(7,1000)-Q(5000,10**8)-Q(105,100000) > Q(5,1000)
    assert Q(9,1000)+Q(5000,10**8) < Q(10,1000)
    assert Q(26,1000)-Q(5000,10**8)-Q(3,10)*Q(105,100000) > Q(25,1000)
    assert Q(26,1000)+Q(5000,10**8) < Q(27,1000)
    F21_min = Q(1,11)-Q(87,1000)-Q(2500,10**8)-Q(3,10)*Q(105,100000)
    F21_max = Q(1,11)-Q(87,1000)+Q(2500,10**8)
    assert Q(35,10000)<F21_min<=F21_max<Q(4,1000)
    assert 6*F21_min-Q(276,10**8)>Q(1,50)
    counts = (1,1,2,1)
    x17 = (0,8,48,31)
    for ray,param in zip(V,counts):
        x17 = add(x17,ray,param)
    clock17 = 40+sum(count*len(CHART[i][0]) for count,i in zip(counts,(0,1,3,4)))+len(CHART[2][0])
    source_control(17,x17,clock17,drift17,7,7,27,23,202,8)
    x21 = source
    for i in (0,2,4):
        x21 = add(x21,route.rays[i])
    clock21 = 4486+sum(intervals[i][1]-intervals[i][0] for i in (0,2,4))
    source_control(21,x21,clock21,drift21,6,9,140,140,1201,12)
    for n, factor, D in ((17,11,10**6),(21,12,10**8)):
        m = (n+1)//2
        assert Q(factorial(n),factor**m*D**(m-1)) > 0
        print(f'n={n} wall original-word liminf >= {n}!/({factor}^{m}*{D}^{m-1})')


if __name__ == '__main__':
    verify()
