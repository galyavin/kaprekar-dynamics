"""Exact pairing, row and original-word weight bounds for the frozen sources."""
import json
import math
from fractions import Fraction


def setup(n, weights):
    m = (n+1)//2
    orbit = [-1]*m
    orbit[0] = 0
    reps = [1]
    for j in range(1,m):
        if orbit[j]>=0:
            continue
        o = len(reps)-1 if j==1 else len(reps)
        if j!=1:
            reps.append(j)
        t = j
        while orbit[t]<0:
            orbit[t] = o
            t = min(2*t,n-2*t)
    offsets = [weights[i] for i in orbit]
    return orbit,reps,offsets


def weight(n, B, L, g):
    """Distinct-word lower and split-count upper weight for the entire box."""
    k = (n-1)//2
    variance = Fraction(L*(L+1),3)
    result = []
    for delta in (-1,1):
        counts = [v+delta for v in g[1:]]
        polynomial = math.prod(counts)*(g[0]-variance*sum(Fraction(1,v) for v in counts))
        numerator = math.factorial(n)*(2*L+1)**k*polynomial
        assert numerator.denominator==1
        mass = numerator/(B**n-B)
        result.append({"numerator": numerator.numerator,"weight": f"{float(mass):.12e}"})
    return result


def rank(matrix):
    rows = [list(map(Fraction, row)) for row in matrix]
    pivot = 0
    for column in range(len(rows[0])):
        index = next((j for j in range(pivot,len(rows)) if rows[j][column]),None)
        if index is None:
            continue
        rows[pivot],rows[index] = rows[index],rows[pivot]
        for j in range(pivot+1,len(rows)):
            if rows[j][column]:
                ratio = rows[j][column]/rows[pivot][column]
                rows[j] = [a-ratio*b for a,b in zip(rows[j],rows[pivot])]
        pivot += 1
    return pivot


def projected_row_rank(n,orbit,row):
    m = len(orbit)
    complement = [j for j in range(m) if orbit[j]]
    phase = list(range(m))
    rotations = []
    while True:
        rotations.append([
            int(phase[j]>=m-row)-int(phase[j]<row) for j in complement
        ])
        phase = [min(2*j,n-2*j) for j in phase]
        if phase == list(range(m)):
            break
    assert all(sum(col)==0 for col in zip(*rotations))
    return rank(rotations)


def geometry(n,B,L,relative_weights):
    k,m = (n-1)//2,(n+1)//2
    orbit,reps,w = setup(n,relative_weights)
    p = [Fraction(1,m)+Fraction(m*v-sum(w),1000*m) for v in w]
    assert sum(p)==1
    assert all(p[j]==p[min(2*j,n-2*j)] for j in range(1,m))
    denominators=math.lcm(*(v.denominator for v in p))
    assert B%denominators==0
    centre=[int(v*B) for v in p]
    q=[None]+[sum(p[i:]) for i in range(1,m)]
    pairs=[(q[i],1-q[m-i]) for i in range(1,m)]
    adjacent = [((1,1),pairs[0])]+list(zip(pairs,pairs[1:]))
    pairing = min(min(a)-max(b) for a,b in adjacent)
    projected_peak = max(abs(q[i]-Fraction(m-i,m)) for i in range(1,m))
    rows=[sum(p[m-i:])-sum(p[:i]) for i in range(1,m//2+1)]
    positive_rows=[a for a in rows[1:] if a]
    assert B*min(positive_rows)>m*L
    assert B*pairing-2*k*L-1>0
    bound = min(centre[0]-k*L, min(centre[1:])-L)
    assert bound>=2
    zero_rows = [i for i,a in enumerate(rows,1) if i>=2 and not a]
    assert len(zero_rows)==1
    pinning_rank=projected_row_rank(n,orbit,zero_rows[0])
    assert pinning_rank == (0 if n==97 else 23)
    return {
        "length":n,"base":B,"half_width":L,"profile_denominator":denominators,
        "orbit_representatives":reps,"orbit_sizes":[orbit.count(i) for i in range(len(reps))],
        "centre":centre,"min_gap_box_lower_bound":bound,
        "centre_homogeneous_pairing":str(pairing),
        "box_exact_pairing_lower_bound":str(B*pairing-2*k*L-1),
        "centre_secondary_rows":[str(a) for a in rows[1:]],
        "box_positive_secondary_lower_bound":str(B*min(positive_rows)-m*L),
        "zero_centre_secondary": zero_rows,
        "rotated_complement_row_rank": pinning_rank,
        "zero_centre_rows_box_min_and_max":[-m*L,m*L],
        "centre_expected_Z":str(projected_peak),"original_weight_bounds":weight(n,B,L,centre),
    }


if __name__=="__main__":
    results=[geometry(97,256025000,256,[0,1]),geometry(119,256020000,256,[0,1,1,1,2])]
    print(json.dumps(results,indent=2))
