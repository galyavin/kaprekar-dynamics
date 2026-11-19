"""Independent recomputation of the signed per-profile weights from the leaf dump.
Reads results/leaves.jsonl (every arithmetic leaf of every residual primary class: lattice, lane, leading profile)
and writes the comparison summary."""
import json, sys, random
from collections import defaultdict, Counter
from fractions import Fraction as Q
from math import gcd
from itertools import product

def det3(c):
    (a, b, d), (e, f, g), (h, i, j) = c
    return a*(f*j - g*i) - b*(e*j - g*h) + d*(e*i - f*h)

leaves = [json.loads(l) for l in open('results/leaves.jsonl')]
print('leaves', len(leaves))
gs = Counter(); phase_mismatch = 0
cls = defaultdict(lambda: defaultdict(Q))     # (label,k,arm,rho,phase) -> profile -> fraction
lanes = defaultdict(set)
for L in leaves:
    cols = L['columns']; D = abs(det3(cols)); g = 0
    for c in cols: g = gcd(g, abs(c[0]))
    gs[g] += 1
    frac = Q(216*g, D)
    phases = tuple(ph for ph in (0, 1) if (L['base'] + 36*ph - L['point'][0]) % g == 0)
    if phases != tuple(L['verifier_phases']): phase_mismatch += 1
    prof = tuple(L['profile'])
    for ph in phases:
        key = (tuple(L['label']), L['k'], L['arm'], L['rho'], ph)
        cls[key][prof] += frac
        lanes[key].add(L['lane'])
print('base step g histogram', dict(gs), '| phase disagreements with verifier', phase_mismatch)
print('classes (entry x phase) =', len(cls), ' (expected 6336*2 = 12672)')
bad = [k for k, d in cls.items() if sum(d.values()) != 1]
print('entries whose fixed-base fibre fractions do not sum to exactly 1:', len(bad))
print('entries with more than one lane:', sum(1 for v in lanes.values() if len(v) != 1))
two = [k for k, d in cls.items() if len(d) > 1]
print('entries with >1 profile:', len(two))
wt = Counter()
for k in two:
    wt[tuple(sorted(cls[k].values()))] += 1
print('weight patterns in two/more-profile entries:', {tuple(map(str, k)): v for k, v in wt.items()})
# lane counts per (k,arm,phase)
cnt = defaultdict(Counter)
for k, v in lanes.items(): cnt[(k[1], k[2], k[4])][next(iter(v))] += 1
print('lane counts per (K mod 6, arm, phase) [S,1,2,L]:')
for key in sorted(cnt): print(' ', key, [cnt[key][l] for l in 'S12L'])
# signed weights
profiles = sorted({p for d in cls.values() for p in d})
lane_of = {}
for k, d in cls.items():
    for p in d: lane_of[p] = next(iter(lanes[k]))
print('distinct profiles', len(profiles))
signed = {}
for kk in (0, 2, 4):
    for ph in (0, 1):
        tot = {arm: defaultdict(Q) for arm in (0, 1)}
        for k, d in cls.items():
            if k[1] == kk and k[4] == ph:
                for p, w in d.items(): tot[k[2]][p] += w
        signed[(kk, ph)] = {p: tot[0][p] - tot[1][p] for p in profiles}
        signed[(kk, ph, 'arm0')] = dict(tot[0]); signed[(kk, ph, 'arm1')] = dict(tot[1])
ref = signed[(0, 0)]
print('signed weights identical in all six (K mod 6, lift) cells:', all(signed[(kk, ph)] == ref for kk in (0, 2, 4) for ph in (0, 1)))
print('\nprofile (s0,sp,sr) | lane | arm5 | arm11 | arm5-arm11   (class units = 1/28512 of box)')
for p in sorted(profiles, key=lambda p: ('S12L'.index(lane_of[p]), [Q(x) for x in p])):
    print(p, lane_of[p], signed[(0, 0, 'arm0')].get(p, 0), signed[(0, 0, 'arm1')].get(p, 0), ref[p])
pos = sum(v for v in ref.values() if v > 0); neg = sum(v for v in ref.values() if v < 0)
print('positive mass', pos, 'negative mass', neg, 'net', pos+neg)
print('sign pattern by lane:', {l: sorted({'+' if ref[p] > 0 else ('-' if ref[p] < 0 else '0') for p in profiles if lane_of[p] == l}) for l in 'S12L'})
longpos = sum(ref[p] for p in profiles if lane_of[p] == 'L')
print('sum over Long profiles', longpos, ' -> contrast', longpos/28512, '=', Q(longpos, 28512), '; 17/7128 =', Q(17, 7128), float(Q(17, 7128)))
print('lane sums of signed weights:', {l: str(sum(ref[p] for p in profiles if lane_of[p] == l)) for l in 'S12L'})
json.dump({'g_histogram': dict(gs), 'entries': len(cls), 'profiles': [[list(p), lane_of[p], str(ref[p]), str(signed[(0,0,'arm0')].get(p,0)), str(signed[(0,0,'arm1')].get(p,0))] for p in profiles],
           'positive': str(pos), 'negative': str(neg), 'identical_all_cells': all(signed[(kk, ph)] == ref for kk in (0, 2, 4) for ph in (0, 1))},
          open('results/weights.json', 'w'), indent=1)

# --- leading chart, clock graphs, plateau slack (independent derivation of P(x), R(x) from the route description)
BOX = ((Q(94040, 10**5), Q(94130, 10**5)), (Q(73770, 10**5), Q(73850, 10**5)), (Q(57940, 10**5), Q(58015, 10**5)),
       (Q(20, 10**5), Q(80, 10**5)), (Q(200, 10**5), Q(300, 10**5)))
def PR(x):
    xa, xb, xc, H4, G5 = x
    n = H4/6; m = (G5 + 5*n)/2
    af = xa - 5*n - 10*m; cf = xc - 3*n
    N = (af + cf - Q(3, 2))/22
    return xc - 3*n - 2*N, xb - 4*n - 2*m - 4*N
U = (Q(940800, 10**6), Q(738087, 10**6), Q(579763, 10**6), Q(484, 10**6), Q(2491, 10**6))
S = (Q(446, 10**6), Q(365, 10**6), Q(342, 10**6), Q(280, 10**6), Q(449, 10**6))
VV = (Q(590, 10**6), Q(-307, 10**6), Q(-447, 10**6), Q(127, 10**6), Q(-296, 10**6))
def theta(x): return Q(121049, 10**6) + sum(v/s*(xi-u) for xi, u, s, v in zip(x, U, S, VV))
verts = list(product(*BOX))
def tau(p, x):
    P, R = PR(x); s0, sp, sr = (Q(t) for t in p); return x[0] + s0 + sp*P + sr*R
slack = {}
for p in profiles:
    if lane_of[p] == 'L': slack[p] = min(tau(p, x) - theta(x) - Q(1, 1000) for x in verts)
    else: slack[p] = min(theta(x) - tau(p, x) - Q(1, 1000) for x in verts)
mn = min(slack.values())
print('\nplateau slack (min over 32 vertices and profiles), my chart:', mn, ' result.json says 53493742523569/96241171257000000 =', Q(53493742523569, 96241171257000000), ' equal:', mn == Q(53493742523569, 96241171257000000))
byl = defaultdict(list)
for p in profiles: byl[lane_of[p]].append(p)
gaps = []
for lo, hi in zip('S12L', '12L'[0:0]) if False else (('S','1'),('1','2'),('2','L')):
    gaps.append(min(tau(b, x) - tau(a, x) for a in byl[lo] for b in byl[hi] for x in verts))
print('adjacent lane gaps (my chart):', [str(g) for g in gaps], ' result.json: 5237/2112000, 5237/3168000, 5237/1683000 ->', [str(Q(5237, 2112000)), str(Q(5237, 3168000)), str(Q(5237, 1683000))])
json.dump({'plateau_slack': str(mn), 'lane_gaps': [str(g) for g in gaps]}, open('results/chart_check.json', 'w'))
