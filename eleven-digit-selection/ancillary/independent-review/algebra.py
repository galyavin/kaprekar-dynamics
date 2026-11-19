"""Own exact re-derivation of the structural twelve-shift constants from the README's
stated inputs: q_R(B,p,r), J, slacks, and the two checkpoint matrices/drifts."""
from fractions import Fraction as Q
out = []
dB = 12
# Delta q_R with N,m,n fixed: q_R = (3(B-1)/2 + r1 - p, r, p, B-1+z-p, r-(B+H-1)/2)
def dq(dp_, dr_): return (Q(3, 2)*dB - dp_, dr_, dp_, dB - dp_, dr_ - Q(dB, 2))
dp_ = 18  # fixes anchor a
assert dq(dp_, 0)[0] == 0
# slacks H=2(b-e)-B+1, G=B-c-d-1, u=q1+q3-3(B-1)/2 invariant for every d_r:
for dr_ in range(-30, 31):
    a, b, c, d, e = dq(dp_, dr_)
    assert 2*(b - e) - dB == 0 and dB - c - d == 0 and a + c - Q(3, 2)*dB == 0
dJ = lambda dr_: 2*dB + 2*dp_ - 4*dr_          # J(q2) = 2B + 2p - 4r + kappa
v = (5, 4, 3, 2, 1)
Mreg = [[Q(17,18),Q(-14,9),Q(10,9)],[Q(-4,9),Q(-4,9),Q(17,9)],[Q(-2,3),Q(1,3),Q(4,3)],
        [Q(4,9),Q(-14,9),Q(10,9)],[Q(-17,18),Q(-4,9),Q(17,9)]]
dreg = [Q(-1),Q(-11,10),Q(-3,5),Q(-1),Q(-11,10)]
Mexc = [[Q(17,18),Q(-14,9),Q(10,9)],[Q(37,18),Q(-31,9),Q(8,9)],[Q(3),Q(-3),Q(-1)],
        [Q(25,9),Q(-20,9),Q(-14,9)],[Q(25,18),Q(-10,9),Q(-7,9)]]
dexc = [Q(-1),Q(-7,5),Q(-3,10),Q(1,5),Q(1,10)]
for name, M, dr in (('reg', Mreg, dreg), ('exc', Mexc, dexc)):
    # unknown d_r: centred lift of row i equals drift_i * s, with s fixed by row 1 (drift -1)
    # row1: M1.(12,18,x) - 2*5 = -s  -> s(x); row2 gives linear equation in x
    def row(i, x): return M[i][0]*12 + M[i][1]*18 + M[i][2]*x - 2*v[i]
    # s(x) = -row(0,x) is affine; solve row(1,x) = dr[1]*s(x)
    a0, a1 = -row(0, 0), -(row(0, 1) - row(0, 0))
    b0, b1 = row(1, 0), row(1, 1) - row(1, 0)
    x = (dr[1]*a0 - b0) / (b1 - dr[1]*a1)
    s = a0 + a1*x
    lift = [row(i, x) for i in range(5)]
    assert all(lift[i] == dr[i]*s for i in range(5))
    L = Q(10)*dJ(x)/36
    out.append((name, x, s, L, lift))
    assert s == L + 10
for o in out: print(o)
x = out[0][1]
w = dq(dp_, x)
print('w =', w, ' dJ =', dJ(x), ' r2 shift mod 72 =', dJ(x) % 72, ' T shift = L - s =', out[0][3] - out[0][2])
print('half-return dJ:', 2*(-10) + 2*(-2) - 2*2 + 4*(-2))
print('lift of p-2r:', 18 - 2*x, ' dp0 = dr2/2 - dB mod 36:', (Q(36, 2) - 12) % 36)
