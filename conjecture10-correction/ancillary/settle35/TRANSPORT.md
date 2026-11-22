# Closed root dynamics and physical transport

Fix an odd length \(n=2k+1\ge7\) and put \(m=k+1\). For the descending sorted differences \(q_i=a_i-a_{n+1-i}\), use \(g_0=B-q_1\), \(g_i=q_i-q_{i+1}\) for \(1\le i<k\), and \(g_k=q_k\). The fold permutation \(P\) fixes 0 and moves \(j>0\) to \(f(j)=\min(2j,n-2j)\). Let \(J\) include 0 and the length-\(r\) orbit of 1. Positive pairing means \(q_k>0\) and decreasing borrowed levels. At level \(i\) the pair is
\[
(q_i-\mathbf1_{i=k},\ B-q_{k+1-i}-\mathbf1_{i<k}),
\]
and the minimum of the earlier pair must be at least the maximum of the later pair. Define
\[
d_i(g)=\sum_{a=0}^{i-1}(g_{k-a}-g_a),\qquad
h_i=(-1-d_i)_+,\qquad e=\operatorname{sgn}(d_1),\qquad 1\le i\le p_0=\lfloor m/2\rfloor.
\]
The original chamber \(C\) consists of positively paired physical states with every secondary \(h_i=0\).

## Reduction and exported loss

The prefix positions \(a\) map to \(2a\) and suffix positions \(k-a\) to \(2a+1\) for \(0\le a<i\). Membership in \(J\) is preserved by \(P\). Thus \(d_i\) involves only roots exactly when all of \(1,\ldots,2i-1\) belong to \(J\). The root-only rows form the initial block \(1,\ldots,I\). At a proper length (\(J\ne\{0,\ldots,k\}\)) the first nonroot slot is \(2I+1\). In fact \(I=p_0\) would cover every slot for even \(m\). For odd \(m\) the only potentially uncovered slot is \(k=2p_0\), which is already in the root orbit.

Use \(C^*\) for positively paired physical states with \(d_i\ge-1\) for all \(i>I\). P10's paired update, obtained by sorting each ordered level, gives on \(C^*\)
\[
g^+=Pg+e(\mathbf e_0-2\mathbf e_1+\mathbf e_2)
+\sum_{i=1}^{I} h_i(-\mathbf e_{2i-2}+\mathbf e_{2i-1}+\mathbf e_{2i}-\mathbf e_{2i+1}).
\]
Indeed the sorting correction at slots \(2i-2,2i-1\) is \((-\Delta h_i,\Delta h_i)\) with \(\Delta h_i=h_i-h_{i-1}\) and \(h_0=0\). Summing over rows yields the displayed formula. Every correction index but \(2I+1\) is a root. The enabled \(h_i\) depend only on differences from \(g_0\), so subtracting the updated anchor gives a closed, base-independent integer map. The root sum increases by \(h_I\). For a complementary label initially at slot \(j\), its position at time \(t\) is \(f^t(j)\) and its gap value is
\[
g_j(0)-\sum_{0\le s<t:\ f^{s+1}(j)=2I+1}h_I(s).
\]
This charges every nonnegative export to exactly the arriving label and holds up to the first failure of the enlarged guard. [P10, paired update](https://doi.org/10.5281/zenodo.24408723)

At length 35 the non-anchor root order is \(1,2,4,8,16,3,6,12,11,13,9,17\). The complementary cycles are \((5,10,15)\) and \((7,14)\), and \(I=2\). In the listed root order let \(x=(x_1,\ldots,x_{12})\) consist of the offsets from the anchor. Put \(y=x_{12}\), \(z=x_5+x_{12}-x_1\), \(e=\operatorname{sgn}(y)\), \(h_1=(-1-y)_+\), \(h_2=(-1-z)_+\) and \(a=h_1-e\). The map \(T_{35}\) applies the right shift \(Sx\) and then these additions.

| Coordinate | 1 | 2 | 3 | 6 | 4, 5, 7, 8, 9, 10, 11, 12 |
|---|---:|---:|---:|---:|---:|
| Addition | \(-3e+2h_1\) | \(2h_1-h_2\) | \(a+h_2\) | \(-e+h_2\) | \(a\) |

The labels from 15, 10 and 5 pay the first three exports, respectively, then repeat this order every third step. The labels at 7 and 14 rotate without loss. The guard checks \(g_{17}>0\), all adjacent borrowed levels and \(d_3,\ldots,d_9\ge-1\). On a zero-export recurrent root, \(d_2\ge-1\) at every recurrent phase. Primary \(h_1\) may remain nonzero.

## Conditional all-base transport

Suppose a proper length admits an invariant normalized root-equal profile \(\pi\) with \(d_i(\pi)>0\) for all \(i>I\). This is complement-strictness. Choose \(\pi\) rational and noncentral. The strict linear inequalities permit rational choices in the invariant root-equal space. Mixing with the constant profile preserves the strict margins, so arrange \(b_*=\min_i\pi_i>1/(m+1)\). If there are no complementary rows, choose a noncentral profile by perturbing a complementary orbit. Set
\[
\Delta=b_*-\frac1{m+1},\qquad
\gamma=\min_{i>I}d_i(\pi),\qquad
A=\max_{1\le i\le k}\left|\sum_{j=i}^k\pi_j-\frac{m-i}{m}\right|>0
\]
with \(\gamma=1\) when there are no complementary rows.

Suppose that for every integer \(L\ge L_0\) there are at least \(a_0L^r\) distinct relative roots with lifetime spread at most \(C_0L\), total export at most \(C_0L\), and eventual periodicity with recurrent \(h_I=0\). Take \(C_0\ge1\), \(E=m+2C_0+1\), any integer \(D\ge4E/\min(\Delta,\gamma/m,A/m)\), and \(L=\lfloor B/D\rfloor\).

For each root \(x\) choose the anchor from the \(L\) integers starting at
\[
\left\lfloor B\pi_0-\frac{\sum_{j=1}^r x_j}{r+1}\right\rfloor .
\]
The root gaps follow by adding the offsets. Choose all but one complementary gap independently from \(L\)-integer intervals starting at their profile floors. The last complementary gap sets the exact integer sum \(B\). There are \(L^{k-r}\) distinct lifts per root. Initially the root mean error is at most \(L\), so the root sum error is at most \((r+1)L\). The independent complementary errors are at most \(L\) and the dependent one is at most \(mL\).

Use the paired formula to generate a candidate integer history. It preserves total sum. The root mean moves from its initial level by at most \(C_0L/(r+1)\), so the spread gives root error at most \((1+2C_0)L\) from the profile. Each complementary label incurs additional error at most \(C_0L\), since the profile is constant on its orbit and the label loses only its own share of the total. Thus every gap throughout the candidate path differs from \(B\pi_j\) by at most \(EL\). The gap minimum exceeds \(B/(m+1)\) by at least \(3B\Delta/4\). Also \(d_i(g)\ge B\gamma-mEL\ge3B\gamma/4\) for \(i>I\).

The gap minimum gives physicality and full height for large \(B\). The paired guard holds if \((m+1)\min_jg_j\ge B+1\). For adjacent borrowed levels the same-column comparisons have lower bound \(\min_jg_j-1\). The two other comparisons are prefix or suffix sums of \(m+1\) gap terms minus \(B\), with integer offset at least \(-1\). Hence the sufficient guard holds once \(3(m+1)B\Delta/4\ge1\). The paired formula identifies the entire candidate history with literal subtraction by induction. At each phase the suffix gaps show radius at least \(A-mEL/B\ge3A/4\).

After export ceases the root mean is constant, so the relative cycle also repeats the absolute root. The complement permutes, and a common multiple of the relative period and the complementary permutation order returns the entire gap state. A gap state determines the next whole word, and a word determines its gap state, yielding the same least period. This proves a macroscopic \(C^*\) terminal. If all recurrent inner secondary losses \(h_2,\ldots,h_{I-1}\) vanish, every cyclic phase is in \(C\). At 35 this additional condition is empty.

Distinct gap sources have disjoint original words. Their exact distinct-digit multiplicity is \(n!g_0\prod_{j=1}^k(g_j-1)\), from the smallest digit and the positive splits of opposite adjacent gaps. For sufficiently large \(B\) this is at least \(n!(B/[2(m+1)])^m\). Since \(L\ge B/(2D)\) for \(B\ge2D\), dividing the word count for the \(a_0L^k\) source states by \(B^n-B\) gives
\[
\Pr_{n,B}\{Z_B>A/2\}
\ge \frac{a_0\, n!}{(2D)^k[2(m+1)]^m}
\]
at every sufficiently large integer base. The exact sum and integer choices use no base congruence.

At 35 a convenient profile is \(3/55\) on roots and \(16/275\) off the root. All nonroot secondary rows are positive. Row \(d_2\) involves only roots, so it vanishes on every invariant root-equal vector and demonstrates non-strictness. [Length-35 constants and access certificate](FAMILY.md)
