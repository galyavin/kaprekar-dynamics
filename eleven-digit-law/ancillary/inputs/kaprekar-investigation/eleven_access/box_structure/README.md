# Twelve-base transport on the eleven-digit box

Fix the whole Task 026 box \(C\), with all three cohorts and \(B\equiv5\pmod 6\). Each conditional population consists of uniform original nonconstant eleven-digit words with their own base's exact word multiplicities. Let \(x\) be the normalized original chart and \(T_B\) the first cyclic-word time. The structural proof below applies at all sufficiently large bases on the common chart domain of the source and its lift, including the cohort-0 alignment. Its paired population consequence is an \(O(B^{-1})\) bounded-Lipschitz discrepancy on \((x,T_B/B)\) for each twelve-base shift in both arms.

**Outcome (a).** The joint law of \((x,T_B/B)\) conditional on the whole box has a proper limit in each of the classes 5 and 11 modulo 12, for every bounded Lipschitz test. The exact graph weights distinguish the limits by the frozen test contrast \(17/7128\).

## Slack-preserving lift and half-return clock

Write \(q_0=(a,b,c,d,e)\) and \(q_R\) for the final late-return section. Set \(p=(q_R)_3\), \(r=(q_R)_2\), and use \(H=2(b-e)-B+1\), \(G=B-c-d-1\), and \(u=q_1+q_3-3(B-1)/2\) at their respective prefix sections. Holding the three slacks fixed under \(\Delta B=12\) and prescribing increments \(d_p,d_r\) at \(q_R\) gives
\[
 \Delta q_0=(18-d_p,\ d_r,\ d_p,\ 12-d_p,\ d_r-6).
\]
Indeed, \(H\) is invariant because \(b-e\) increases by six, \(G\) is invariant because \(c+d\) increases by twelve, and \(u\) is invariant because \(q_1+q_3\) increases by eighteen. The cohort, overshoot \(z\), pin \(r_1\), and counts \(n,m,N\) therefore stay fixed. All corresponding updates through the short switch take exactly the same time in the two bases. Work with a lift that fixes the original anchor \(a\), so \(d_p=18\). Cohort 0 uses \(q_0+\alpha_0\) at the same base with its five-step clock difference. That translation commutes with the lift.

At the first paired block, \(J(q_2)=2B+2p-4r+\kappa\), where \(\kappa\) is fixed by the unchanged early label. Thus
\[
 \Delta J=60-4d_r.
\]
The six possible twenty-step words are squares of ten-step half-words. Direct composition of the borrowed digit rows gives displacement \(\Delta_2/2=(-10,-2,0,2,-2)\) for **each half** on its pinned section. It is enough to check the six half-words on the 88 short-switch rows at \(q_2\) and the two extreme section remainders. The ten-step map is affine and the latter two states differ only along \(\Delta_2\). This proves the displacement at all interpolating remainders. The existing whole-block guards give both halves throughout all actual twenty-step iterations. Beyond the original stop, the checked first half from each upper representative \(r_2\in\{637,639,\ldots,671\}\) takes it to the corresponding lower representative with \(J\) decreased by 36. [Prefix certificate](../box_mod12/README.md#routes-from-the-first-block-exits) · [Half-return checker](structural.py)

Use the canonical representative \(\bar r_2\in\{601,603,\ldots,635\}\) of \(J\bmod36\). The ten-step count to that section is \((J-\bar r_2)/36\). Any lift with \(36\mid\Delta J\) keeps the same canonical row. Its prefix through this section and into the oriented checkpoint is longer by \(L=10\Delta J/36\) steps, since the bounded continuation has the same chart cell and discrete labels. The oriented checkpoint on the common row is affine in \((B,p,r)\). At the oriented section use the centered variable \(Y=q-Bv/6\), where \(v=(5,4,3,2,1)\). The centered first coordinate decreases by one per update in an oriented return. The anchor at the canonical second-block endpoint decreases by \(L\). Both oriented matrices below have the same first row as that endpoint, so the bounded passage preserves this anchor difference. The shift in \(Bv/6\) adds ten to the decrease. Thus an equal-anchor return comparison must advance the source by \(L+10\) steps.

## Two oriented return slopes determine the lift

The checked canonical sections yield just two matrices of linear coefficients, in the row order \(q_1,\ldots,q_5\) and column order \(B,p,r\). Their oriented return velocities are in difference-state units per step.

\[
 M_{\rm reg}=\begin{pmatrix}
 17/18&-14/9&10/9\\ -4/9&-4/9&17/9\\ -2/3&1/3&4/3\\
 4/9&-14/9&10/9\\ -17/18&-4/9&17/9
 \end{pmatrix},\quad
 d_{\rm reg}=(-1,-11/10,-3/5,-1,-11/10),
\]
\[
 M_{\rm exc}=\begin{pmatrix}
 17/18&-14/9&10/9\\ 37/18&-31/9&8/9\\ 3&-3&-1\\
 25/9&-20/9&-14/9\\ 25/18&-10/9&-7/9
 \end{pmatrix},\quad
 d_{\rm exc}=(-1,-7/5,-3/10,1/5,1/10).
\]

The regular second-coordinate equation is
\[
 -\frac{64}{3}+\frac{17}{9}d_r
 =-\frac{11}{10}\left(10+\frac{10}{36}(60-4d_r)\right),
\]
which yields \(d_r=-12\). The exceptional second-coordinate equation gives the same answer. Substitution in all five rows gives, respectively,
\[
 M(12,18,-12)^{\mathsf T}-2v =(-40,-44,-24,-40,-44)
 \quad\text{and}\quad (-40,-56,-12,8,4).
\]
These are \(40d_{\rm reg}\) and \(40d_{\rm exc}\). The regular word lengths are 10 and 20, and the exceptional word length is 10. The original checkpoint guards provide at least ten oriented returns, including the two or four used here, and first-coordinate transience through them. The finite calculation checks the route slopes and guards at all 1,584 canonical sections from the two independent cohorts. The canonical remainders represent every original section and no base-bucket witness matching is required. [Oriented guards](../cohort1_reconstruction/README.md#checkpoint-and-matching) · [Exact slope record](structure.json)

It follows that \(\Delta J=108\), so the original representative gets \(r_2\mapsto r_2+36\bmod72\) and the canonical checkpoint takes \(L=30\) extra steps. Its centered state is the old checkpoint advanced by 40 steps. Both paired suffixes have the same autonomous dynamics from this common state, by permanent pairing and base covariance. The matched points therefore satisfy
\[
 B'=B+12,\qquad q_0'=q_0+(0,-12,18,-6,-18),\qquad T_{B'}(q_0')=T_B(q_0)-10.
\]
The period and first cyclic centered difference agree. The first word time includes the same one-update conversion in both bases. [Permanence](../../eleven_centered_clock.md#4-actual-paired-permanence-and-capture) · [Base covariance](../../odd_base_covariance.md#universal-odd-digit-base-translation)

The original checkpoint lattice has \(p=p_0+36h+2r\), where \(p_0=((r_2-\kappa)/2-B)\bmod36\). The lift raises \(p-2r\) by 42, and \(p_0' -p_0\equiv (r_2'-r_2)/2-12\equiv6\pmod{36}\). Thus \(h'=h+(42-p_0'+p_0)/36\) is an integer. The unchanged residue \(r\bmod6\) and the unchanged three incoming counts give a bijection of the original label lattices. The inverse subtracts the same input shift. In the canonical half-return convention the congruence is \(p-2r\equiv (\bar r_2-\kappa)/2-B\pmod{18}\), which keeps both relabelled branches in one lattice class.

The translation discards at most a bounded-thickness original chart collar, including the cohort-0 alignment collar. Each collar has \(O(B^4)\) difference states and conditional mass \(O(B^{-1})\). On the overlap use each arm's own ordered-word multiplicities \(11!\prod g_i+O(B^5)\). Each gap changes by \(O(1)\), so the total variation of the transported weights is \(O(B^{-1})\) after normalization. Also \(T_B=O(B)\) uniformly, so the normalized clocks and inputs differ by \(O(B^{-1})\). For every fixed bounded Lipschitz test on the original position and scaled first cyclic-word time, the difference of twelve-shift expectations is \(O(B^{-1})\), in both box-conditioned arms. [Original weights](../box_lane_table/README.md#decomposition-and-original-weights)

## Locality

The archived threefold enlargement witness at \(B=60{,}005\) is \(q=(56516,44321,34724,25107,14328)\). Its cohort is 2 and \(H=-18\), giving the inadmissible initial count \(n=(H+6)/6=-2\). The concrete sorted-order guard occurs one step into the putative \(S^5\) word. At phase 1 the comparison of digit 2 with digit 8 is \(2(b-e)-B+1=H\), since under the first \(S\) update the third coordinate is \(b-e\). The comparison is negative and the actual order is \(U\) there. This is an internal digit-order guard in the prefix. It explains the missing lift proof even where pairing itself may remain intact.

The normalized slack \(2(b-e)/B-1=-19/60005\) is strictly negative along the entire ray. At the archived point the whole-word clocks are \(6980\) and \(6974\). Scaling the input and base by 1,201 gives \(B=72{,}066{,}005\), \(H=-22{,}818\), and again a negative digit 2 versus digit 8 comparison at phase 1 in both lifts. The independently replayed clocks are \(8{,}356{,}189\) and \(9{,}000{,}522\). This rules out interpreting the locality failure only as a finite-base guard offset. [Exact values and reconstruction](locality.json)

<a id="the-finite-route-criterion-for-the-common-component"></a>

## Finite-route convergence

The guarded routes below partition each canonical section by congruence conditions. A leaf with a nonempty pullback at a given base spans the whole spatial box up to the fixed-thickness alignment collar. Its fixed-base lattice index gives its limiting mass. The normalized word time equals one affine graph uniformly to \(O(B^{-1})\) on the leaf. The original multiplicity asymptotic \(11!\prod_i g_i(q_0)+O(B^5)\) therefore gives a weighted Riemann sum for each graph and its spatial density \(\nu_C\). [Original weights and chronology](../box_lane_table/README.md#decomposition-and-original-weights)

On each class of a common multiple of 12 and the pulled-back base steps, the leaf masses are fixed and the joint law converges. The twelve-shift estimate identifies the limits of any two successive fine classes separated by 12. Iteration identifies all fine-class limits in each arm, and finite interleaving gives the limit along the entire class mod 12. This argument also applies to sublaws indexed by unions of canonical labels. For a given canonical label, both original remainders \(r_2\) and \(r_2+36\) enter together. [Base period collapse](../../odd_base_covariance.md#what-covariance-doesand-does-notimply-for-a-probability-law)

### Coverage algorithm and clock bound

`sections.prefixes()` reconstructs \(q_2\), its affine \(J\), and the paired block from the accepted switch records for cohorts 1 and 2. `structural.sections()` sends the upper representative to the canonical one by its guarded half-return, and returns the label, chart bounds, affine oriented checkpoint, paired word, integer displacement, and elapsed time from the canonical section. These are the inputs to `coverage.intrinsic()` and `coverage.groups()`. With \(X=6q-Bv\), take \(A=X_1\) and the first \(D=X_j\) for which the \((p,r)\) minor is invertible. Solve for \(p,r\) and rewrite all five rows in \((B,A,D)\). Group rows with the same pivot index, five rewritten rows, leading pivot maps, paired word, and centered displacement. `physical_lattice()` intersects the lattice for \(B\equiv5\pmod6\) with the five congruences \(X_i\equiv5i\pmod6\). Some physical starts may lie outside the original source lattices. Permanent pairing for the actual sources comes from the accepted box route, and the continuation uses the autonomous paired map. `CheckpointGuards` tests the leading part of each inequality at the box extrema and takes the smallest intercept across the source offsets. A positive leading minimum determines a base cutoff, and an inequality with zero leading minimum requires a nonnegative minimum intercept. [Prefix reader](sections.py#L11) · [Structural sections](structural.py#L14) · [Pivot selection](coverage.py#L24) · [Group assembly](coverage.py#L128) · [Physical lattice](coverage.py#L111) · [Guard bounds](coverage.py#L49)

The continuation is the accepted Task 066 engine. The independent reconstruction reproduced its archived output. A separate simulator tested each primary residue class in eleven complete runs below the proof cutoff, as well as 720 random points above it. [Task 081 independent reconstruction](independent-review/RECONSTRUCTION.md#3-coverage)

- `verify.walk_macro()` takes the five affine rows, the clock initialized to \(-q_{4,1}\), the paired word and centered drift, the current affine lattice, the guard interface, and the initially empty return history. `repeat_guards()` verifies the translation and provides each affine phase comparison \(G+kD\prime\) under repeated returns. The engine chooses the smallest candidate exit quotient at a lattice point obtained by rounding the box midpoint in the lattice basis. `Lattice.residue()` and `Lattice.split()` make the fractional part of \(G/(-D\prime)\) constant on each child coset. On each child the engine checks every guard at repetitions 0 and \(\lfloor G/(-D\prime)\rfloor\), updates the state and clock by \(\lfloor G/(-D\prime)\rfloor+1\) repetitions, and then checks proposed single updates. It searches for two consecutive copies of a word of length at most 100 with the same nonzero constant displacement. It continues with that word when the trial point gives more than 100 repetitions and rejects a switch that finds no continuation within 400 updates. Repetition counts and guards are exact, and each split is arithmetic. [Continuation](../box_lane_table/verify.py#L372) · [Repeat guards](../box_lane_table/verify.py#L308) · [Lattice operations](../box_lane_table/verify.py#L151)
- `verify.terminal_time_bounds()` handles the case in which the state at the real root of the declining guard is constant. It certifies the truncated return \(m_*=\lfloor G/(-D\prime)\rfloor-k\) for a fixed offset \(k\) and requires \(G/(-D\prime)-k\ge5042\). It enumerates all physical integer endpoints in the resulting bounded line segment and replays the centered clock to obtain each endpoint's delay to its first cyclic difference. Each leaf returned by `walk_macro()` contains the affine clock at the real truncation, the return history with the truncation offset, the terminal replays, and the resulting coset in \((B,A,D)\). A failed guard or missing continuation rejects the certificate. [Terminal bounds](../box_lane_table/verify.py#L328) · [Centered replay](../channel_continuation/certify.py#L292)

The exact pass closes 99 checkpoint families in 6,555 intrinsic leaves, with at most three repeated words per leaf in the continuation. The cutoff \(B\ge86{,}496{,}843\) certifies the guarded continuation on the closed cohort-1 and cohort-2 charts. In cohort 0 the accepted five updates after the initial exit take the state to the exit of the same-base cohort-1 input \(q_0+\alpha_0\), where \(\alpha_0=(-5,0,1,2,-1)\). The two inputs have the same \(n\), and their trajectories agree from that exit with \(T_0(q_0)=T_1(q_0+\alpha_0)+5\). This applies wherever both inputs lie in the box and keeps \((z,r_1,r_2)\) in the cohort-1 convention. The excluded cohort-0 inputs lie in an original chart collar of bounded lattice thickness and conditional mass \(O(B^{-1})\). The cutoff for cohort 0 is used on this overlap via the cohort-1 guards. The chart translation fixes the leading values of \(BP(x)\) and \(BR(x)\). Its \(-5\) change in \(a\) cancels the five extra updates in the graph error. [Full coverage record](coverage.json) · [Same-base alignment](../box_mod12/README.md#original-weights-and-the-effect-of-cohort-0)

The elapsed prefix is exactly \(a-q_{4,1}\) in the two guarded cohorts. The companion checks the anchor decrements through the short switch and the passage from the canonical section to the checkpoint. The integer truncation is earlier than its real affine clock by less than the final return length \(\ell\). A point 5,041 returns before the integer endpoint and the endpoint differ by 5,041 times a nonzero integer drift. They cannot both lie on a centered cycle, whose coordinates lie in \([-2520,2520]\). Thus the first cyclic difference is later than the earlier point. If \(d_*\) is the greatest delay among the physical endpoint replays for that leaf, the error \(T_B-(a+E)\) from the returned affine clock \(E\) lies in \([-5042\ell,\,1+d_*]\). The extra update converts a cyclic difference time to a cyclic word time.

The uniform terminal window over all physical leaves is \([-201{,}680,167]\) updates. After substituting each source chart and including its affine intercept, every covered input satisfies \(-201{,}734\le T_B-B\tau_s(x)\le136\), giving the uniform absolute bound \(201{,}734\). This bound includes cohort 0 on the alignment overlap. The exact rational bounds by cohort are in the companion record. [Clock bounds](limit_weights.json)

### Original-lattice pullback

An independent count starts with the determinant of the inverse original prefix map. At fixed \(B,c_0,z,r_1\) its columns in \((p,r,N,n,m)\) have determinant of absolute value 132. For example the unimodular coordinates \(q_{0,3},q_{0,2},q_{0,1}+q_{0,3},q_{0,3}+q_{0,4},q_{0,2}-q_{0,5}\) reduce that determinant to the determinant of the last three coordinates of the \(N,n,m\) columns, which are \((22,0,0),(8,5,3),(10,-2,0)\). The independent \(r_2\)-congruence splits each such coset into 36 equal-index cosets, so the total index is \(132\cdot36=4{,}752\). The further \(r\bmod6\) refinement gives \(28{,}512\). A separate set check expands the actual guarded canonical labels, includes the upper half-return labels and the aligned copy in cohort 0, and verifies that each original label occurs once. [Original inverse map](../box_lane_table/README.md#decomposition-and-original-weights) · [Independent census check](limit_weights.py)

For a canonical \(\bar r_2\) and its switch constant \(\kappa\), the parent coset in \((B,p,r)\) has point \((5,(\bar r_2-\kappa)/2-5,0)\) and columns \((6,-6,0),(0,18,0),(0,2,1)\). For a checkpoint leaf \(y+L\mathbb Z^3\) in \((B,A,D)\) with column matrix \(L\), the three Cramer numerators for \(L^{-1}((B,A,D)-y)\) must vanish modulo \(|\det L|\). Substitute each source's \(A,D\) and clear denominators. Intersect the resulting integer congruences one at a time with the parent. For a congruence \(u\cdot z+c\equiv0\pmod m\) on \(z=z_0+\sum_i k_i v_i\), existence requires \(\gcd(m,u\cdot v_1,u\cdot v_2,u\cdot v_3)\) to divide \(u\cdot z_0+c\). Integer column operations combine the three increments, after which solving a single congruence gives the new point and the kernel columns.

Let \(z_*+M\mathbb Z^3\) be a nonempty pullback, \(D_*=|\det M|\), and \(g_*\) the gcd of the base components of its columns. It is present exactly when \(B\equiv(z_*)_B\pmod{g_*}\). At such a base its \((p,r)\) fibre has covolume \(D_*/g_*\). The incoming counts \(N,n,m\) remain free in the original lattice parametrization. Hence the index in the original five coordinates is \(132D_*/g_*\) and its limiting conditional mass is \(g_* /(132D_*)\) for that source. The canonical parent has fixed-base covolume 18 and comprises twelve primary classes, two full remainders times six \(r\)-residues. The sum of compatible leaf masses is checked against these twelve units at every base residue for every source.

The exact pass finds 99,076 nonempty source-leaf intersections. Every nonempty pullback has base step 12 and belongs to exactly one of the two arms. The companion records the base steps and fixed-base covolumes, and its `--leaves` option emits each compatible base congruence and pullback determinant. [Exact pullbacks](limit_weights.json) · [Pullback calculation](limit_weights.py)

### The two limit laws

In the chart \(x=(a,b,c,2(b-e)-B,B-c-d)/B\), the box is
\[
 C=10^{-5}([94040,94130]\times[73770,73850]\times[57940,58015]
             \times[20,80]\times[200,300]).
\]
The leading return count and late-return coordinates are
\[
 A_N(x)=\frac{x_1+x_3-3/2-(11/2)x_4-5x_5}{22},\qquad
 P(x)=x_3-\tfrac12x_4-2A_N(x),\qquad
 R(x)=x_2-\tfrac32x_4-x_5-4A_N(x).
\]
For the spatial law recover the original normalized differences
\[
 y=(x_1,x_2,x_3,1-x_3-x_5,x_2-(1+x_4)/2),
 \quad (g_0,\ldots,g_5)=(1-y_1,y_1-y_2,y_2-y_3,y_3-y_4,y_4-y_5,y_5).
\]
Then \(d\nu_C(x)=w(x)\,dx/\int_Cw(u)\,du\) on \(C\), with \(w(x)=11!\prod_{i=0}^5g_i\). The chart has constant Jacobian, which cancels in this normalization.

For \(\beta=5,11\), the limit is
\[
 \mu_\beta(dx,dt)=\nu_C(dx)\sum_s\frac{w_{\beta,s}}{28{,}512}\,\delta_{\tau_s(x)}(dt),
 \qquad \tau_s(x)=x_1+s_0+s_pP(x)+s_rR(x).
\]
Each table entry is in units of \(1/28{,}512\).

| Lane | \((s_0,s_p,s_r)\) | \(B\equiv5\pmod{12}\) | \(B\equiv11\pmod{12}\) | 5 minus 11 |
|---|---|---:|---:|---:|
| S | \((-11/324, -133/162, -35/81)\) | 7202 | 7202 | 0 |
| S | \((-7/408, -57/68, -15/34)\) | 2940 | 2961 | -21 |
| S | \((1/576, -247/288, -65/144)\) | 4900 | 4935 | -35 |
| M1 | \((31/96, -19/16, -5/8)\) | 1328 | 1328 | 0 |
| M2 | \((29/54, -38/27, -20/27)\) | 2688 | 2700 | -12 |
| L | \((139/144, -133/72, -35/36)\) | 8741 | 26032/3 | 191/3 |
| L | \((863/918, -836/459, -440/459)\) | 713 | 2126/3 | 13/3 |
| **Total** | | **28512** | **28512** | **0** |

The row differences agree with Task 071, including positive and negative totals \(+68\) and \(-68\). The frozen test is 1 on both long graphs and 0 on every other graph, as verified by exact vertex bounds. This gives \(\mu_5(F)-\mu_{11}(F)=68/28{,}512=17/7128\) directly from the full table. [Frozen test](../box_lane_table/frozen-test.md) · [Task 071 signed profile audit](../../../docs/drafts/research-direction/replay-task071/REPLAY.md#2-claim-ii-signed-profile-weights-68-maximality-fibre-fraction)

### Corollary on the common sublaws

The residual label set consists of 176 original labels in each arm and synchronized bucket, and is invariant under \(r_2\mapsto r_2+36\bmod72\). Its complement defines a sublaw of mass \(26/27\) in each arm. The same finite closure and twelve-shift argument give limits for both common sublaws. Lemma O matches their labels across the six-shift and thus identifies the two common sublimits. The full-box convergence and separation follow from the guarded closure, the twelve-shift, and the full table above, independently of Lemma O. [Residual labels and six-shift matching](../box_mod12/README.md#signed-whole-box-residual)

## Repeated pairing failures

The accepted prefix has exactly \(N\) repetitions of one late word. On each repetition there is one failure at phase zero and, for \(z=1\), a second at phase ten. For \(z=2\) the phase-ten equality sorts while paired. The remaining failures are in the switch of at most 21 updates. Thus the exact failure count is \((1+\mathbf1_{z=1})N\) plus the listed switch count. There are no failures once the certified paired block begins. The coefficient
\[
 N/B=A_N(x)+O(B^{-1}),\qquad
 A_N(x)=\frac{x_1+x_3-3/2-(11/2)x_4-5x_5}{22}
\]
has minimum \(1/55000\) on the closed box. In particular any fixed cutoff on the *number* of failures eventually excludes the entire box-conditioned population. The compressed history accounts for every failure, so its unhandled-failure mass is zero at all sufficiently large bases. [All-cohort failure words](../box_mod12/README.md#routes-from-the-first-block-exits)

If the observed failure mark is the integer list of failure times through the last failure, the twelve-shift preserves that list pointwise on the matched boxes. The preservation follows from the unchanged \((c_0,z,r_1,n,m,N)\) and the same switch. It holds for every fixed bounded test Lipschitz in \((x,T_B/B)\) uniformly in the failure list, with \(O(B^{-1})\) discrepancy. In particular no finite-modulus limit family for those tests can distinguish a finer residue within either class mod 12. Separately, the scaled failure point measure \(B^{-1}\sum_{t\,\mathrm{fails}}\delta_{t/B}\) converges jointly with the normalized input through all \(B\equiv5\pmod6\), when tested against bounded Lipschitz functions of time. Given \(x,z\), its leading measure is Lebesgue measure on \([\tau_F(x),\tau_F(x)+20A_N(x)]\) multiplied by \((1+\mathbf1_{z=1})/20\), where \(\tau_F(x)=5x_5+5x_4\). The finitely many switch failures vanish after division by \(B\). The two overshoots each have limiting input weight \(1/2\) with density \(\nu_C\). This leading failure measure follows by a Riemann sum over repetitions in each overshoot. It is a statement about failure times and spatial input before adding the clock.

There is additional arithmetic if a **post-failure checkpoint mark** is retained. On the fixed repeated-failure cohort \((c_0,z,r_1)=(0,1,17)\), take \(r\equiv2\pmod6\) and \(r_2=603\). This primary class has mass \(1/28512+O(B^{-1})\) under whole-box conditioning. In the arm \(B\equiv5\pmod{12}\), its normalized clock is on M2 for \(B\equiv5\pmod{24}\) and on L for \(B\equiv17\pmod{24}\). In the other arm it is on L for \(B\equiv11\pmod{24}\) and on S for \(B\equiv23\pmod{24}\). The long clock keeps both refined profiles with fibre weights \(2/3\) and \(1/3\). The fixed bounded Lipschitz Task 066 test, integrated against the marked class's sublaw, therefore alternates between 0 and \(1/28512\) in each arm. The partner \(r_2=639\) reverses these phases and restores the sum. The diagnostic isolates arithmetic on a positive-weight repeated-failure channel, with the alternating discrimination supplied by preserving a mark *after* the failures. It gives no such discrimination for the same two remainders pooled together or for the unmarked box law. Eight fresh whole-word clock controls at 90 million passed the frozen plateau predictions. [Exact marked weights](marked-weights.json) · [Frozen controls](freeze-marked-clock.md) · [Control clocks](marked-clock.json)

## Revisions (19 Nov, Task 089)

- Specified the guarded continuation's inputs, leaf outputs, and source functions and identified the scope of its independent reproduction.
- Recorded the uniform terminal and full graph clock bounds.
- Made the scope and collar of the cohort-0 alignment explicit.
- Counted the original labels by the inverse-prefix lattice index and verified their bijection with the guarded sections.
- Stated the Lemma O comparison as a corollary on the common sublaws.
- Computed each canonical source's leaf pullbacks, fixed-base indices, and base congruences.
- Tabulated all seven graph coefficients in both arms and checked the signed differences against Task 071.

## Reproduction

Run from the repository root. The snapshot wrapper copies all executable inputs into a private directory and takes the exclusive compute lock by default. The local compact checks were bounded to two minutes per command. The marked-weight reader reads its fixed archived input directly. The original SHA-256 manifest pins Task 081's outputs and two prediction files.

```sh
bundle=kaprekar-investigation/eleven_access/box_structure
python3.13 -B "$bundle/snapshot.py" structural.py | cmp - "$bundle/structure.json"
python3.13 -B "$bundle/snapshot.py" coverage.py | cmp - "$bundle/coverage.json"
python3.13 -B "$bundle/snapshot.py" locality.py | cmp - "$bundle/locality.json"
python3.13 -B "$bundle/snapshot.py" marked_clock.py | cmp - "$bundle/marked-clock.json"
python3.13 -B "$bundle/marked_weights.py" | cmp - "$bundle/marked-weights.json"
```

The Task 089 companion reuses the guarded continuation inside the same locked snapshot.
The optional `--leaves` flag includes the per-source leaf ledger in the output.

```sh
python3.13 -B "$bundle/snapshot.py" limit_weights.py "$PWD/$bundle/coverage.json" |
  cmp - "$bundle/limit_weights.json"
```
