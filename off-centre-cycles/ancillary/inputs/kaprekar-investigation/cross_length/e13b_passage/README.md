# Thirteen-digit passage from chamber exit to the marked outer section

A full-dimensional original source cell in class 0 reaches the normalized gate after first chamber exit. Its first \(R\le2000\) state is a previously certified key leading to least period 177. The cell's exact limiting original weight is
\[
\mu=\frac{3773279180259731755300995131345652558791147811352345556322576361}{3392252763571048539774051887827968000000000000000000000000000000000000000000000000000}
\quad(\text{approximately }1.11232\cdot10^{-21}).
\]
It is disjoint from Task 158's source cell. An independent two-word certificate supplies positive original weight to the normalized gate in each class. The proof currently certifies no positive class-6 limiting weight through to a mark. With Task 158, \(L_0=\eta+\mu\) and \(L_6=0\), where \(\eta=3575/1350641342287205909485672071168\). These lower bounds refer to denominator \(B^{13}-B\); the missing weights are \(1-L_0\) and \(1-L_6\). The any-radius contrast certified by these bounds is \(L_0+L_6-1\). The ten-point target remains open. [Task 158](../tl3_survival/README.md) · [Class-0 certificate](marked_cell.py)

## Premise

The strict phase-zero, centre-pair chamber cell has a primitive wall \(L_3\), with a rank-five kernel in the six free integer gaps. At \(B=114736\) the literal slide first fails at time 126 with \(J=21\), \(\eta_{\rm exit}=5\), and \(L_3y=3=6+5-8\). The ray direction has strict positive leading slacks with the duplicated centre rows uniquely minimal. The exact inverse \(n=y-Jv\) recovers the inlet. Three Task 163 outer keys independently replay to periods 177, 1 and 33 on their archived bases and on the least eligible bases above 32774 on their cosets. The class-0 period-1 key shares the anchor sheet of the favorable key, so the label must still be read. [Literal premise checks](premise.txt) · [Task 163](../e13_selection/README.md)

## Return and weight calculation

Use \(X_i=7q_i-(7-i)B\) for \(i=1,\ldots,6\) and \(G=(-X_1,X_1-X_2,\ldots,X_5-X_6,X_6)\). A stopped history carries its original variable \(q_0\), its arithmetic lane \(q_0\in p_B+\Lambda\), and its exact affine state and time. Write its state as \(X=Mq_0+B b+c\). The directions are \(M\Lambda\). In a paired physical phase each within-level and adjacent cross-level comparison is affine. For a general subtraction branch, fix the borrow bits and the complete digit order with a deterministic priority on ties. The raw digits then have the form \(\pm q_i+B b_j+c_j\), including the central digit, so the state update and its guards are affine as well. A bounded switch multiplies the chart by the branch derivative and carries the exact time and original lattice.

A return word of \(p\) phases acts on an admitted chart by translation \(d\), with its composite linear part fixing \(M\Lambda\) and fixing \(d\). The offset also has to agree over the base progression. At phase \(e\) of copy \(j\) each guard has slack \(h+j\dot h\), with \(h=aX+B b_0+c_0\) and exact ties resolved in the phase guard. For a row losing \(s=-\dot h>0\), the first failing copy is \(1+\lfloor h/s\rfloor\). Compare these integers and read the earliest failing phase of the winning copy. A residue \(r=h\bmod s\) gives the affine count \(1+(h-r)/s\), updating \(X\mapsto X+d[1+(h-r)/s]\) and time by \(p[1+(h-r)/s]\) before the bounded phase stem. The new direction map is \(I+d a/s\). Split the incoming lattice, before applying that map, by each feasible residue. If \(d_0=\gcd(s,aM\Lambda)\), the realized residues in the original lane form a coset modulo \(d_0\), and each further lane has index \(s/d_0\). At tied winning ratios compare their derivatives on the incoming tangent and their leading base terms. When these coincide the offset decides on the lane. A genuine leading equality between distinct forms pulls back to an original hyperplane. A numerical tie in the sample allows either row only on the appropriately guarded branch of the race.

To stop at an actual section, for each phase in the repeated word intersect the integer copy interval where all seven gaps meet the threshold with the interval where the prefix to that phase is valid. For the normalized threshold use \(64G_i+B\ge0\). For the fixed outer threshold use \(G_i+2000\ge0\). A positive gap drift gives a ceiling lower bound and a negative drift gives a floor upper bound. A zero drift either leaves the interval unchanged or makes it empty. Take the earliest feasible phase and copy, and split the active lower bound by its remainder. For a rising active gap with row \(a\) and drift \(s>0\) this gives the direction correction \(I-d a/s\), followed by the phase derivative, or equivalently the correction in the phase coordinates. At fixed base the stopped state can have rank zero while still depending on \(B\). This distinction matters at the normalized stop. A fixed finite suffix can use rank zero for a fixed key only after its base dependence also vanishes.

For disjoint original cells and one fixed compatible base refinement, a residue-refined route of original lattice index \(I\) on a guarded open polytope \(P\) contributes \(I^{-1}\int_P13!\prod_{i=0}^6u_i\,du_1\cdots du_6\). Here \(u_i=g_i^{\rm original}/B\). At a fixed base one can instead sum exact original word multiplicities over the admitted lattice points. Pool the weights of histories that merge physically while keeping their source sets separate. Use the original weights through every projection. [Task 136](../../eleven_access/r2_routes/README.md) · [Task 137](../../eleven_access/z1_stages/README.md) · [Word multiplicities](../../odd_gap_transport.md)

For any one fixed finite affine route, if the map to the fixed outer section has direction rank \(r>0\), a finite list of exact keys absorbs at most \(O(B^{-r})\) original weight through that route. Each key imposes \(r\) independent equations in the six source coordinates. At most \(O(B^{6-r})\) original gap states satisfy them, and any one state has \(O(B^7)\) original words. A finite union of these routes has the same vanishing conclusion with the minimum positive rank. A population lift through such charts therefore needs a further rank collapse, a number of split charts growing with base and an accompanying original-weight estimate, or a continuation certificate that handles the variable state at another section. This bound concerns fixed routes. The retrospective rank counts supply no limiting upper bound on any mark.

## A full-dimensional marked route

Use the physical chamber-entry state from archived draw 502 at \(B_0=140000\) as a direct original source. Its centered differences are \(S=(90594,64043,49868,17220,10073,13923)\). Define
\[
n=63\left\lfloor\frac{B-B_0}{63 B_0}\right\rfloor,\quad
X=(1+n)S+\sum_{i=1}^6 z_i D_i,\qquad
|z_i|\le\left\lfloor n a_i/12\right\rfloor,
\]
where the columns \(D_i\) are in the verifier and
\[
(a_1,\ldots,a_6)=(9584/7,169/21,902/21,1192/21,826/9,2063/63).
\]

Every base \(B\equiv0\pmod7\) with \(B\ge8960000\) is eligible. There are six returns, of respective lengths \(6,12,30,18,36,18\). The bounded bridges after the first five returns have lengths \(0,25,2,100,56\). The six copy counts are the sample vector \((1370,8,42,56,89,6)\) plus the independent count increments \(n a_i+z_i\). The last repeat then uses a 16-step remainder phase. The verifier reconstructs every mask by literal subtraction at the seed, checks the chosen exact wall remainder and the six slopes, and checks each physical inequality for all parameters at both ends of every repeat. The source lies at strict leading chamber entry at time zero. Along all-left motion \(X_1\) decreases by 7 per step, and at first exit it is still \(>B/64\). The computations explicitly certify every step before the first fixed-outer arrival outside \(R\le2000\).

The first five quotients lower the direction rank \(6,5,4,3,2,1\). At the first fixed-outer crossing the last free direction is parallel to the return translation. Its exact quotient cancels all six parameter directions including the common \(nS\) direction. Thus the whole source reaches the exact key \(G=(-1988,1484,1491,-1974,315,322,350)\). It follows the archived suffix with preperiod 437, least period 177, and peak \(R=1988\). The normalized gate has necessarily been reached by this time. The proof charges every unresolved path outside the source box.

The six source columns are independent because their induced increments in the six consecutive exact counts are the six coordinate vectors. They are integer multiples of 7 and \(n a_i\) are integers, so every sufficiently large base in the class is compatible. The original gap \(g_0\) satisfies \(g_0<B/8\), whereas Task 158's original source has \(g_0\ge B/8\). For \(w_i=a_i/12\) and \(v_i=z_i/n\), the exact limiting weight in this box is
\[
\frac{13!}{B_0^6}\int_{\prod_i[-w_i,w_i]}
\prod_{j=0}^6\frac{B_0+G_j(S+\sum_i v_i D_i)}{7B_0}\,dv_1\cdots dv_6.
\]
The direct count uses one original word state per integer sextuple \(z\) and \(n/B\to1/B_0\). The distinct-digit term leaves a uniform \(O(B^6)\) correction per state for inputs with at least one tied digit. Polynomial integration, using vanishing odd moments and even mean moments \(w^k/(k+1)\), gives \(\mu\) above. [Exact route and weight log](marked-cell.txt)

## A positive gate-only route in both classes

The two seed exits have \(X_3+X_4=-14\). The first word is the chamber permutation, followed through phase 5. The next 12-step word has two occurrences of mask 12 with five intervening zero masks and displacement \((-84,-70,-7,7,-28,-14)\). Its phase-5 exit coordinates vary in five independent integer directions on the wall, with an independent earlier chamber count. Paired physical guards hold throughout both words. Inside the second word \(X_1\) decreases by 7 every update and the other gaps are above the normalized threshold. The exact first hit occurs at post-exit time \(\lceil(64X_1^{\rm exit}-B)/448\rceil\). The verifier checks all inequalities on the product box at both ends of each repeat and tests the exact gate time by literal subtraction in each class.

For every \(B\ge4000000\) in either class the gate-only certified original weight is at least \(34749/1443710306877440000000000000000000000000000000000000000000\). This follows by counting the independent integer variables and bounding each initial physical gap below by \(B/500\). At fixed section residue and phase this route still has four tangent dimensions. The gate-only weight supplies zero credit to \(L_0,L_6\) without a marked continuation. [Verifier](gate_cell.py) · [Exact gate log](gate-cell.txt)

## Archived passage structure

The physical replay uses all six frozen Task 163 panels of 512 original draws. It derives the first chamber exit again, checks the exact first gate and fixed outer times, then fits paired returns of at most 64 steps from three observed copies. A return must have a nonzero translation, fix that translation, and fix the rational span of the inherited tangent. The minimum accepted total run length is \(\lceil\sqrt{B}/8\rceil\). The finite-base checker applies the affine rank update at the first binding guard or, when the section truncates a word, at its active first-hitting threshold. Primitive steps between returns are actual full-digit subtraction. The starting rank is that of the full first-exit sheet. Subsequent ranks describe selected affine charts before extra restrictions from narrower original domains and residues. The return search and observed frequencies supply no deterministic population weight.

Each line uses 512 original draws. The six rank entries count the first fixed outer section ranks \(0\) through \(5\). The two early hits occurred at or before the chamber exit. The 32 other omitted draws lacked a leading entry or six physical all-left phases at that entry. Both amounts are still in the original denominator.

| Base | Read | Favorable mark among read | Outer ranks 0 / 1 / 2 / 3 / 4 / 5 | Median / max returns to gate | Longest bridge to gate |
|---:|---:|---:|---|---:|---:|
| 140000 | 500 | 494 | 1 / 26 / 355 / 102 / 16 / 0 | 4 / 9 | 9010 |
| 140006 | 508 | 508 | 0 / 46 / 332 / 107 / 21 / 2 | 4 / 10 | 8480 |
| 280007 | 505 | 503 | 0 / 62 / 346 / 84 / 13 / 0 | 4 / 8 | 19752 |
| 280013 | 506 | 506 | 0 / 69 / 350 / 72 / 12 / 3 | 4 / 9 | 19186 |
| 1120000 | 510 | 507 | 0 / 88 / 365 / 48 / 9 / 0 | 4 / 10 | 64933 |
| 1120006 | 509 | 509 | 0 / 98 / 361 / 42 / 8 / 0 | 4 / 10 | 76537 |

On the 3038 read paths the 13094 recognized returns cover exactly \(51790731/64814714\) of exit-to-gate updates. Of the complete returns, 12390 have two or more guard rows failing in the same copy. Of these, 12389 contain an equal tangent ratio and 1047 contain a distinct-ratio tie at integer precision. These two descriptions can overlap in one return. Every recognized return comes within 196 raw units of at least one reflected pair in its first word, and 1360 approach all three, possibly at different phases. These windows use \(|X_i+X_{7-i}|\) for \(i=1,2,3\). The reader assigns 8596 direction-rank drops to maximal guard returns and 231 to normalized crossings inside a return, with zero primitive drops. The word cap is a diagnostic limit. A numerically recurring word can also fail the inherited-chart identity test, as it did at least once on 1157 of these paths.

The gate-to-outer paths make 15 outward crossings of the normalized threshold. They spend \(618098/3006775\) of their steps in 1379 recognized returns. All distinct outer keys and their terminal recrossings remain those certified by Task 163. The sole chart of rank zero at first outer arrival is the class-0 seed used in the exact marked certificate. Draw 221 at \(B=1120006\) takes 111 unpaired steps before its first normalized gate, with maximum unpaired \(R=390669\). The other read paths remain paired on that passage. This exact exception rules out an all-paired replay even on the archived population.

The observed uncompressed growth lies in the paired bridges before the normalized gate, with repeated near contact of the second and third reflected pairs. At the largest bases, on the longest bridges, the respective counts in the second and third comparison windows are \(18942,21766\) of 64933 steps in class 0 and \(25217,26106\) of 76537 in class 6. This observation identifies a place to improve the reader. It gives no uniform bound on the unresolved weight in those bridges, on rank-positive outer images or on energy stalls. [Exact per-panel ledger](stage-summary.txt) · [Per-draw replay](replay.tsv.gz)

The earlier frame passes are [Task 093](../../../docs/drafts/research-direction/o2l-task093.md), [TL selection](../../../docs/drafts/research-direction/o2l-tl-selection.md) and [U-1 access](../../../docs/drafts/research-direction/o2l-u1-access.md).

## Reproduction

Run from this directory. The independent premise and exact box checks are short. The archived-panel replay uses a private build and the shared compute lock. Set `GMP_PREFIX` if GMP is installed outside your compiler's default paths.

An initial diagnostic exceeded its short-check estimate and was interrupted. The complete replay was paused and then resumed with the shared lock. The reproducible runner holds the lock for the entire panel calculation.

```sh
python3 -B premise.py ../e13_selection
python3 -B gate_cell.py
python3 -B marked_cell.py ../e13_selection
GMP_PREFIX="$(brew --prefix gmp)" python3 -B run.py | gzip -n > /tmp/e13b-replay.tsv.gz
python3 -B summarize.py /tmp/e13b-replay.tsv.gz
```
