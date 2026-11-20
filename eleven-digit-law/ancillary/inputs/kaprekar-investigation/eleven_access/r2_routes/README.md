# Finite physical routes for the eleven-digit entrant law

Let $y$ be the original normalized differences, $D$ their finite leading entrance event into strict $\Pi$, and $X=6q-B(5,4,3,2,1)$ the centered physical differences. The target observable is $(y,T_B/B,Z_B)$ on $D$, with $T_B$ the first cyclic **word** time and $Z_B$ its normalized differences. The finite-route transfer below covers strict paired, tied and unpaired records. [Truncation theorem](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md)

## The original lattice calculation

**Dimension-independent.** On a fixed original open polyhedral cell $P$, retain its source variable $q_0$, its limiting density $w(y)$, and its original congruence coset $q_0\in p_B+\Lambda$. On a refined base class an admitted history at its current section has an exact representation

$$
X=Mq_0+Bb+c,\qquad t=uq_0+Bv_0+e,
$$

for fixed rational coefficients, with $X,t$ integral on the physical coset. Its leading state and time are $x(y)=My+b$ and $S(y)=uy+v_0$. The image lattice has directions $M\Lambda$. This is the affine image of the entire coset, before restricting to $P$ and its physical guards. Bounded switches multiply the chart by the actual affine branch and add their integer cost. They may collapse directions or identify different incoming points.

**Dimension-independent.** A certified repetition is a word of $p$ primitive affine phases. Its composite $R_B(X)=A X+Bf+g$ fixes the leading chart pointwise, so $(A-I)M=0$ and $(A-I)b+f=0$. Its drift on the sheet is $d=(A-I)c+g$, a fixed integer vector, since the sheet consists of integral physical states and $R_B$ is integral. Require $d\ne0$ and $A d=d$ to use translation along all repeats. The latter condition also covers $d$ transverse to the sheet. Pull every phase guard back to the start of the word, with equality resolved by the physical convention. For a row $a(X)=a_0X+B a_B+a_c\ge0$, the guard at repeat $j$ is $a(X)+j a_0d$. Nondeclining guards stay valid. Among declining guards choose a row with $s=-a_0d>0$ minimizing the exact ratio $a(X)/s$. The first repeat on which the recorded word fails is

$$
n=1+\lfloor a(X)/s\rfloor,\qquad
X'=X+n d,\qquad t'=t+p n.
$$

These are the section state and elapsed time after all complete repetitions. For a mark at the first failing phase within the next repeat, read the phase guards in chronological order from $X'$, retain the valid intervening primitive maps and add their actual costs. If several rows first fail on that repeat, any minimizing row gives the correct $n$. Identifying the first failing phase still requires the chronological reader. Changes of tied sort labels can leave the numerical update unchanged. A pairing loss is read from its own physical guards.

First compare leading slopes of the guard ratios on the source cell. A nonidentical comparison changes sign only across its original affine hyperplane. When two leading ratios agree everywhere on $P$, their cleared difference annihilates $M\Lambda$ and its leading base term cancels. Compare their exact offsets on that history. In particular, a leading tie on a lower-dimensional **image** can carry positive original mass. Its exact offset is constant on the inherited coset and may vary between histories or base residues. A leading tie with nonzero variation on $M\Lambda$ pulls back to a proper original hyperplane and has original limiting mass zero. The same argument governs individual sorting, height and pairing comparisons and the pointwise return test on a leading cell. Positive leading exit distance $D_a(y)=a_0x(y)+a_B>0$ gives a macroscopic stage. Zero distance throughout the chart gives bounded work for that candidate.

**Dimension-independent.** For the winning guard put

$$
h=\gcd(s,a_0M\Lambda),\qquad j=a(X)\bmod s.
$$

Here the gcd uses the values on a generating set for the direction lattice. The realizable $j$ comprise one coset modulo $h$. Their preimages are cosets of index $s/h$ in $\Lambda$. For a lane $j$, the exact affine map uses $n=1+(a(X)-j)/s$. It changes the leading marks to

$$
x'=x+\frac{D_a}{s}d,\qquad
S'=S+\frac{pD_a}{s}.
$$

The new lattice is the exact image of the residue sublattice under the affine exit map with derivative $I+d a_0/s$. A later residue test pulls back equally indexed cosets to the old lattice because the homomorphism onto its image is surjective. This applies to rank-dropping and rank-preserving returns and to all bounded switches. Pool histories with the same current state and continuation by retaining their separate original densities, accumulated elapsed times and admitted labels in the pooled record.

Clearing denominators in the chart gives an integer congruence on $(B,q_0)$ for each lane refinement. A fixed finite record therefore uses a finite base refinement. On a compatible base class the original index $I=[\mathbb Z^5:\Lambda]$ and the limiting affine maps are fixed. Incompatible bases give weight zero. For a compatible closed record and any continuous test $\varphi$, its limiting weighted pushforward is

$$
\frac{1}{I}\int_P w(y)\,\varphi\bigl(y,S(y),Z(y),\mathcal M(y)\bigr)\,dy,
$$

where $\mathcal M$ denotes any additional admitted marks. For a stopped record use the stopped mark $x$ and the physical elapsed time. An exact microscopic mark needs a lane-resolved readout whose limit has been verified if the mark is included in a limiting statement. The integral follows from weighted Riemann sums on the original coset. The equal-index weight $1/I$ is exact in the limit. At finite $B$ sum the actual word multiplicities of the original points satisfying the guards and congruences, using the exact integer maps above for the marks. The bounded affine offsets affect membership near proper original walls only in collars of width $O(B^{-1})$ and vanishing original weight. The collar estimate also permits coupling the finite original words to eleven common uniform digits by flooring.

## Physical hypotheses by sheet type

**Eleven-specific entry and weights.** Here $w(y)=11!\prod_{i=0}^5 g_i(y)$ in ordered difference coordinates, and the denominator at each base is $B^{11}-B$. The ordered-word multiplicity is $B^6 w(q_0/B)+O(B^5)$ on the strict source. Thus the original Riemann sum over its five-dimensional lattice has exactly the scale $B^{11}$. Strict fixed-length leading entry prefixes have integral, unimodular physical branches off null original walls. On those branches the entry lattice has the full physical congruences $X_i\equiv(i-6)B\pmod6$ for one-based $i$. Keep the prefixes and their original marks separately even when their leading entries overlap. Each fixed prefix costs a bounded number of updates, which vanishes after scaling by $B$. [Original prefix transfer](../en11b/README.md) · [Difference weights](../../odd_gap_transport.md#exact-states-and-their-input-weights)

For positive height $r$, the rightmost negative signed difference $-q_1$ starts the borrow and gives $B-q_1$. Every subsequent negative column gives $B-q_i-1$, and the intervening zero columns give $B-1$. The first positive column consumes the borrow and gives $q_r-1$. The remaining positive columns keep $q_i$. The graph sorts these actual digits, so the height, adjacent sort and interlevel guards apply even on a tied or unpaired sheet.

| Incoming sheet | Eleven-specific hypothesis checked on a route | Dimension-independent calculation |
|---|---|---|
| Strict paired | Positive height and the actual borrowed sort and interlevel guards at every phase. The full sort and paired update then agree. The accepted first two stages and the cube provide examples. | Pull back all guards, compare declining ratios and push the original residue cosets. |
| Tied | The earlier borrowed digit label wins equality in the full sort. An interlevel tie at zero slack is paired. | Test the exact comparison on each inherited coset. A tie on the image can carry full conditional weight. Resolving it preserves that weight. |
| Unpaired or shorter height | Sort all eleven borrowed digits. If $r<5$ is the positive height, there are $11-2r$ central copies of $B-1$ and the end offsets use $r$. Guard the zero trailing differences and the positive last one. A return may include paired and unpaired phases. | Apply the same affine guard and first-exit proof. A single admitted block can account for arbitrarily many physical pairing failures. |

On the physical lattice a strict guard can be written as weak after subtracting six. The graph checks weak interlevel flags separately. Its ambiguous-comparison stop on a genuine finite chart can occur only at a proper original wall, because a leading identity on the chart annihilates $M\Lambda$. Such walls have zero limiting original weight. The fixed-grid sampler can hit them and records them as stops. Short-switch and stage caps likewise leave charged stops.

Each finite graph history can be certified by pulling back the comparisons made in its first $m$ decisions, with $m$ fixed. On its positive-weight original cells those comparisons have a strict leading sign or the exact constant sign above. The return test on the lattice directions and the test $Ad=d$ support the translation. Minimizing the physical guard ratio gives the first failed repeat. A selected primitive switch costs exactly one update with its own guard and affine map. Induction gives the physical realization and original-weight limit of each admitted finite history, including stopped histories.

For a fixed finite observation after leading entry, first restrict to entrance by $K$ and prove the finite-history identity on compatible fine base sequences. Then send $K$ to infinity. The omitted original weight is $\mu(D\setminus D_K)$, where $D_K$ increases to $D$. The final clock requires the separate stopping and refinement hypotheses.

The count also covers rank zero, where the guard has one realizable remainder and the image is a singleton. Task 097's implementation asserts a nonempty direction lattice when splitting a macroscopic return. A reachable rank-zero macroscopic section at positive radius would require implementing the singleton case. Reachability of such a section remains open.

**Eleven-specific terminal hypothesis.** A leaf with identically zero leading centered mark has bounded actual $X$ on its fixed route. For sufficiently large bases it lies in the central pair-permanent region. The centered clock theorem bounds its remaining time, and the first cyclic word has $Z=(5,4,3,2,1)/6$ in the limit. First cyclic difference and first cyclic word times differ by at most one. An alternative is a proved positive-radius collar, using the actual first hit or a certified difference in stopped times. The centered remaining-clock estimate gives at most $16\eta/3+o(1)$ at normalized radius $\eta$. On a completed central record, the central cycle envelope bounds every centered cyclic difference even before the terminal section. If a translation block contains first capture, its same-phase states then progress by nonzero integral $d$ within that bounded envelope. At most $O(1)$ repeats of that block and each later block can lie after capture, so a finite record preserves the bounded time error. [Paired permanence and centered clock](../../eleven_centered_clock.md#4-actual-paired-permanence-and-capture)

## Mass accounting and period

On each compatible base refinement, disjoint terminal records cover $\sum_\ell I_\ell^{-1}\int_{P_\ell}w(y)dy$ of the original denominator. Subtract this sum from $m(D)$ to charge the residual. All three sheet types transfer on admitted finite records. A population clock law still needs almost-sure progress to each collar and agreement of the compatible refinements after projection.

The accepted first two stages transfer on the fixed-entry domains. Task 130 proves an original mass of at least 13.3940217690% through four leading steps. Task 138's accepted fixed-entry certificates through step 19 raise the proved lower bound to 16.2038159009%. The original direct-entry route through the second-stage interlevel loss has an exact rational mass of approximately $4.66615569\times10^{-8}$ in every class. Its onward residual beyond the certified continuations is open. [First unpaired section](../stage_graph/README.md#first-two-stages-and-first-unpaired-section) · [Early original entrance weights](../en11b/README.md) · [Extended entrance certificate](../e1_basins/README.md)

The paired cube has original mass approximately $3.5119795377818166\times10^{-25}$ in every class. The repeated-failure whole box has original mass approximately $8.51238978209143\times10^{-14}$ in classes 5 and 11. These are disjoint closed controls with exact rational integrals. The bounded alignment collars and the chart edges have vanishing limiting original weight after all three box cohorts are included. A fixed cap on the number of pairing failures eventually removes the entire box, since its 20-step unpaired return occurs $N\ge B/55000-O(1)$ times. The finite route compresses those failures and carries the full limiting original box mass. [Exact mass ledger](control-weights.json)

For the cube the six tables repeat by base mod 6, using positive-pairing covariance and the certified physical guards. They project to a common law on each mod-12 class. The repeated-failure box requires an original lattice pullback before pooling. On $(c_0,z,r_1,\rho)=(0,1,17,2)$, the checkpoint labels $r_2=603,639$ each have conditional box mass $1/28512$, approximately $2.9855463601611357\times10^{-18}$ of all original inputs. The two profiles within L have exact weights $2/3,1/3$.

| $B\bmod24$ | $r_2=603$ | $r_2=639$ | Equal-weight projection |
|---:|---|---|---|
| 5 | M2 | L | M2 and L |
| 17 | L | M2 | M2 and L |
| 11 | L | S | L and S |
| 23 | S | L | L and S |

The checkpoints retain mod 24 when the label is kept. Pooling the partners gives the same projected clock jointly with the original normalized differences. The twelve-shift is a bounded translation of the original chart lattice that matches partners and preserves the clock up to a bounded shift. The full-box pullback shows all seven final affine graph weights are constant inside each arm mod 12. The box has terminal mark $Z=(5,4,3,2,1)/6$. [Box chronology, translation and full table](../box_structure/README.md)

The box table gives a signed difference on the **clock marginal**. For class 5 minus 11 the pointwise difference of lane-averaged clocks is

$$
\frac{5233465}{2512705536}
-\frac{2687455}{1256352768}P(x)
-\frac{707225}{628176384}R(x),
$$

where $P,R$ are the box chart's late-return forms. The exact minimum over its 32 vertices is $148149493/9213253632000>0$. All seven clocks lie in $(0,1)$ on the box, so the bounded reader $F(t)=\min(1,\max(0,t))$ has the same mean difference there. Integration gives a box-conditional difference of approximately $1.79659361890\times10^{-5}$ and an original-weight contribution of approximately $1.52933051641\times10^{-18}$. The ledger gives both exact rationals. The remaining original mass can contribute with either sign. [Signed marginal calculation](control-weights.json)

The independent determinant reader reproduces the checkpoint weights at every eligible residue mod 72. At all eight physical witnesses near $90{,}000{,}000$, the frozen graph follows the physical route when each lane is selected from the actual incoming remainder. For these $z=1$ witnesses it uses a 10-step failing word that covers 43,976 of the 43,977 failures. It preserves exact lattice membership through the exits. Its rational leading clock equals the predicted affine profile at each witness, and it reaches the bounded central section 20 to 62 updates before the independently replayed first cyclic word. The box certificate uses a 20-step word and covers the same chronology. [Exact weight rederivation](control_weights.py) · [Exact selected lanes](box-route.txt)

For a general finite refinement, subdivide the original charts in common and group complete leaves by their final affine clock and terminal mark. Keep the original density and index weight. Equal grouped weights prove equality of the projected joint laws and hence of the clock marginals. Conversely, distinct affine graphs coincide only on original null sets unless they are identical on a cell. On a cell of positive density the grouped sums thus give a local necessity for the joint law. The population-scale period remains a separate proof obligation.

## Class-3 diagnostic

The 60,000-vector common-rank panel resolves a shrinking class-3 minus class-5 discrepancy between the physical and model clocks. Pooling the two subresidue pairs, the strict reader estimates $0.00201033$, $0.00054978$ and $0.00012841$ at the three successive scales. Its 95% simultaneous half-widths are $0.00007280$, $0.00006417$ and $0.00006365$. The five-phase reader gives the same conclusion. The highest-scale residual is about 6.4% of the smallest-scale gap. The measured contrast is consistent with a finite-base effect, while the analytic results above govern the scope of the limit. The complete ledger gives both individual arms, both mod-24 subresidues, paired uncertainty, original-denominator exclusions and their overlaps. [Frozen protocol](PROTOCOL.md) · [Source hashes](FREEZE.sha256) · [Diagnostic results](diagnostic-results.json)

## Reproduction

Run from this directory. The cube and marked clock wrappers build in private snapshots. The diagnostic wrapper verifies the frozen files, cross-checks the physical clock against column subtraction, and holds the shared compute lock throughout the panel.

~~~sh
python3 -B ../stage_graph/exact_cells.py | cmp - ../stage_graph/cell-lotteries.json
python3 -B ../box_structure/marked_weights.py | cmp - ../box_structure/marked-weights.json
python3 -B ../box_structure/snapshot.py --small marked_clock.py | cmp - ../box_structure/marked-clock.json
python3 -B control_weights.py | cmp - control-weights.json
python3 -B box_route.py | cmp - box-route.txt
shasum -a 256 -c FREEZE.sha256
python3 -B run.py | cmp - diagnostic-results.json
~~~
