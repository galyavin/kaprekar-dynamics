# Exact continuation from the strict regular chamber

The population is the original eleven-digit nonconstant word law, restricted to the event $D$ of finite leading entry into the strict chamber \(\Pi\). Arithmetic lanes are the compatible cosets inherited through each physical exit. The chart calculation separates leading geometry, exact physical offsets and the varying base.

**The two-keep bound is false.** One strict original source cell has positive-weight history $543222110$ in six classes modulo 12. The three keeps occur at ranks two, two and one. All six witnessed lanes reach the centre in eight macroscopic stages. The rank-three exclusion and the uniform bound to rank one hold on the whole entrant law. Almost-sure collar completion from all rank-one arrivals remains open.

## Replay premise

The Task 135 replay ran from a private snapshot with the shared compute lock. The output matched the archived rank replay byte for byte. The exact sample counts are 119,457 with no keep, 101,316 with one rank-one keep, 7,651 with two rank-one keeps, 8,467 with one rank-two keep, 3,106 with one keep at each of ranks two and one, and three censored before entry. These are sample counts. [Replay output](premise.out) · [Frozen replay driver](../../../docs/drafts/research-direction/astra-strategy-13/ranks.py)

## Rank three and above

**Proved by exact chart coverage.** Almost every leading entrant, on every compatible inherited lane, follows rank history $5\to4\to3\to2$ through its first three macroscopic exits. The first and second stages are available without short switches. The third is available after at most 25 short updates in even bases and 45 in odd bases. Its labelled words have lengths 5, 10 or 20, including words with an unpaired phase. Every accepted arm of the first-exit partition loses rank at the third exit.

| Root residues modulo six | Rank-three word cells per root | Rank-three first-exit cells per root | Word cells containing an unpaired phase per root | Maximum shorts before the third stage |
|---|---:|---:|---:|---:|
| 0, 2, 4 | 1,525 | 3,515 | 25 | 25 |
| 1, 3, 5 | 1,665 | 3,855 | 30 | 45 |

These exact counts refer to the verifier's partition and arithmetic sheets. The test starts with the six full entrance lattices and uses the base refinements described below. Finer lane weights remain inherited throughout. Null leading walls of the finite partition are charged at the original input. [Exact coverage output](rank-three.out)

## Access to rank one

**Proved by an exact arithmetic overcover.** Every inherited lane reaches rank one in at most six macroscopic stages, with at most 230 total short updates in even bases and 246 in odd bases. The stage bound consists of the first three exits and at most three rank-two exits. The rank-two check enlarges each fixed-base rank-two lattice to the physical integer lattice in its real tangent plane. The saturated coset contains the inherited one. The stage's initial label choices, drift and rank test depend only on its real tangent and its offsets modulo that tangent. Saturation can add subsequent exit remainders, and the checker verifies all of them. A bound on this enlarged search therefore applies on the inherited lane.

| Base parity | Shorts to stage three | Shorts to the first rank-two stage | Shorts after one rank-two keep | Shorts after two rank-two keeps | Sum |
|---|---:|---:|---:|---:|---:|
| Even | 25 | 85 | 80 | 40 | 230 |
| Odd | 45 | 78 | 73 | 50 | 246 |

The separate inherited-lane run closed with the same maxima. It changes two-dimensional coordinates and retains the exact incoming lattice, offsets and base progression. Every allowed history through the first rank-one arrival is $54321$, $543221$ or $5432221$. All three occur in each base root on strict nonempty projected cells with compatible original preimages. The last pattern contains two successive rank-two keeps, and was absent from the 240,000 archived draws. Its strict projected cells have positive original volume because the rank-two projection is full rank onto their affine image. Each compatible path has positive lane weight. [Inherited-lane coverage](rank-two-inherited.out) · [Saturated coverage](rank-two-saturated.out)

## Geometry and physical applicability

Use $X=6q-B(5,4,3,2,1)$. Each record stores the leading affine map to $X/B$ and a fixed-base lattice. On a fixed root chart the original leading density stays marked. The initial \(\Pi\) region is the intersection of the six-gap simplex, the within-level comparisons and the interlevel comparisons over one five-phase rotation. Its closure has 77 vertices. Exact evaluation of the convex gap energy at those vertices gives the maximum $1/30$, attained only at $(0,1/5,1/5,1/5,1/5,1/5)$. Strict entry has energy below $1/30$. A finite physical lift of a compact subset of a strict entrance chart therefore satisfies Theorem B′ with a positive margin. Up to $S_B$, every gap has a positive lower bound proportional to $B$. After $S_B$, the centered permanence theorem supplies positive pairing. This excludes zero-inner updates on these physical continuations. B′ also bounds physical time. Further work on the stage rules is needed to control the number and kind of recorded switches. [Strict chamber](../stage_graph/README.md#first-two-stages-and-first-unpaired-section) · [Theorem B′](../theorem_bprime/README.md#1-statements-and-conventions) · [Centered permanence](../../eleven_centered_clock.md#4-actual-paired-permanence-and-capture)

On a strict finite entry prefix, the physical affine derivative is the unimodular derivative of the leading branch. The centered output lattice at entrance is thus $X_i\equiv iB\pmod 6$, intersected with the transported source chart. The prefix costs $o(B)$ and transports an absolutely continuous original density. Therefore a null set of entrance positions for each fixed arithmetic branch is null for that prefix. A countable union covers all finite prefixes and finite records. [Finite-route counting and transfer](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md#1-the-restricted-law-and-its-arithmetic-parameter)

## Guard partition and arithmetic

At a chart, the verifier tests one exact rational interior representative and pulls the primitive digit comparisons back to the chart. A strict comparison is certified by its positive leading form throughout the cell interior. A structural tie makes the leading form identically zero. Its row annihilates both the fixed-base lattice directions and the joint base-shift direction. The primitive uses the exact affine intercept, preserving the physical $>0$ versus $\ge0$ distinction on the lattice. Exact equality of two borrowed digits leaves the numerical update invariant, and increasing borrowed-list label fixes the sorting convention. A coincident quotient with an identical leading form uses the least exact affine intercept. A proper leading tie of distinct forms is a hyperplane in the original chart and has zero source mass.

The convex source cell is partitioned by the pulled-back primitive comparisons. For each repeated word, the checker collects every declining guard. Normalize each by its positive decrement and group it by its leading time form. Partition the cell by the minimum of all distinct leading time forms at once. Every simultaneously defined cell uses the same word, with a first exit given by its own winning guard. All decisions are rational. When a candidate word is used, the verifier checks that its leading return is the identity on the entire current affine chart, and the stage reader checks translation on all the fixed-base lattice directions and that its linear part fixes the drift. The physical induction of all guards up to the first exit follows from their affine change per return. Bounded short updates use the full borrowed sort and its exact affine map. No sampled geometry decision is used as a certificate.

At a fixed-base lattice with direction space $V$, the exit derivative is $I+da/s$, where $d$ is the drift, $a$ the exit row and $s=-ad>0$. Its ambient kernel is \(\mathbb R d\). Compute the outgoing rank from the exact columns of this map on the residue sublattice. This applies the span criterion, since the residue sublattice has the same real span as the incoming one. Positivity and exact tie resolution are checked before applying the criterion. [Stage transport](../stage_graph/README.md#section-and-renewal-calculation) · [Rank criterion](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md#1-the-restricted-law-and-its-arithmetic-parameter)

For a varying base, represent the physical incoming lattice at $B=b+kP$ by $p+kt+L\mathbb Z^r$. Under a short step, push $p,L$ through the affine update and $t$ together with the base shift $P$ through its linear part. At a macroscopic exit, put $h=\gcd(s,aL_1,\ldots,aL_r)$ and \(\delta=at+a_B P\), where $a_B$ is the base coefficient of its guard. The remainders for a fixed base form one coset modulo $h$, with $s/h$ choices of equal weight. The coset itself varies with the base index $k$, through \(k\delta\bmod h\), so refine $P$ by $m=h/\gcd(h,\delta)$. The extended Euclidean algorithm finds an $L$-displacement \(\ell\) with $a\ell+m\delta\equiv0\pmod s$. An exit in any of the $m$ refined base subclasses preserves its remainder under the combined shift, and the exiting slope is

\[
  mt+\ell+\frac{m\delta+a\ell}{s}d.
\]

The program uses this refinement even on an interlevel exit, so an integer representative stands for a proved progression in $B$. It covers compatible fine classes without assuming projected clock invariance modulo six or twelve.

## Low-rank chart images

For existence and universal bounds, different histories with the same rank list and exactly the same normalized location and physical sheet have the same available continuations. With line compression enabled, use the first nonconstant centered coordinate as the parameter of a rank-one image. The image of a convex source cell is the interval between its two exact vertex extremals. Intervals already assigned work are covered by their union. With plane compression enabled, use two independent coordinates as parameters of a rank-two image. The image is the exact convex hull of the projected vertices. Reparametrizations retain the same affine lattice, base progression and rank history. They do not enlarge the lane space or assign a new probability to the projected geometry. Any proposed weighted counterexample found on a compressed chart would need its original preimage to supply its original mass.

In either compressed chart, the code reuses a simplex for affine guard intersection. A proper tie in the active parameters remains a proper tie on the original chart because the projection has full rank in those parameters. The synthetic simplex is only a domain representation.

The optional saturation flags enlarge a fixed-base direction lattice to its intersection with $6\mathbb Z^5$ in the same real span. The affine coset still contains all inherited physical lifts. Any universal bound certified on this enlargement applies to the original lanes. A witness on a newly admitted coset needs reconstruction on an inherited lane. The central-ray option widens a homogeneous interval to the same half-ray within the permanence radius. It likewise enlarges the set checked for a universal bound.

The short-update lookahead rejects a proposed macroscopic period only if its leading position cannot return or its next two labelled phase words differ. A positive macroscopic duration would force two identical copies for all sufficiently large bases on the current chart. Every primitive comparison used to reject these periods is kept in the partition. Where the test leaves a possible period, the checker uses the original stage search.

## What a bounded short orbit would prove

The leading position at a positive-length stage exit is periodic. The certified word has affine leading return $A y+b=y$ on its incoming chart and $Ad=d$. Its inequalities remain weakly valid at $y+td$ through the leading first exit. The full positive-inner leading map is continuous, so the same word sends the endpoint back to itself. Subsequent bounded primitive updates follow that leading cycle. Along it the physical sort may keep changing because of exact offsets at structural ties.

There is a sufficient condition for a next macroscopic stage. Start at a noncentral periodic phase on a fixed strict chart with tangent space $V$. Suppose a finite sequence of primitive updates returns to that same leading affine phase and also returns to the same physical affine coset modulo $V$. Write the primitive word's exact map as $AX+bB+c$. The leading return gives $A|_V=I$ and $Ay+b=y$ on the chart. The same affine coset means the drift $d$ is in $V$. These statements also give $Ad=d$ and make the drift independent of advancing the compatible base progression. It cannot be zero, since then the physical lifts would already cycle outside the central collar in contradiction with B′. Every comparison whose leading value vanishes identically annihilates $V$ after its phase pullback and hence has zero drift. There must be a declining comparison, since otherwise the word could translate by the nonzero $d$ indefinitely in the finite physical state space. The first such exit has positive leading time and drops rank by one.

In particular, an unbounded run of short updates at positive radius under a search allowing arbitrary finite word lengths must avoid recurrence of the physical normal coset at its leading phase. The quotient of the integer lattice by its saturated tangent sublattice is discrete, so a bounded orbit of these offsets is sufficient for finding a macroscopic return. For the fixed 120-word search, the same argument applies when the recurring displacement is found within its search length. B′ bounds physical time to a collar and excludes zero drift. The remaining short-update obligation starts on rank one.

On a rank-one section in the permanent paired region whose leading line passes through the centre, a rank-dropping stage ends at the centre. The active within-level comparison and the return are homogeneous in centered leading coordinates. The nonzero tangent drift lies along that line, so the rank-one projection to the active exit kills the entire line. This covers any such section once its offset orbit yields the recurrent word just described.

## Positive original weight of the third keep

The inherited-lane search on descendants of two rank-two keeps uses neither saturation nor ray widening. An uncompressed reconstruction puts the six witnesses on one common strict convex source cell $C$, whose closure has 12 vertices in the original six-gap simplex. The exact Dirichlet mass of $C$ is in the replay file. The starting population is direct entry into $\Pi$ on $C$. The inherited joint lane weights in the table count one witnessed path in each stated base class, and the rounded bounds are on the original denominator $B^{11}-B$.

| $B\bmod12$ | Inherited lane weight on $C$ | Proved lower bound on original mass | Short updates through the eighth stage |
|---:|---:|---:|---:|
| 6 | $1/155520$ | $2.375446984\cdot10^{-15}$ | 141 |
| 7 | $1/124416$ | $2.969308730\cdot10^{-15}$ | 51 |
| 8 | $1/311040$ | $1.187723492\cdot10^{-15}$ | 86 |
| 9 | $1/124416$ | $2.969308730\cdot10^{-15}$ | 106 |
| 10 | $1/217728$ | $1.696747845\cdot10^{-15}$ | 131 |
| 11 | $1/139968$ | $2.639385538\cdot10^{-15}$ | 111 |

For $B\equiv6\pmod{12}$ the periods are $5,10,20,10,20,10,30,20$. The third word contains one unpaired phase. The keeps are stages four, five and seven. In centered integer coordinates each row below repeats a verified word with drift $d$. A stage with winning guard $G(X,B)$ and decrement $s$ exits after $1+(G-r)/s$ returns, at the listed remainder $r$.

| Stage | Rank | $d$ | $G(X,B)$ | $s$ | $r$ |
|---:|---:|---|---|---:|---:|
| 4 | 2 | $(-60,-48,0,-24,-12)$ | $2X_1-X_2-X_3+2X_4-X_5-60$ | 108 | 54 |
| 5 | 2 | $(-120,-96,-36,-120,-204)$ | $2X_1-X_2-X_3+2X_4-X_5-96$ | 144 | 96 |
| 7 | 1 | $(-96,-60,0,-96,-60)$ | $7X_1+2X_3-7X_4+3X_5-78$ | 180 | 36 |

The eighth stage drops to rank zero and its leading output is identically zero. Its endpoint is physically bounded on the inherited lattice, uniformly in the base. The analogous direct replay closes each of the six fixed witnessed prefixes at the eighth stage. Thus these counterexamples reach every positive terminal collar. [Exact class replays, source vertices and mass fractions](third-keep-classes.out) · [All eight stage maps at the class-6 witness](third-keep-route.out)

## Direct collar estimate on permanent paired rays

The centered clock gives a strict radius decrease through any positive elapsed macroscopic time in the permanent paired region. Until the first $R\le120$, there is a lower hit costing at most 32 updates and reducing integer $R$ by at least six. The intervening increase is bounded by 384. After the first $R\le120$, the physical bound is 504. Thus at time $T$,

\[
R(X_T)\le \max\left(504,\ R(X_0)-6\left\lfloor T/32\right\rfloor+384\right).
\]

On a certified record with scaled duration $\Delta$, the leading radii therefore satisfy $r_{\mathrm{out}}\le\max(0,r_{\mathrm{in}}-3\Delta/16)$. On any fixed homogeneous ray in this region, both $\Delta$ and $r_{\mathrm{in}}$ scale linearly. Each certified positive-duration arrow then has a strict contraction ratio independent of the amplitude. A finite closed graph of such ray and offset states, with finite short continuations, would prove access to all positive collars even if a path had infinitely many stages. The remaining proof is to control the reachable rank-one offset orbits and the possible nonhomogeneous or unpaired rank-one sections before using this contraction.

## Five positive-weight histories

Every observed signature in the archived replay occurs on a strict open original-input cell. The rational source centers and guards in the witness ledger certify a whole gap box for each row. For a box of half-width $\delta$ in the first five gaps, the sixth gap changes by at most $5\delta$. The exact lower bound on original law is $11!(2\delta)^5 \prod_i(g_i/2)$ times the inherited lane probability, because the selected box leaves every gap at least half its central value. The bound counts direct starts in the chamber on the denominator $B^{11}-B$.

| Rank history | Compatible bases for the witnessed lane | Proved lower bound on original mass |
|---|---:|---:|
| 543210 | $B\equiv0\pmod6$ | $4.4\cdot10^{-21}$ |
| 5432110 | $B\equiv0\pmod6$ | $1.89\cdot10^{-20}$ |
| 54321110 | $B\equiv0\pmod6$ | $2.66\cdot10^{-24}$ |
| 5432210 | $B\equiv0\pmod6$ | $9.57\cdot10^{-28}$ |
| 54322110 | $B\equiv2\pmod{12}$ | $4.19\cdot10^{-26}$ |

For both rank-two keep witnesses, the ten-update third-stage word has one unpaired phase. Both return to paired words, followed by the rank-two keep and a drop to rank one. The second witness then has a rank-one keep before the central exit. The exact fractions, drifts, guard rows, remainders and affine lattice columns appear in the ledger. [Witness reader](witnesses.py) · [Exact witness ledger](witnesses.out)

## Reproduction

Run from the repository root with Python 3.9 or later and a C++17 compiler. The runner compiles a snapshot of the two accepted stage-graph sources and the chart reader, acquires the project compute lock and writes results to stdout. No shared input is modified.

```sh
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 3 --jobs 6 --limit 40000
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 2 --jobs 6 --limit 42000 --compress-planes --saturate-planes
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 2 --jobs 6 --limit 50000 --compress-planes
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 0 --jobs 6 --limit 80000 --compress-planes --compress-lines --paired-fast --rank-one-min-rank-two-keeps 2
```

The last command locates third-keep prefixes and returns status 2. The paired-sort shortcut orders each borrowed level separately when the leading point is in the permanent paired region. The full adjacent digit and interlevel guards are still retained for verification.

The direct class-6 reconstruction uses the original geometric chart. The first command returns status 2 with the weighted witness. The second follows it through the central exit.

```sh
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 0 --offsets 0 --paired-fast --prefix kaprekar-investigation/eleven_access/z1_stages/third-keep-prefix.txt
python3 -B kaprekar-investigation/eleven_access/z1_stages/cover.py --rank 0 --offsets 0 --paired-fast --prefix kaprekar-investigation/eleven_access/z1_stages/third-keep-prefix.txt --keep-limit 12
```

For the five witnesses, compile `chart.cpp` with copies of `lottery.cpp` and `graph.cpp` from `stage_graph/` in the same private directory. Give the resulting executable path to `witnesses.py` and compare stdout with `witnesses.out`.
Pass `--third-keep` for the class-6 route replay and compare with `third-keep-route.out`.

To reproduce the archived premise independently, copy `docs/drafts/research-direction/astra-strategy-13/ranks.py`, `docs/drafts/research-direction/lens9-breakthrough/` and `kaprekar-investigation/eleven_access/stage_graph/graph.cpp` to a private snapshot with their relative paths intact. In that snapshot, link `.task-monitor/compute.lock` to the shared lock and run `python3 -B docs/drafts/research-direction/astra-strategy-13/ranks.py`. Compare stdout with `premise.out`.
