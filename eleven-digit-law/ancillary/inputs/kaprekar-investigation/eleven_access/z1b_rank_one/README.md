# Exact rank-one continuation

The input is the set of inherited rank-one sheets reached from strict leading entry into \(\Pi\). The weights in any finite route are those of the original nonconstant eleven-digit words on denominator \(B^{11}-B\), with the conditional lane splits of Task 097. The checker proves universal coverage by enlarging each arithmetic sheet. Weight calculations always use the inherited source sheets.

## Result

Almost every original leading entrant into \(\Pi\), on every compatible inherited lane, reaches each positive terminal collar in finitely many certified records from its rank-one arrival. The graph is finite, its outer subgraph is acyclic, and its permanent paired rays have a uniform contraction. The limiting unresolved original weight on the entrance event is zero.

The six exact premise replays reproduce all three rank histories `54321`, `543221`, and `5432221`, including the even and odd short-update bounds. The uncompressed class-6 control reproduces `543222110`, with its unpaired third word, all eight drift and exit records, and a central endpoint. [Exact ingress and closure summaries](verification.out) · [Class-6 physical replay](class6-replay.out)

| Base root mod 6 | Inherited rank-one arrivals | Initial canonical sections | Reads to graph closure | Final sections | Distinct edges between sections |
|---:|---:|---:|---:|---:|---:|
| 0 | 178,133 | 25,214 | 26,513 | 26,409 | 164,139 |
| 1 | 252,289 | 37,818 | 39,380 | 38,997 | 210,713 |
| 2 | 178,160 | 25,214 | 26,513 | 26,409 | 164,355 |
| 3 | 252,700 | 37,804 | 39,373 | 38,994 | 210,756 |
| 4 | 178,173 | 25,214 | 26,503 | 26,403 | 164,103 |
| 5 | 252,560 | 37,786 | 39,352 | 38,977 | 210,528 |

These are exact counts in the verification partition. All 1,292,015 first rank-one arrivals have homogeneous leading lines. The initial sheets have 14 distinct primitive tangent directions in each even root and 15 in each odd root, with normal-sheet base period six. Every reached continuation section still has a homogeneous rank-one line, apart from terminal central exits. Thus nonhomogeneous sections have zero limiting original weight on this rank-one continuation. The exact reader checks all borrowed comparisons at the outer sections and their short switches. Every verified rank-one macroscopic word was paired. The outer DAG proof controls the regions before the permanent-pairing condition applies.

For every root the exact minimum on permanent ray arrows is \(\Delta/r=5/48\). The physical radius inequality therefore gives a contraction factor at most \(251/256\). A path leaves the acyclic outer graph for a permanent ray or reaches a central exit. Every subsequent nonterminal inner arrow contracts, so a finite record reaches any specified positive collar.

## The affine normal quotient

Use centered coordinates \(X=6q-B(5,4,3,2,1)\). An inherited rank-one sheet at bases \(B=b+kP\) has the form

\[
X=p+kT+nL,\qquad k,n\in\mathbb Z.
\]

Let \(u\) be the primitive integer direction of \(L\) with positive first nonzero coordinate. Choose an integral row \(c\) with \(cu=1\), and let \(N=I-uc\). This is an integral projection with kernel \(\mathbb R u\). Its restriction to \(6\mathbb Z^5\) has image \(6(\ker c\cap\mathbb Z^5)\), and the physical tangent lattice is \(6u\mathbb Z\). Define

\[
a=NT/P,\qquad d=Np-ba.
\]

Every inherited physical point satisfies \(NX=Ba+d\). Fix \(\beta=B\bmod6\) and choose an integral representative \(p_\beta\) of \(X_i\equiv i\beta\pmod6\). For fixed \(B\), all physical solutions of this exact normal equation constitute one coset of \(6u\mathbb Z\) when a solution exists. Starting at one compatible base, the allowable base shifts are precisely multiples of

\[
P_0=6\operatorname{lcm}\bigl(\operatorname{den}(a_1),\ldots,\operatorname{den}(a_5)\bigr).
\]

Indeed a shift must preserve \(\beta\) and satisfy \(N(X'-X)=a(B'-B)\) with \(X'-X\in6\mathbb Z^5\). The right-hand side already lies in \(\ker c\), so coordinatewise divisibility by six is sufficient. With \(b_0=b\bmod P_0\), a representative at \(b_0\) is \(p_\beta+N(p-p_\beta)+(b_0-b)a\). The new slope is \(P_0a\). These formulas remove redundant base refinements and tangent offsets while retaining the full physical normal offset, including the leading affine intercept. At rank zero use \(N=I\) and the same formula, with an empty tangent lattice.

A rank-one graph section comprises this canonical sheet and an open interval along the normalized line, with the first varying coordinate as parameter. Each union of intervals is maintained exactly. The normalized line is through the centre if and only if \(a=0\). Such a line is cut at the centre and at the two boundaries of \(R(X)/B<1/7\). If a strict interval on an inner half-ray is reached, the checker admits that entire open inner half-ray. On that domain, the leading paired map and its within-level comparisons are homogeneous in \(X/B\). The full paired map remains valid for all sufficiently large bases on each compact subinterval. Every leading point added to a ray is a positive scaling of its periodic seed within the same paired domain. Intervals on nonhomogeneous or outer sections are retained as given. A rank-zero noncentral point is read as a singleton.

## Exact graph certification

The ingress pass runs the accepted Task 137 inherited-lattice checker through rank two. It records each candidate rank-one arrival before inserting it into the continuation graph and compares the six root summaries exactly with the accepted report. Thus all three verified predecessor histories and all inherited remainder choices feed the new graph. The projected geometry has full rank in its active parameter, so discarded proper walls pull back to null original cells. A guard tied identically on a chart is decided with its exact normal offset.

Each continuation read first uses the accepted primitive borrowed sorter and the stage reader. At an interior rational parameter it records the exact guards for its bounded initial switches and repeating word. It partitions the one-dimensional interval exactly by the reported guards, and then by the race among all structurally distinct declining guards. All compatible outgoing remainder lanes are included. A duration is the winning guard's leading affine value times the verified word length divided by its decrement. The same partition retains the exact rank and outgoing chart. Full digit sorting is used whenever the permanence test does not supply the paired shortcut. A primitive interlevel failure, a base-dependent leading exit, a nonhomogeneous line and a singleton are therefore read by the same physical rules, including simultaneous structural ties. Graph closure requires the exhaustive case calculation on these canonical sheets.

For a finite graph with no recurrent outer component, each path either ends centrally or eventually stays on homogeneous permanent rays. Each positive-duration certified arrow on an inner ray has elapsed time proportional to radius. The checker compares this ratio exactly and checks that every nonterminal inner arrow stays on an inner ray. Task 137's paired radius inequality then gives

\[
r_{n+1}\le \max(0,r_n-3\Delta_n/16).
\]

A positive minimum of \(\Delta_n/r_n\) over the finite graph implies a uniform strict contraction along inner arrows. Arbitrarily many arrows are allowed. The noncentral short updates in an arrow have bounded physical length. Each starts at a phase of the same leading periodic orbit, because the limit at a first exit of a repeating word is a fixed point of its leading composite by continuity. Thus the bounded shorts preserve the leading radius until the next macroscopic block.

For a finite outer graph with cycles there is a second sufficient test. Every finite physical prefix that avoids the target collar has uniformly bounded normalized elapsed time. Theorem B′ gives this from the strict \(\Pi\) entrance sublevel, followed by the centered physical clock. Equivalently, apply that bound to a physical reconstruction of each finite prefix on its compatible base progression before taking the limit. If the prefix has already reached first capture, each fixed nonzero-drift block can spend only a bounded number of repeats on the bounded centered cycle, and each short part is bounded. An infinite path must consequently have \(\Delta_n\to0\). In a finite graph, each edge used infinitely often must have a zero duration at an endpoint of its allowed interval, since its duration is affine and strictly positive in its interior. For such edges record the one-sided approach parameter, its output coordinate and the output approach side. Acyclicity of this finite zero-duration boundary graph excludes an infinite outer path. One can ignore transitions into rank zero in this last test. The rank cannot rise again, and each positive stage on a fixed point has a fixed strictly positive duration.

## Clock tolerance

For any fixed \(0<\eta<1/7\), a certified visit to the strict collar \(r<\eta\) gives membership in the permanent paired source for all large bases. The remaining cyclic-word time is at most \((16/3)\eta B+C+33\), where \(C=\binom{509}{5}\) is the accepted uniform difference-cycle bound. If a finite recorded visit follows first capture, a nonzero integral drift cannot persist for more than a bounded number of repeats within the centered cycle envelope. All remaining shorts on the fixed route also cost a bounded number of updates. Thus the recorded elapsed time approximates actual first cyclic-word time divided by \(B\) to within \(16\eta/3+o(1)\). The normalized terminal difference mark tends to \((5,4,3,2,1)/6\). Compatibility of the clock marginal across finer progressions within each class modulo 12 remains a separate obligation.

[Task 137 reproduction and bound](../z1_stages/README.md) · [Task 136 physical transfer](../r2_routes/README.md) · [Theorem B′](../theorem_bprime/README.md) · [Centered permanence and clock](../../eleven_centered_clock.md)

## Reproduction

Run from the repository root with Python 3.9 or later and a C++17 compiler. The verifier compiles the accepted C++ sources in a private directory and holds the project compute lock during the six-root calculation. `--state-dir` may point to a private directory for saving validated ingress checkpoints. `--resume` uses those checkpoints for a continuation rerun in the same task environment.

```sh
python3 -B kaprekar-investigation/eleven_access/z1b_rank_one/verify.py
```

The six `closed [root]` ingress summaries match `rank-two-inherited.out`. Each continuation `rank-one closed` summary shows zero outer and zero boundary topological remainder and the exact ratio (5/48). `verification.out` associates the completion summaries with their roots using the unique inherited arrival counts. For the physical class-6 control, build `chart.cpp` with the accepted `graph.cpp` and `lottery.cpp` in a private snapshot and run `z1_stages/witnesses.py` with `--third-keep` as in Task 137's reproduction instructions. Its stdout matches `class6-replay.out` byte for byte.
