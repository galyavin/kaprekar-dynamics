Independent re-derivation of the reductions from the text, with separately written encodings

**Verdict:** PASS WITH ISSUES

No disagreement was found between the paper's claims and the checks. All issues found are editorial or concern how the supplement's checks are scoped. None changes a theorem, a constant or a table entry. Some evidence is statistical or covers only small bases, not symbolic: the 165 kernels, the 33 route words and the nine-step descent for large B. §5 and §7 give the exact coverage.

- Paper: `earlier draft, paper.tex` (SHA-256 `15b10def…f060`), PDF `4bdb55af…d063`.
- Supplement data used for the predictions: `nine_kernel_results.json` (`b8432d44…732d`) and `nine_clock_certificates.json` (`ceb0ca4f…c429`). Both hashes match `RERUN.md`.
- Code and logs: `independent-review/` (file list in §8).
- The paper and the ancillary directory were not modified. The supplement was imported only from a scratch copy, with bytecode writing disabled.

## 1. Method

1. I read the paper and restated every computer-assisted step as a closed statement (§2). I then traced how the theorems use each step.
2. Before opening any ancillary file, I wrote my own transition model and encodings:
   - **Integer map.** The definition is column subtraction of the ascending arrangement from the descending one (`km.py`, `kap9.c`). A z3 lemma (Z1) shows that this equals the borrow-formula digit list for every rank and every B ≥ 2. The SMT encodings then sort all nine digits with a Batcher max/min network, checked by the 0-1 principle, and take symmetric differences (`zmodel.py`). There are no guarded branches and no sorting permutations.
   - **Real map.** The leading map f is a 9-element sorting network. The paper's radius map H is used only after the lemma C⁻¹fC = H is proved (L1).
3. Only then did I read the authors' encodings (§4).
4. I ran fresh-base tests with independent code, at bases the paper does not use (§5).

Substantial runs held the exclusive `fcntl` lock through `withlock.py`.

## 2. Computer-assisted steps

**Result codes.**
- UNSAT: my own encoding's negation is unsatisfiable in z3 for all integers or rationals, with no bound on B.
- verified: exact computation.
- sampled: checked by direct subtraction on a sample.
- exhaustive-small-B: every case checked for the stated bases.

Each UNSAT query also passed a non-vacuity check: its hypotheses alone are SAT. A mutation check (`z_vacuity*.py`) shows that the encoding can detect a violation. All SMT results were also replayed with z3's QF_LRA and QF_LIA solvers under two other random seeds (`z_reseed.py`); all 78 verdicts (39 queries × 2 seeds) were unchanged.

| # | Location | Precise statement (all variables real unless stated) | Use in the proof | Result |
|---|---|---|---|---|
| 1 | App. A (i); §3 ¶2 | Let 1>a≥b≥c≥d≥0 and r_i=H^i r. Suppose (r_i)_1=a for i≤4, (r_4)_1<1, and sort↓D(r_i)=sort↓D(r) for i≤4. Then r_4 ∈ 𝒞 ∪ ℰ ∪ {(1,1,1,1)/3, (3,3,3,1)/5, (3,1,1,1)/5}. | Classifies limit cycles: the first radius and the sorted gaps are constant on them. | UNSAT (1.9 s). Mutation SAT (ℰ starting at 0.72). |
| 2 | App. A (ii) | Let 0<η≤10⁻⁴ and r be ordered as above, with 5/7−η≤a≤3/4, 0≤b−c≤η, \|c−2a+1\|≤η, 100η≤d≤1−a+η. Then (H⁶r)_1=(H³r)_1. | a is eventually constant near ℰ. | UNSAT. Mutation SAT (η≤1/10, d≥0). |
| 3 | App. A (iii) | Under the hypotheses of (ii), together with a≥5/7 and (H⁶r)_1=a: ν(H⁶r)≥ν(H³r), where ν=b−c+\|c−2a+1\|. | ν is monotone on late sections, so it eventually vanishes. | UNSAT. Strict increase is SAT, so the statement is not trivial. |
| 4 | App. A (iv) | Let η>0, W(r)≤η, min D(r)≥100η and H⁹r∉𝒞. Then H⁹r−H⁶r=H⁶r−H³r. | Gives an arithmetic progression, which is impossible without entry. | UNSAT (7.1 s). It becomes SAT if the min D or the W hypothesis is dropped. |
| 5 | §3 ¶1; end of App. A ("direct substitution") | f(x)=C(sort\|2x−1\|) for ordered x∈[0,1]⁴. f(Cr)=C(Hr) for ordered r∈[0,1]⁴. On C(𝒞), f=R. H(𝒞)⊆𝒞. H³=id on 𝒞. The three ℰ phases cycle 1→2→3→1. The exceptional points form a fixed point and a 2-cycle, each with a zero gap. | Exact periodicity once the orbit is in 𝒞 or ℰ. | UNSAT for all (L1b, L1, L2a–c, L3). Points verified exactly. |
| 6 | §4.1 bulk | For all integers B and q with all 12 integer core slacks ≥16: F_B³q = q+(−3,−2,−1,−1). | Drift toward the faces. | UNSAT (Z5). With margin 0 it is SAT. |
| 7 | §4.1 faces | Of the 12 core facets in y coordinates, exactly six have slack decreasing under δ=(−3,−2,−1,−1): A₀, A₁, A₂ (−4) and B₀, B₁, B₂ (−3). R permutes them (A_i∘R=A_{i+1}, B_i∘R=B_{i+1}). Q(4/5) lies on all six. | Face list and the rotation rule "A_i→R^{i+1}, B_i→R^i". | verified (exact) |
| 8 | Table 1, "Exit" column; App. B ¶"Facet and corner" | For each canonical face or intersection and its stated (ℓ,δ): the decreasing slacks that can be the strict first hit (min h/s) from the relative interior. Also each δ keeps its defining slacks at 0; B₂ changes by +4 under the A₂ return and +18 under the AB return; s=30 or 5 at AA; 24d₁+5δ_AB=18d₀. | First-exit structure. | verified. All six exit lists reproduced exactly from the core geometry (`t_exits.py`), as are all the slack identities. |
| 9 | App. B: sorting words, releases (B₀=A₀, B₀=B₁=A₁), slack maps at AA/AB/BB, terminal words at AAA/AAB/BBB, first capture | For large B the integer orbit from a strict compact cell of each route follows the listed words. T_B/B and Z_B then converge to the route's atoms. | Clock and terminal atoms. | sampled: end-to-end direct subtraction on all 33 routes (item 10). Not re-encoded symbolically. |
| 10 | §4.2; App. B "Weights"; `nine_kernel_results.json` | 33 route domains (12 AA, 15 AB, 6 BBB, with 5, 4 and 1 atoms). For each route and β: the weights p_{jβl}, clock rows t_{jl} and terminal matrices U_{jl}. The weights depend only on B mod 5. | The law μ_β. | sampled: 165 kernels × 2 base families of opposite parity, 3,000 orbits per kernel and family (§5.3). |
| 11 | §4.4 table | Route A₀,B₁,AB, radius (18,14,6,3)/25. Times L+c(s−4/5) for c∈{7/2, 9/2, 50/9, 37/6}, with the stated probabilities. Terminals Q((s+4)/6) for the first atom and Q(4/5) for the other three. | Distinctness of the five measures; the lower bound in Thm 3. | sampled at high precision: 1.5 M orbits, 10 fresh bases, all 16 inverse branches (§5.1). Row sums are exactly 1 (verified). |
| 12 | App. B "covariance" | Every core sorting word respects the paired levels. | Explains the period 5 (Thm 4). | sampled: no orbit in items 10 and 11 left the paired region after its first image. |
| 13 | App. C "First access" | Integers. For B≥40 (the paper states odd B=2a+1 with a≥16): if q₀…q₃ are positive (q₄-coordinate >0) and q_{i+1}=F_B q_i, then F_B q₃ is positive. | A positive run that ends at zero has at most 3 states. | UNSAT for both parities with B≥40, and for odd B≥33 (Z3, Z3b). The text's even-B argument (positive→positive→zero is impossible) is UNSAT (Z3c). Both gate characterisations were brute-forced over 18.6 M sources: 0 mismatches. |
| 14 | App. C "Cover" | Integers, B≥40. Take a zero state q₀ of rank 2 or 3 and j∈{1,2,3,4}, with F^i q₀ positive for 0<i<j and F^j q₀ a zero state. Then q₀ lies in the stated union. | At most 308B returning zero states. | UNSAT for all 16 queries (Z4). Dropping the row \|4u−3B\|≤8 gives SAT (witness B=43). Note: j=4 is infeasible for ranks 2 and 3. |
| 15 | App. C one-step; §5 | Integers. For B≥2 and q₄≥1: S_i(g_B(F_Bq)) ≤ S_i(g_B(q))+1 for i=1…4. Hence G(F_Bq) ≤ G(q)+8. | The dispersion bound. | UNSAT, 4 queries (Z2). With +0 each is SAT, so the +1 is attained. |
| 16 | App. C descent; eq. (descent) | Integers. Let B≥40, q positive with P(q)>100, and F^i q positive for i≤9. Then min_{1≤i≤9} G(F^i q) < G(q). | Uniform linear clock. | exhaustive-small-B: every base 101–160 and 173, 191, 199, 211, 223, which is 811,157,145 qualifying sources, 0 violations. Also sampled: 3.03 M stratified sources with B up to 2·10⁹ and 10⁵ near the slow family, 0 violations. The tree reduction itself was re-derived (§3). |
| 17 | §5 constant | T_B ≤ 1232B+7 + 9(4B+8M) + M + 9 + 1 = 1268B+17+365(B+1)·C(104,4), with M=5(B+1)·C(104,4) and cover size 1+24+283 = 308 (times B). | Thm 3. | verified, including each interval count of the cover. |

Hand-proved steps were tested numerically as well (`t1_majorization.py`, `t_misc.py`):
- Thm 1, real: 21,000 exact rational inputs, k=1…7, 0 violations. The circle inequality was tested for general C≤L.
- Thm 1, integer +2: exhaustive for n=3, 5, 7, 9 at small B, plus random inputs for n≤15 with B≤10⁹. The largest excess observed was +1 everywhere, so the theorem's +2 holds; for n=9 it is +1, as App. C states.
- Rational orbits of f: eventual periods ≤ (k+1)! for k≤6. For k=4 with positive gaps, 4,000 of 4,000 have period 1 or 3 and end in 𝒞; none ends in ℰ (a measure-zero set).
- Thm 4 identity: 3,000 paired sources for each k=1…6, 0 failures.
- Injectivity of the number map: verified.
- h≤T_B≤h+1, and T_B=h+1 for transient q: all 9-digit strings at B=4 and B=5.
- Multiplicity leading term: checked.

## 3. Reductions re-derived from the text

- **§3 (finite entry).** Items 1–5 suffice, together with Thm 1, the monotonicity of (Hr)₁=max(a,1−c−d) and the monotonicity of min D. The exact-hit argument (integer slopes) is sound. So is the arithmetic-progression contradiction for the core case. Every hypothesis of (ii)–(iv) is met by late sections for one fixed small η. No upper bound on B or on the entry length is used.
- **§4.**
  - The lattice-count rule is a correct description: accessible remainders e′ carry weight g/s, the outgoing lattice is [M|δ]·ker, and B mod 5 suffices if g divides the coefficient of b at every split.
  - The unimodular-prefix claim holds, because the linear part of C∘H∘C⁻¹ has det ±1 with integer entries.
  - T_B=h+1 for transient inputs, and f(Q(s))=Q(s) on [4/5, 5/6] (checked by hand).
- **§5 and App. C.**
  - First access plus the cover give at most 4 steps per returning zero visit and at most 308B such visits before capture. So by time 1232B+7 the orbit is cyclic or permanently positive. The +7 is generous; 1232B+4 suffices.
  - Exceptional states (P≤100) number at most 5(B+1)·C(104,4), because the four non-maximal gaps have excess Σe_i ≤ 100.
  - Telescoping G∈[0,4B] gives at most 4B+8M decreasing blocks. The arithmetic gives the stated constant.
- **Descent tree** (`nine_clock_audit.py`, read after my own runs). The roots cover all 120 weak gap orders × the two residual cases, which together are exactly P>100 with B≥40 and no upper bound. A node holds the ancestral sorting guards, positivity, and G(q₀) ≤ G(q_i) ≤ G(q₀)+8i. The coverage query computes the true successor with an independent 19-comparator network, which I checked by the 0-1 principle. It asks for a positive successor with G≥G(q₀) outside every child cell. The depth limit is <9. UNSAT at every node therefore implies eq. (descent) by induction. The +8 cap is justified by item 15. The horizon is sharp: on B=34t−9, q=(26t−12, 20t−9, 14t−6, 8t−3), P>100 for t≥51 and the first decrease is at step 9, for t=100 … 10⁵.
- **Thm 4.** The identity F_{B+jm}(q+jv) = F_B(q)+jv follows from the level increments j(m,k,k,…,1,1). The limit argument by matching is sound. For nine digits it applies to core populations; the period-5 claim of Thm 2 rests on the kernel divisibility check, not on Thm 4.

## 4. Comparison with the authors' encodings

| Check | Authors | Mine | Difference |
|---|---|---|---|
| Entry (i)–(iv) (`nine_global_entry_audit.py`) | Radius map H with a 5-comparator 4-sorter; insertion-sort gaps; Z3+cvc5 on QF_LRA | f by a 9-sorter, H only after L1; 5-sorter for the gaps; z3 | None in statement, domain (a<1, ordered, d≥0) or tie handling. Same 12 facets and the same ℰ predicates. |
| Zero returns, one-step (`nine_zero_return_audit.py`) | B=2a+β with a≥20; borrow digits with the central B−1 removed; 19-comparator 8-sorter | Borrow digits for every rank, including q=0; 9-sorter; B≥40 | Their step constraints require every iterate to have q₁≥1, so returns to the zero vector are excluded. This is harmless: F_B(q)=0 only for q=0 when B≥3, since B−1 and B−q₁ cannot both equal the other digits. My encoding allows q=0 and is UNSAT. |
| First access (`odd_zero_access_audit.py`) | Odd B, a≥16; the gate written as adjacent equal pairs | Direct condition "F_B q₃ has q₄=0"; both parities | Same statement for odd B. The gate is equivalent (brute force). Even B is a hand argument in the paper; I proved it in z3. |
| Descent tree (`nine_clock_audit.py`) | 240 roots, 29,534 nodes, Z3+Yices | Exhaustive brute force at small B, plus sampling | The tree logic is sound (§3). Not replayed with a solver here. |
| Fresh test (`separating_route.cpp`) | One branch, 4,000 per base, bases 300035–300039, tolerance 0.035 | 16 branches, 100–200 k per base, 10 other bases, χ² tests | My bases are disjoint from theirs and I use all 16 branches, not one; my test has 1.5 M orbits against their 20,000, about 75× more. |
| Kernels (`nine_lattice.py`, `nine_kernel.py`) | Exact column Hermite reduction | Statistical end-to-end test | The method matches §4.2. Not re-derived symbolically. |

Nothing in the authors' code uses a stronger statement than the paper states. I found no hidden upper bound on B.

## 5. Fresh-base tests (direct column subtraction, `kap9.c`)

### 5.1 The §4.4 distinguishing route

Setup:
- x ranges over the 16 inverse branches of r=(18,14,6,3)/25.
- q = round(Bx) + U[−0.001B, 0.001B]⁴.
- Per-sample predictions use the sample's own y=f(q/B).
- Bases: B=120010+β with 200,000 orbits each, and B=290015+β with 100,000 each. The two families differ in parity and in B mod 3.

Every orbit matched exactly one atom (0 unmatched), and no orbit left the paired region.

| β | n | observed (7/2, 9/2, 50/9, 37/6) | paper | z-scores | χ²(3), p |
|---|---|---|---|---|---|
| 0 | 300,000 | .2288 .1045 .4766 .1901 | .2292 .1042 .4757 .1910 | −0.5 0.6 1.0 −1.2 | 2.20, 0.53 |
| 1 | 300,000 | .2253 .1050 .4895 .1801 | .2266 .1042 .4896 .1797 | −1.7 1.6 −0.1 0.6 | 4.63, 0.20 |
| 2 | 300,000 | .2021 .1045 .5264 .1669 | .2014 .1042 .5269 .1675 | 1.0 0.7 −0.5 −0.9 | 2.04, 0.57 |
| 3 | 300,000 | .1857 .1036 .5408 .1698 | .1858 .1042 .5391 .1710 | −0.1 −1.0 2.0 −1.7 | 5.01, 0.17 |
| 4 | 300,000 | .2023 .1048 .5062 .1867 | .2023 .1042 .5069 .1866 | 0.1 1.1 −0.9 0.1 | 1.56, 0.67 |

Summary statistics:
- All classes: χ²=15.4 on 15 df (p=0.42).
- Homogeneity of the two families (same B mod 5, opposite parity): χ²=17.3 on 15 df (p=0.30), so there is no dependence on B mod 10 or mod 15.
- Per-branch fit (16 branches × 10 bases): χ²=471.7 on 480 df (p=0.60).
- Rate:
  - Largest clock error: 4.74·10⁻⁴ at B≈1.2·10⁵ and 1.99·10⁻⁴ at B≈2.9·10⁵. Both are ≈57/B, a bounded number of steps.
  - Largest terminal error: 6.4·10⁻⁵ and 2.0·10⁻⁵, about 6–8/B.
  - Both are O(B⁻¹) as expected.

### 5.2 Uniform non-constant strings

Setup:
- B=60010…60014 and 60035…60039, with 160,000 strings each.
- Each residue class has two bases of opposite parity and different B mod 3.
- Prediction: μ_β from the supplement, using the general-entry evaluator and the 33 route domains.
  - It was integrated over 48,000 x drawn from ρ₉ (sorted differences of nine uniform digits at base 10⁷).
  - None of the 48,000 draws was rejected; entry lengths were 1–288 steps, and all 33 routes were hit.
- "P45" means P(Z_B = Q(4/5)) under the classifier \|s−4/5\|·B ≤ 3.

| β | B | P45 observed (se .0010) | P45 predicted (±.0003) | E[T/B] observed | E[T/B] predicted (±.0002) | sup over grid \|CDF−pred\| | max dist(Z_B/B, segment)·B |
|---|---|---|---|---|---|---|---|
| 0 | 60010 | .7777 | .7760 | .18724 | .18708 | .0021 | 3.9 |
| 0 | 60035 | .7756 | .7760 | .18721 | .18708 | .0013 | 3.9 |
| 1 | 60011 | .7799 | .7777 | .18694 | .18707 | .0026 | 2.4 |
| 1 | 60036 | .7796 | .7777 | .18718 | .18707 | .0026 | 2.4 |
| 2 | 60012 | .7876 | .7865 | .18732 | .18723 | .0021 | 4.4 |
| 2 | 60037 | .7878 | .7865 | .18720 | .18723 | .0025 | 4.4 |
| 3 | 60013 | .7918 | .7913 | .18724 | .18735 | .0016 | 2.9 |
| 3 | 60038 | .7924 | .7913 | .18737 | .18735 | .0017 | 2.9 |
| 4 | 60014 | .7864 | .7853 | .18716 | .18725 | .0025 | 2.8 |
| 4 | 60039 | .7854 | .7853 | .18740 | .18725 | .0022 | 2.8 |

Findings:
- **Same residue, different parity.** All five pairs agree: \|z\| ≤ 1.4 for P45, and the KS test on T_B/B gives p ≥ 0.31.
- **Differences between classes.** Observed versus predicted differences in P45, relative to β=0:

  | β | observed | predicted |
  |---|---|---|
  | 1 | +.0031 | +.0017 |
  | 2 | +.0110 | +.0105 |
  | 3 | +.0154 | +.0153 |
  | 4 | +.0093 | +.0092 |

  The largest separations are at 15σ, so the five classes are clearly distinct.
- **The small positive offset of observed P45 (≈ +0.001) is explained.** It grows linearly with the classifier tolerance: +0.0006, +0.0011 and +0.0036 at tolerances 1/B, 3/B and 12/B. So it comes from moving terminals that lie within the tolerance of Q(4/5), an O(1/B) effect, and not from the law.
- **Terminal segment.** Every terminal Z_B/B lies within 4.4/B of the segment Q([4/5, 5/6]), as Thm 2 states.

### 5.3 All 165 kernels

Setup:
- For each route, the centre is the ℓ∞ Chebyshev centre of the route domain (LP). It is re-centred because the supplement's centres have small normalised margins.
- The sampling box is certified per sample: r = C⁻¹(f(q/B)) must lie strictly inside the domain.
- All 16 branches are used, with 3,000 orbits per kernel in each of two families: B≈3.0·10⁵ (Bbase ≡ 0 mod 10) and B≈4.5·10⁵ (Bbase ≡ 5 mod 10).
- Atoms that come within 2E of each other anywhere in the box (E=120/B) cannot be resolved at finite B, so they are merged into one group.

| Run | Bases | Orbits | Unmatched / outside box / left paired region | Largest error | Groups tested | χ² (df), p |
|---|---|---|---|---|---|---|
| Family A | 301000–301064 | 495,000 | 0 / 0 / 0 | 47.4 steps | 488 (398 single atoms) | 291.4 (323), 0.90 |
| Family B | 451005–451069 | 495,000 | 0 / 0 / 0 | 47.3 steps | 505 (422 single atoms) | 382.0 (340), 0.058 |
| Pooled | — | 990,000 | — | — | — | 673.4 (663), 0.38 |
| Homogeneity, A vs B | same β, opposite parity | — | — | — | — | 314.0 (311), 0.44 |
| Replication of the 4 family-B outliers: same bases, new seed | 451009–451067 | 48,000 | 0 / 0 / 0 | 43.0 steps | 16 | 12.3 (12), 0.42 |
| Replication of the 4 family-B outliers: B ≥ 10⁶ | 1200019–1200043 | 48,000 | 0 / 0 / 0 | 44.3 steps | 16 | 7.8 (12), 0.80 |

Findings:
- **Every orbit matched an atom of its route's kernel** to within 48 steps (0 unmatched). Every sample stayed in its certified box, and every orbit stayed paired after its first image.
- **The family-B outliers were a statistical fluctuation.** In family B the AA kernels gave χ²=234.2 on 184 df (p=0.007). This came from four groups with \|z\| = 3.2–3.7, in A2/A:A0/AA β=2, A0/A:A0/AA β=4, B1/B:A0/A:A1/AA β=3 and B1/B:A0/A:A0/AA β=3; family A had no group above \|z\|=2.6. With 12,000 new orbits per kernel, the four did not recur, at the same bases or above 10⁶. For example, the deficient group in A2/A:A0/AA β=2 (weight 259/1200 = .2158) came back at .2162 and .2114.
- **A first pass was biased and is superseded.** It discarded samples that were close to two atoms instead of grouping them, which biased near-coincident atoms (χ²=549 on 465 df). The cause is the test, not the kernels; the bias disappears with error-bounded grouping.
- **Resolution.** 398–422 of the 630 atom weights were tested individually. The rest were tested as sums within merged groups, which is unavoidable where two atoms nearly coincide inside the box. The 30 BBB kernels have one atom each.

### 5.4 Uniform O(B) bound

Exact max over all N of T_B, by a memoised orbit graph over all C(B+3,4) difference states:

| B | 40 | 41 | 42 | 44 | 47 | 50 | 64 | 97 | 128 | 150 | 173 | 199 | 223 | 251 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| max T_B | 40 | 44 | 48 | 44 | 53 | 48 | 55 | 73 | 73 | 98 | 96 | 122 | 119 | 145 |
| max T_B / B | 1.00 | 1.07 | 1.14 | 1.00 | 1.13 | 0.96 | 0.86 | 0.75 | 0.57 | 0.65 | 0.55 | 0.61 | 0.53 | 0.58 |

(All of B = 40…50 are in `logs_sims.txt`.) The largest sampled T_B/B at B≈6·10⁴ was 0.51 over 1.6 M strings. These numbers are consistent with Θ(B) and lie far below the proved bound, T_B/B ≲ 1.7·10⁹ (the bound is dominated by 365·C(104,4) = 1,678,315,990).

## 6. Issues

| # | Severity | Location | Evidence | Smallest fix |
|---|---|---|---|---|
| 1 | Minor; must be fixed before release | paper.tex, end of "Verification" (lines 740–742); `RERUN.md`, last line | Both still say that an independent re-derivation with separate top-level encodings "remain pending" or "remains open". | Replace both with a sentence saying that an independent re-derivation of the reductions, with separately written encodings, was completed and agreed with the paper. |
| 2 | Low (scope of an encoding; no effect) | App. C "Cover"; `nine_zero_return_audit.py` `prefix()` | The text claims the cover for every first return. The encoding adds q₁≥1 at each iterate, which drops returns to the zero vector. The gap cannot occur, because F_B q=0 forces q=0 when B≥3. My encoding includes q=0 and is UNSAT for all 16 queries. | Optional: drop the q₁≥1 constraint for the final iterate, or add one clause to the text noting that no nonzero state maps to 0. |
| 3 | Low (robustness of an argument; not an error) | §4.4, "the first column has five distinct entries" | For β=2 and β=4 the first-column entries 29/144 and 233/1152 differ by only 1/1152. Distinctness of these two (X, Z)-marginals therefore rests on two exact rationals and cannot be checked by simulation. The joint laws are separated robustly through the clock: the 50/9 column has 607/1152 and 73/144, a difference of 0.020, confirmed at >20σ (§5.1). | Optional: add that the clock columns also separate every pair, for example 50/9. |
| 4 | Trivial (typos) | §5, first paragraph ("we have … holds"); App. A ("Encode a absolute value") | Text as written. | "At a positive source in nine digits, S_i(…) ≤ S_i(g)+1 for 1≤i≤4"; "an absolute value". |

Observations that need no change:
- The +2 in Thm 1 is never attained in tests; +1 is the largest excess seen, for n ≤ 15.
- The j=4 cover queries are vacuous for ranks 2 and 3.
- The "+7" in the preparation bound is conservative.

## 7. Not checked, or checked only partly

- **Second-solver replays.** cvc5 and Yices are not installed here. Every universal claim I checked used z3 5.4.1 on my own encodings; the authors' cvc5 and Yices runs were not re-executed.
- **Symbolic re-derivation of the route words, releases, terminal offset dynamics and the 165 rational kernels.** These were not redone. They were tested end-to-end by simulation at the resolution given in §5.3. The weights are tested to about ±0.02–0.03 per group at 3σ.
- **The exact 1/1152 first-column gap between β=2 and β=4.** It cannot be resolved by sampling (Issue 3).
- **The descent tree for B>223.** It was not replayed with a solver. For large B the evidence is the reduction re-derivation, 3.1 M sampled sources and the authors' certificate.
- **Analytic steps** (circle lemma, convergence in Thm 1, the limit argument in Thm 4, compact exhaustion and the lattice-sum equidistribution). These were checked by reading and numerical tests, not formally.

## 8. Tools

- macOS 27.2.1 (arm64), Python 3.9.6, z3 5.4.1.0 (z3-solver, commit 145c0c5), NumPy 2.0.2, SciPy 1.13.1 (HiGHS LP), Apple clang 21.0.0 (`-O2`).

Main files in `independent-review/`:
- Model: `km.py`, `zmodel.py`, `kap9.c`.
- SMT: `z_entry.py`, `z_clock.py`, `z_vacuity.py`, `z_vacuity_iv.py`.
- Geometry: `t_faces.py`, `t_exits.py`.
- Brute-force and numerical tests: `t0_*.py`, `t1_majorization.py`, `t_misc.py`, `t_gates.py`, `t_eight.py`, `t_descent_family.py`, `run_descent.sh`.
- Fresh-base runs: `sim_route.py`, `sim_uniform.py`, `sim_allroutes.py`, `analyze_allroutes.py`, `sim_targeted.py`, `pred_global.py`, `compare_uniform.py`, `summarize_route.py`. Each has a lock wrapper `run_*.sh`.
- Raw all-route samples: `allroutes_raw_{300k,450k}.npz` and `targeted_{sameB,bigB}.npz`, with analyses in `allroutes_analysis_families.json` and `targeted_analysis.json`.
- The SMT robustness replay is `z_reseed.py`.
- Logs: `logs_*.txt`.
- `logs_allroutes_150k_ungrouped_superseded.txt` is a first all-route pass. It excluded ambiguous samples instead of grouping them, and so biased the near-coincident atoms (χ²=549 on 465 df). §5.3 replaces it.
- The 115 MB of raw uniform samples were deleted. The seeds are fixed (1000·B+i), so they can be regenerated.
