Independent re-derivation of the reductions from the text, with separately written encodings

**Verdict: PASS WITH ISSUES.** Every computer-assisted step in Appendix A.1–A.5 was re-derived from the paper's text and re-proved with separately written code. Every result agrees with the paper. Each step's statement is exactly what its proof uses. No upper bound on B is hidden, and no case is missing. The only issue is one low-severity editorial discrepancy, between the appendix's description of one query and the bundled query (Issue 1). It does not affect any result. The headline laws hold at fresh bases.

Paper reviewed: the 18 November 2026 draft of *Seven-digit Kaprekar dynamics: global law, two clocks, parity*. The encodings and logs accompany this report. The release incorporates the wording correction in §3.

## 1. Method

1. I read the paper and restated each finite assertion with its quantifiers (table in §2). For each one, I traced where the proof uses it and checked the margin arithmetic at every point where a check feeds a proposition.
2. I wrote all encodings before opening the original reproduction files. The transition model was derived by hand from column subtraction. It was checked against integer subtraction of L_B(q) on every state for B ≤ 40 and on 20,000 random states with B < 10^7 (`model.py`). In C, it was checked on every state for 18 ≤ B ≤ 200 (`exhaust.c`, `int128mismatch=0`). My sorting networks differ from the authors': a 12-comparator 6-input network from Dobbelaere's list, and a 21-comparator 7-input odd–even transposition network. Both were validated on all permutations and all 0-1 inputs.
3. Only after that did I compare with the authors' encodings (§4).
4. Fresh-base tests use a separate C program. It applies κ_B to whole words by sorting and column subtraction with borrow, and finds T_B with Brent's cycle algorithm on words. It never uses difference vectors.

## 2. Computer-assisted steps

"UNSAT" means that the counterexample query (hypotheses plus the negated conclusion) has no model, so the implication holds. "SAT" marks a deliberate sanity query. The model of every such query was replayed exactly.

| # | Location | Precise statement (all quantifiers) | Proof use; scope check | My check | Result |
|---|---|---|---|---|---|
| 1 | Lemma A.1 (i)–(vi), 8 claims | For all reals 0 ≤ c ≤ b ≤ a < 1, with q_j = H^j(a,b,c): (i) a, a₁ < 1/2 ⇒ 2a₁ + (q₁)₂ ≥ 1; a < 1/2 and 2a + b ≥ 1 ⇒ q₁ = (a + (−s)₊, a − s₊, 1 − a − b). (ii) a ≥ 1/2 ⇒ W(q₃) ≤ W(q₀). (iii) a > 3/5 ⇒ a₃ = a₂; if also a₆ = a, then W(q₆) ≤ max(0, W₀ − min(5a−3, 1−a)). (iv) a > 1/2, a₆ < 1, W₆ = W₀ > 0 ⇒ q₀ ∈ Z. (v) a > 1/2, W₀ < ½ min(2a−1, 1−a) ⇒ W₆ = 0. (vi) 1/2 < a ≤ 3/5, (b,c) ≠ (2−3a, 2a−1), max(\|b−2+3a\|, \|c−2a+1\|) < (2a−1)/100 ⇒ some j ≤ 4 has a_j > 3/5 or W_j < 2a_j − 1 | Prop. 3. Every use meets the stated hypotheses: (iv) at accumulation points with a = A ∈ (1/2, 3/5], so a₆ = A < 1; (vi) near the first point of Z, with crossing 3/5 excluded; (ii) iterated with a ≥ 1/2 | Own QF_LRA encoding (median as max(min, min(max, ·)); a₃ ≠ a₂ negated directly). Run in Z3 and in Yices on the exported SMT-LIB | All 8 UNSAT in both solvers. Four false variants SAT, with exact rational replays: strict three-step decrease, (iv) without excluding Z, (vi) with only j ≤ 1, (v) with a 4× threshold |
| 2 | Lemma A.2, first part | For β = 0, 1, every node reachable from \|m\|, \|S\| ≤ 20 (m ≡ β): for all B ≡ β (mod 2) with B ≥ 100000 and all integers u with 3B/4 + 128 ≤ u ≤ 5B/6 − 128, the representative three-step rule gives the exact offsets and u-shift, and every partial step moves u by ≤ 20. 884 and 845 nodes, the listed cycles, transients ≤ 10 and 11 blocks, cumulative shift ≤ 24 and 30 | Prop. 4. The margins check out: block starts stay ≥ 3B/4 + 1000 − 30 − 20 and ≤ 5B/6 − 1000 + 50; in the final circuits they stay ≥ H_B − 60; skipped states have u − 3B/4 ≥ 176 ≥ 128 | Own symbolic block composition in (u, b, 1) (`edge.py`). At B = 100000 + β, u = 80000, the digits are sorted by value. Each adjacent comparison (non-strict suffices, because only the sorted values enter the differences), each w ≥ 1 and each state inequality is certified exactly on the real polyhedron of the strip (both vertices plus both recession rays). Offsets and shifts are checked to be constant forms | Zero guard failures. 884/845 nodes. Cycles, ℓ, D and r exactly as in Prop. 4's table. Maxima: transient 10 and 11, shift 24 and 30, partial shift 20 |
| 3 | Lemma A.2, second part | For each cycle phase, each B mod 4, and u_end = ⌈3B/4⌉ + 256 − j (1 ≤ j ≤ D): the path to first repetition is the (3,2,1)(B − B₀)/4 translate of the path at B₀ = 100000 + (B mod 4); terminal cycles have u − 3B/4 < 128; the discrepancy 3k + re + h_end − r(u_end − 3B/4) lies in [−43/2, 343/8] | Prop. 4. The discrepancy interval is ⊂ [−200, 200] | Own code (`edge_terminal.py`). For every path state: w ≥ 1; slopes non-increasing in the (value, slope)-lexicographic order, so the order persists for all real t ≥ 0; symmetric slope differences equal (3,2,1). Spot-checked by direct simulation at t = 1, 7 and 250001 | 140 families (52 even, 88 odd). All certificates hold. Discrepancy range exactly [−43/2, 343/8] (even bases [−19, 131/4]) |
| 4 | Lemma A.3, return | For all integers with B ≥ 100000, u ≤ B − 1000, w ≥ 1000, m ≥ −6, S, S−m, E, E−m ≥ 16, J, K₁, L ≥ 8: three steps give (u−3, w−1, \|m+1\| − 5) | Prop. 5, for every j < N. I re-derived that the source margins imply these guards: S_j, E_j ≥ 16; S_j − m_j, E_j − m_j constant or ≥ S_j; J_j, L_j non-decreasing; (K₁)_j − K₁ ≥ 3j; w_N > B/6 | Unbounded-integer Z3 query with the general step (cases w > 0, v > 0, both zero), so intermediate positivity is not assumed | UNSAT. Sanity SAT: guards feasible at m = −6 and at m ≥ 6; replacing J ≥ 8 by J ≥ 0 breaks the rule |
| 5 | Lemma A.3, exit | B ≥ 100000, −d ≥ 8000, B + 6d ≥ 24000, E = e ∈ {12, …, 15}, m ∈ {−5, …, −1} with m ≡ e: the next state has u − 1 and offsets (eq. neg), and lies in Prop. 4's domain | Prop. 5, negative side | 10 Z3 queries; the conclusion includes \|m'\|, \|S'\| ≤ 20 and 3B/4 + 1000 ≤ u' ≤ 5B/6 − 1000 | 10 UNSAT. All 10 hypothesis sets SAT (not vacuous) |
| 6 | Lemma A.3, vertices | The closed canonical sign regions have vertices {A, U, G, D, V} (+) and {A, G, F, V} (−) | §5 "Which sources are reached" and the support argument | Exact vertex enumeration from (a) C(K) + canonical rotation + sign and (b) the homogeneous source inequalities + ordered coordinates + sign | Both descriptions give exactly the listed vertices. d ∈ [0, 1/3] (+) and [−1/6, 0] (−). Clocks at the vertices lie in [0, 1/4] |
| 7 | Lemma A.3, residue tables | Fast columns as printed in Prop. 5 | Phase weights p_β | Own computation: N, s/e and m_N from (S or E, m) mod 8, then (eq. neg), then the rate from my edge graph. Checked against several representatives per residue | All 8 × 4 printed cells match, covering 96 residue pairs. Fast fractions 16/32 and 12/32 (+), 12/16 and 6/16 (−) |
| 8 | Lemma A.3, integrals / Cor. 6 | M_Q = 1339759/43200000, M_P = 937297/5400000 | Cor. 6 | Own decomposition, without the authors' tetrahedra: P is split into the 48 sign/radius-order pieces of 𝒟; Q comes directly from C⁻¹. Exact SymPy integration over triangulated polytopes. Sanity: ∫_𝒟 ρ = 1 | Both masses match exactly. Monte Carlo (2·10⁶ points): 0.03101 and 0.17363 |
| 9 | Lemma A.4 (i) | For all integers B ≥ 18: q, Fq, F²q positive ⇒ F³q positive | §6: a positive start either stays positive or reaches zero within two steps | Z3, unbounded B | UNSAT. Without B ≥ 18: SAT at B = 4, replayed |
| 10 | Lemma A.4 (ii) | B = 2b + β, b ≥ 9, q = (u, v, 0), v > 0: Fq positive; F²q has w = 0 iff q is in the listed set, and then F²q = (b+1, 1, 0); F²q positive ⇒ F³q positive | §6: every start is permanently positive by step 6 | Z3: five queries per parity (a, both directions of the iff, the value, c), plus "(b+1,1,0) is not listed for b ≥ 9" | All 12 UNSAT |
| 11 | Lemma A.4 (iii) | F(u, 0, 0) = (max(B−u, u−1), min(B−u, u−1), 0) | §6 | Z3 for all B ≥ 2 | UNSAT |
| 12 | Lemma A.5, one step | For all integers B ≥ 2 and every positive q: V(Fq) ≤ V(q) + 8 | Thm. 2 | Z3 with both networks | UNSAT (both). V(Fq) > V(q) + 7: SAT at B = 5, replayed |
| 13 | Lemma A.5, nine steps | For all integers B ≥ 2: q, …, F⁹q positive and V(q) > 200 ⇒ V(F⁹q) < V(q) | Thm. 2 (54B; the ≤ 272 cap on cycles) | Unpartitioned Z3 query, plus the auxiliary bounds V_j ≤ V_{j−1} + 8 and V_j ≥ V₀ − 8(9−j). These are sound consequences of #12 and the negated goal | Z3: UNSAT in 116 s for the whole unpartitioned query. Yices on my encoding is only partial: of the 16 source-sign partitions, 8 are UNSAT and 8 hit the time limit (600 s, then 450 s with implied ordering facts added); none is SAT. The unpartitioned Yices run timed out at 900 s. Supported by the exhaustive and sampled checks below |

The constant in Theorem 2 is 54B + 614125 + 6 + 1, as the text says. The bound B ≥ 400 is never used: every step needs only B ≥ 18. This is an observation, not an error.

**Supporting exhaustive checks** (`exhaust.c`, every state at every base listed). At all 283 bases 18 ≤ B ≤ 300, and at 129, 130, 400, 401, 500–507, 1000 and 1001: zero failures of A.4 (i), (ii) and (iii), of the one-step bound, or of the nine-step decrease (largest nine-step change −8). The last zero-face visit is at step ≤ 3, inside the claimed "by step six". A stratified sample of 40,000 states at B up to 10⁹ (near the attractor with V just above 200, boundary faces, and gaps near B/4) also gave zero failures.

**Residue-table check.** The working notes, read after all checks were finished, record a correction: in the negative odd table, row E ≡ 3 (mod 8), column 3 is fast. My independently computed table agrees with that cell and with the paper's printed table.

**Reproduced §7 claims.** 22,315,002 nonconstant difference vectors in bases 129, 130, 400, 401. Maximal first cyclic index h: 56, 49, 113, 127.

## 3. Issues

| # | Severity | Location | Evidence | Smallest fix |
|---|---|---|---|---|
| 1 | Low (editorial) | Appendix, text after Lemma A.4: "It checks the formula in (iii) through all bases." | In the bundled query (`rank_one_formula` in `seven_zero_transience_audit.py`) the hypothesis is `B>=18`. The lemma itself is stated only for B ≥ 18, so the proof is unaffected. My own query proves (iii) for all B ≥ 2 (UNSAT) | Write "for every B ≥ 18", or remove `B>=18` from that query |

There are no other discrepancies. Each statement's quantifiers match its use. No proposition relies on a case that the check omits, and no check assumes a bound on B that the theorem does not state.

## 4. Comparison with the authors' encodings (after my checks)

- **Transition model.** Both use the same hand-derived borrowed digits. For positive steps, both remove the maximal digit B−1 and sort the other six; the image is (B−1−d₅, d₀−d₄, d₁−d₃) in both. The networks differ. My zero-face steps sort all seven digits; theirs remove one copy of B−1. Same model.
- **Domains.** The potential, zero-face, entry and edge domains are identical to mine. I matched the facet `RETURN_DOMAIN`, `NEGATIVE_EXIT_DOMAIN` and `SOURCE_*` row by row against Lemma A.3, eq. (source) and the sign bounds. The only mismatch is the rank-one query's B ≥ 18 (Issue 1). The authors' potential prefixes also assert the ordering of each image state and nq₀ < B. Both follow from sorting, so they do not narrow the domain; my queries omit them and are still UNSAT.
- **Tie handling.** In edge and facet words, the authors require a margin of one at a reversed-index tie. This is stronger than necessary; I use non-strict comparisons, which suffice. In terminal lifts they order ties by index, and I order them by (value, slope). Both conventions certify all 140 families.
- **Bounds and constants.** Identical: B ≥ 100000; strip margins 128; threshold ⌈3B/4⌉ + 256; source margins 1000/8000/24000; the discrepancy interval; V > 200 and +8. Their 328 and 330 distinct edge guards (counted in a scratch copy) are consistent with the text.
- **Method differences.** For the facet guards, the authors replay 490 Farkas identities over ℚ; I use unbounded-integer SMT. For the masses, they integrate nine tetrahedra and their 72 preimages; I integrate a direct decomposition of 𝒟. The independent methods agree.

## 5. Fresh-base tests (bases not used in the paper)

Limit law ν_β, from 2·10⁷ samples ξ ∼ ρ through my own implementation of Theorem 1's construction. One sample did not enter C(K) within the 10⁴-step cap: a slow plateau at a ≈ 0.99994. Exact rational iteration showed that it enters at step 10212.

| β | P(fast) | mean | var | E e^{4t} | min t_f | max t_s |
|---|---|---|---|---|---|---|
| 0 (even) | 0.56436 | 0.140588 | 0.003080 | 1.79771 | 0.0003 | 0.24998 |
| 1 (odd) | 0.375 | 0.142617 | 0.003089 | 1.81236 | | |

Predicted mean gap (odd − even): 0.002029. In every sample, 0 < t_f < t_s < 1/4.

**A. Two clocks and parity weights.** I drew points ξ ∼ ρ with \|x − 3z\| ≥ 0.03. For each point, I took the 512 words whose difference vectors fill q₀ + [0,8)³, so that every residue mod 8 occurs once, and computed T_B for every word by word-level subtraction.

| B | B mod 8 | boxes | fast fraction per box (side: value × count) | max \|T − B t_atom\| | min clock gap ·B | times far from both atoms |
|---|---|---|---|---|---|---|
| 200000 | 0 | 100 | −: 0.75 ×21, +: 0.5 ×79 | 59.3 | 789 | 0 |
| 200001 | 1 | 100 | −: 0.375 ×15, +: 0.375 ×85 | 76.0 | 801 | 0 |
| 200002 | 2 | 100 | −: 0.75 ×15, +: 0.5 ×85 | 73.1 | 763 | 0 |
| 200003 | 3 | 100 | −: 0.375 ×17, +: 0.375 ×83 | 83.7 | 777 | 0 |
| 200004 | 4 | 100 | −: 0.75 ×26, +: 0.5 ×74 | 62.6 | 804 | 0 |
| 200005 | 5 | 100 | −: 0.375 ×17, +: 0.375 ×83 | 67.0 | 778 | 0 |
| 200006 | 6 | 100 | −: 0.75 ×15, +: 0.5 ×85 | 71.0 | 777 | 0 |
| 200007 | 7 | 100 | −: 0.375 ×13, +: 0.375 ×87 | 50.2 | 769 | 0 |
| 1000003 | 3 | 30 | −: 0.375 ×8, +: 0.375 ×22 | 98.8 | 4267 | 0 |
| 1000004 | 4 | 30 | −: 0.75 ×3, +: 0.5 ×27 | 72.9 | 3987 | 0 |

In all 860 boxes, every time lies within 99 steps of one of the two predicted atoms. The fast share equals p_β exactly: 1/2 or 3/4 in even bases, 3/8 in odd bases, the same for all eight residues mod 8.

**B. Population law, uniform random words.**

| B | N | mean T/B | limit | var | limit | E e^{4T/B} | limit | sup grid \|ΔF\| | max T/B |
|---|---|---|---|---|---|---|---|---|---|
| 20011 | 10⁵ | 0.14294 | 0.14262 | 0.003075 | 0.003089 | 1.8145 | 1.8124 | 0.0033 | 0.2519 |
| 20012 | 10⁵ | 0.14064 | 0.14059 | 0.003075 | 0.003080 | 1.7981 | 1.7977 | 0.0020 | 0.2518 |
| 20013 | 10⁵ | 0.14284 | 0.14262 | 0.003095 | 0.003089 | 1.8141 | 1.8124 | 0.0038 | 0.2511 |
| 20014 | 10⁵ | 0.14109 | 0.14059 | 0.003074 | 0.003080 | 1.8012 | 1.7977 | 0.0043 | 0.2513 |
| 100005 | 2·10⁴ | 0.14274 | 0.14262 | 0.003077 | 0.003089 | 1.8131 | 1.8124 | 0.0040 | 0.2485 |
| 100006 | 2·10⁴ | 0.14074 | 0.14059 | 0.003088 | 0.003080 | 1.7989 | 1.7977 | 0.0046 | 0.2495 |

The standard error of the mean is 1.8·10⁻⁴ for N = 10⁵ and 3.9·10⁻⁴ for N = 2·10⁴. The observed odd − even gap, about 0.0020, matches the prediction.

**C. Exact weighted law over all B⁷ − B words** (exhaustive difference-level DP with exact word counts; the totals equal B⁷ − B).

| B | 500 | 501 | 502 | 503 | 504 | 505 | 506 | 507 | 1000 | 1001 |
|---|---|---|---|---|---|---|---|---|---|---|
| mean T/B | 0.1466 | 0.1574 | 0.1567 | 0.1577 | 0.1466 | 0.1573 | 0.1566 | 0.1576 | 0.1435 | 0.1499 |
| mean − limit | 0.0061 | 0.0148 | 0.0161 | 0.0151 | 0.0060 | 0.0147 | 0.0160 | 0.0150 | 0.0029 | 0.0072 |

The bias is about c/B, with c ≈ 3 for B ≡ 0 (mod 4) and c ≈ 7.5 otherwise. It halves from B ≈ 500 to B ≈ 1000, and at B ≈ 2·10⁴ it is consistent with Table B.

**D. Destination and linear bound** (every difference vector, every base listed).

| B | 400 | 401 | 500–507 | 1000 | 1001 |
|---|---|---|---|---|---|
| max T_B | 114 | 128 | 141–161 | 272 | 282 |
| max over cyclic states of ‖q − (3B/4, B/2, B/4)‖∞ | 3.0 | 2.5 | ≤ 3.5 | 3.0 | 2.5 |

The bound 54B + 614132 holds with a very wide margin. Max T_B/B is 0.27–0.32 at B ≈ 500–1000 and 0.249–0.252 in the random samples at B ≥ 2·10⁴, consistent with Θ(B) and liminf ≥ 1/4. Every cycle lies within 3.5 ≤ 34 of (3B/4, B/2, B/4), so the bound 34/B holds; the same holds at all bases 18–300 (maximum 3.5).

**E. Propositions 4 and 5 at fresh bases** (`prop_tests.py`).

- Prop. 4: 324 states, with B random in [10⁵, 4·10⁵) plus every corner offset at both ends of u. The observed h − r(u − 3B/4) range is [−18.5, 42.875], inside Lemma A.2's [−21.5, 42.875]; the upper end is attained.
- Prop. 5: 300 random source states. Max \|T − B t_r\| is 23.25 (positive side) and 20.5 (negative side), against the claimed 216. The other clock was always more than 1000 steps away, so the residue rule's choice of rate is observable and correct.

## 6. Tools

Python 3.9.6, z3 5.4.1 (Python API), Yices 2.7.0 (`yices-smt2`), SymPy 1.14.0, NumPy 2.0.2, Apple clang 21.0.0 (C, `-O2`), macOS on a 10-core machine with 16 GB of memory. cvc5 is not installed, so the paper's cvc5 runs were not repeated; Yices was the second solver for A.1 and A.5.

## 7. Independent review files

| File | What it does |
|---|---|
| `model.py` | Reference transition model and its self-checks |
| `exhaust.c` | Exhaustive lemma checks, maxima and the exact weighted law |
| `traj.c` | Word-level tests: clocks, population, limit law |
| `smt_common.py` | Shared SMT encoding (sorting networks, steps, V) |
| `smt_entry.py` | Lemma A.1 |
| `smt_zero_clock1.py` | Lemma A.4; Lemma A.5 one-step bound |
| `smt_clock9.py` | Lemma A.5 nine-step decrease |
| `smt_facet.py` | Lemma A.3: return and exits |
| `edge.py`, `edge_terminal.py` | Lemma A.2 |
| `facet_geometry.py` | Vertices, residue tables, masses |
| `prop_tests.py` | Propositions 4–5 at fresh bases |
| `sample_lemmas.py` | Stratified lemma sampling at large bases |
| `analyze_fresh.py` | Summary of the fresh-base campaign |
| `logs_*.txt`, `logs/` | Computation logs and summary statistics. The six raw population samples are regenerable as described in the bundle README |
| `smt2/` | SMT-LIB exports sent to Yices |

## 8. Not checked

- cvc5 (not installed).
- A complete second-solver proof of the A.5 nine-step claim on my own encoding. It rests on Z3 (whole query UNSAT), 8 of 16 Yices partitions, and the exhaustive and sampled checks. The authors report all 16 partitions UNSAT in both Z3 and Yices on their encoding; I did not re-run that.
- The authors' 490 and 338 Farkas-identity replays. By design, I did not run the authors' scripts; I replaced those replays with my own unbounded SMT proofs.
- The measure-theoretic parts were checked by reading, not by computation: the Riemann sums, the exhaustion argument and absolute continuity.
