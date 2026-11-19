# Independent verification of *Linear-time access for eleven-digit Kaprekar dynamics*

Verification run: 18 November 2026. The independent checks were written from the draft before the original code and notes were consulted.
Target: the 18 November 2026 draft (17 pp.).
Code and logs: the sources beside this report and [logs/](logs/).

## Verdict: **PASS WITH ISSUES** (no blocking discrepancy; four non-blocking exposition items)

**Summary**

- I traced each reduction, from finite check to lemma to Proposition 6 to Theorem 1. All the steps hold. I found no gap, quantifier error, bad limit interchange or misused hypothesis.
- My own code and a solver the authors did not use (cvc5) reproduce every decision-relevant constant and identity. These include the 235-state basin, the 391 count, the sharp work 0, 2/7, −3/7 (so 2/7 is attained), the immediate exit at B = 7t, the triangle energy range (5/384, 1/48, 1/21), the paired-return identity, eq. (edge), the one-unit transfer, the ray factor 295/77 and the isolated-cycle predicate. Nothing failed.
- The headline laws hold at 55 bases of both parities and every residue mod 5. That is 461,060 random words, with B from 3.5·10³ to 1.2·10⁸. The largest observed σ_B/B was 0.161, far below 4509, and the largest T_B/B was 0.208. B·‖Z_B − (5,4,3,2,1)/6‖∞ was at most 6, against the proved bound 420. All structural invariants held on every trajectory. These are finite outcomes. They agree with the limit theorems but do not prove them.
- The issues are typesetting and wording only (§4). None of them changes a statement.

Labels: **[me]** = checked by me with my own code or derivation. **[inf]** = inference from reading the text. **[not]** = not checked.

---

## 0. Independence and tools

- I read only the paper until all my checks were frozen (recorded in `logs/`). After that I opened the ancillary README and the `independent/` checker only to run the certificate cross-check (§5). I never needed the authors' sources to diagnose a discrepancy, because there was none.
- Encodings are mine and come from the paper's definitions. F_B is **positional subtraction of the sorted word** (`kap.h`), so eq. (digits) is tested, not assumed. K, ρ (radius map), P_B and Ψ are written from their definitions. Arithmetic is exact integer or `Fraction` throughout. Floats are used only to *report* work values that were found as integer ratios.
- Solvers: **cvc5 1.4.2**, which the authors did not use (Python API), and **Yices 2.7.0** (CLI). No Z3 CLI was available, and I did not use the Z3 Python module, to keep the solvers independent.
- Compute: the fresh-base campaign and the authors' checker ran under an exclusive fcntl compute lock, from separate working copies. Everything else was short.

---

## 1. Reductions (finite check ⇒ lemma ⇒ Proposition 6 ⇒ Theorem 1)

| Step | Re-derivation | Status |
|---|---|---|
| eq. (digits), (track), (outer) | Borrow chain re-derived. Errors are in [−1,0], so after sorting \|F_B−BK\|∞ ≤ 1 and the gap error is ≤ 2. Smallest digit min(q₅−1, B−q₁−[q₁=q₂]). Also the q₅=0 digit list used by the inverse rule. | [me] 14M random states at 7 bases, 0 violations. Max \|F−BK\| = 1, max gap error 2. |
| Energy allowance | e is exactly 2-Lipschitz in gap sup-norm on the simplex (midpoint argument), so the increment is ≤ −δ+4/B. The paper's 20/B is a safe over-allowance. | [me] derivation, plus measured max (E(Fq)−E(Kq))/B = 1.94 < 4 |
| Lemma 2 | From G4, re-read. Majorization tested. | [me] numeric, 0 violations (≈1.4M exact states + 4000 rational orbits). [inf] proof |
| Lemma 3 ⇐ Lemma 7 + Lemma 9 | A non-entering a.e. orbit has an interior limit cycle with e ≥ s*. It must lie in the 5 cells. Ray: F_S re-derived from ρ, and F_S¹⁰ expands by more than 2, so a convergent offset is eventually 0 and lies in finitely many preimages of a line (null). Triangle: I re-derived the cone invariance, the L>0 branch, monotonicity and the "u-difference ≥ 3× after two steps" claim by hand, which gives ≤1 stable point per line, hence Fubini null. Isolated: r₁ is nondecreasing ⇒ a ≤ 0. e ≥ s* plus homogeneity ⇒ predicate 9(4). | [me] all algebra. Predicates 9(1), 9(3) by exact random tests with my own ρ. 9(2), 9(4) by SMT (below). |
| Lemma 4 ⇐ Lemma 8 | One-step increment/B ≤ A_p − (1−θ)J + O(‖x−p‖+1/B). A_p − J ≤ w(p) by torus recurrence and continuity of sorting. Summing h steps with ε ≤ 3‖y‖/7 turns 4/7 into 1/7. | [inf] sound. Constants re-derived. |
| Lemma 5 ⇐ Lemmas 10, 11 | Paired-return identity. The Δ formulas are algebraic consequences of eq. (paired-return): I re-derived both by hand. Greedy prefixes give ⌊(N−1)/19⌋ − 19 = 20 at N = 760, minus 1 for alignment, so 19. Δ-shift monotonicity: guards and H_f = 0 are preserved when w decreases, (ℓ,p,s,r) do not depend on w, and q₄−q₅ grows and q₃−q₄ shrinks with w. | [me] identities. [inf] composition. |
| Prop. 6 | Occupation measures, then K-invariance, then a finite periodic support. Case α<0: Lemma 4 with d/J = (1/7)(√e/√6)/720, so γ = (7√6·720)⁻¹. Case α≥0: Lemma 5 gives c′ ≥ 19/1521 and b′ ≤ 0. The identity e = V(b,c) on two-cycles, with e′ = 2(b−a)b′ − 2(5a+b)c′ and e ≤ 6(4a+b)², gives coefficient 1/(50√6) ≥ γ. Extinction at 2√s*/γ = 10080/√5 < 4508. At the horizon e ≤ 1/2940 < 1/1470, so every gap is > B/7 + margin, which puts the state in S_B. | [inf] every constant re-computed and confirmed (1/1470, 10080/√5 ≈ 4507.87, 19/1521 > 1/100) |
| §6.1 access | Dirichlet(2,…,2) density 11!∏gᵢ re-derived (unit Jacobian, split lengths gᵢ). Portmanteau, shadowing error (2^j−1)/B, 4-Lipschitz energy, then Prop. 6 with T = 4508 + h ≤ B, so 4509B. | [inf] |
| §6.2 pairing / clock | 7m ≥ B+1 ⇒ levels ordered. Re-derived: cross sums of 7 and 5 gaps give margin ≥ m − Σh ≥ 1. R ≤ (B−6)/7 suffices, and at entry R+384 ≤ (B−6)/7. 504 ≤ (B−6)/7 iff B ≥ 3534. The envelope is applied to the pure Ψ-orbit first, then transferred by induction (valid). Stars-and-bars count C(509,5) re-derived. \|X_i\| ≤ 5R ⇒ 420/B. | [me] criterion tested (100k states at the threshold m, 0 failures). [inf] rest |
| Lemma 13 / §7 | Count G₀∏(Gᵢ+1) re-derived. Bound (min G+1)B⁵. Σ over the basin of (min G+1) = 141 + 2·48 + 3·30 + 4·16 = 391. | [me] |

**Result: no gap found in the reductions.** One wording ambiguity (item 4.3) concerns *which* constant has equality at the triangle vertex. My computation resolves it in the paper's favour.

---

## 2. Own encodings of finite obligations (all [me])

| Obligation | Method | Result |
|---|---|---|
| **Four-cycle at B = 5m, complete inverse basin** | Brute force over **all** C(B+4,5) difference states, including zero-innermost ones. Forward F_B orbits with memoization, no depth bound, no inverse rule. | m = 18, 19, 20, 21, 22, 23, 25, 29 (B up to 145, 5.7·10⁸ states): **235 states, layers 4,27,31,38,57,36,25,17, min-gap multiplicities 141/48/30/16, Σ(minG+1) = 391**. Exact. |
| Four-cycle closure and phases | Exact iteration, m = 20…399 and 6 large m up to 10⁹+7 | Closes every time. Min gap ≤ 2 in every phase, q₅ > 0. (It also closes for 3 ≤ m ≤ 19. The basin is 235 already at m = 18, so the m ≥ 20 threshold is conservative.) |
| B ≡ 1, 2 (mod 5) cycles (4 families) | Exact iteration from each threshold to m = 599, plus large m | Periods 8, 8, 10, 8 as stated. Every phase has min gap ≤ 5 (max 5) and q₅ > 0. None of them cycle below the stated thresholds. |
| **Immediate exit, B = 7t** | Exact, t = 400…3000 and 4 large t | Start ∈ S_B. Next gaps are exactly (M, A+1, M−1, M, M, M), so the state leaves S_B. 0 failures. (The exact formula also holds from t = 395. At t = 390 it fails, consistent with t ≥ 400.) |
| **Sharp work orbit (12,8,6,6,3)/14** | w(p) computed from its definition: ΔE_B = ΣG′²−ΣG² in integers, sup over ‖q−Bp‖ ≤ R, one base ~1.4·10⁹ in each of the 14 residues mod 14 | Orbit (12,8,6,6,3) → (12,9,5,2,2) → (12,10,10,7,4). **Phase works 0, 0.285714, −0.428571** (= 0, 2/7, −3/7), the same for every R = 1…6. W = −1/7 = −(2/7)·h·‖y‖∞, with h‖y‖∞ = 1/2. **Equality: 2/7 is attained.** |
| Lemma 8, 4/7 in the ball | Same brute force on 101 distinct low-energy periodic cycles (all-positive gaps, 0<g₀<1/6, e<s*) from lattice denominators 12–42, R = 3, the worst three re-run at R = 5. Also triangle interior and edge points at denominator 224/896. | Max W/(h‖y‖∞) = **−0.889** < −4/7 (stable at R = 5). On the triangle and its edges, values are exactly linear toward V₁ and extrapolate to −2/7 at h‖y‖ = 1/2, i.e. **ratio −4/7: equality in the 4/7 full-closure test at the excluded vertex**. That is consistent with the paper. Ray and isolated-orbit points give −4 and −3. |
| **Catalogue cells, energy ranges** | Exact K on rationals | Ray is fixed by K, e = 6f². Triangle: vertices have period 3, and K, K², K³ are affine on the closed triangle with K³ = id. **Energy inf: interior 5/384 (attained only on the excluded edge V₂V₃, at its midpoint), edges 1/48 (only at V₂/V₃), common sup 1/21 (only at V₁): all unattained.** Isolated orbit: period 3, e = 1/30 exactly. |
| Catalogue sanity (Lemma 7 consequences) | Complete periodic census of K on (1/N)ℤ⁵∩D for N = 6…90, 112, 120, 126, 140 (up to 5·10⁸ points per N) | **All 274** interior periodic points with e ≥ s* lie in the five listed cells (136 ray, 51 triangle interior, 27+27 edges, 33 isolated). **All 1149** periodic points with g₀ ≥ 1/6, e < s* are the center or eq. (two) with b ≥ a ≥ 0, 3a+2b ≤ 1/6. Unexplained: **0**. |
| Two-cycle eq. (two) | Exact K, 1054 rational (a,b) | K swaps the phases everywhere on 3a+2b ≤ 1/6. The displayed adjacent slacks are exact. All 270 in-ball samples with 3a+2b > 1/6 are *not* 2-periodic, so the edge is real. |
| Lemma 3, empirical | 4000 random rational orbits (denominator 1000003) | All of them enter e < s*. Majorization never fails. |
| **Two-edge identities** (`paired.py`) | Integer tests with my own Ψ, P_B, F_B | (A) Ψ = centered P_B: 200k states, 0 mismatches. (B) eq. (paired-return) = Ψ² under the two guards: 185k states, 0. (C) the displayed first output: 0. (D) both gap-increment formulas: 0. **(E) eq. (edge) F_B = P_B + Δe₄ near the other phase on the open edge and at a = 0** (6 edge points, B ≈ 10⁶, offsets up to B/2000): 180k states (86k with Δ>0), 0. Canonical phase F_B = P_B: 180k, 0. **(F) actual canonical return = formula with Δ added to w′ only**: 180k (77k with Δ>0), 0. (G) one-unit transfer identity: 200k, 0. |
| SMT proofs (cvc5 + Yices, QF_LRA) | Own SMT-LIB generator (`smt_checks.py`) | S1 paired-return identity over the **reals** under the guards: unsat/unsat. S2 transfer identity: unsat/unsat. S3 R(ΨX) ≤ R(X)+6 (and +12): unsat/unsat, while the control +5 is sat/sat, so 6 is sharp. S4 ray predicate 9(2): unsat/unsat, and bisection gives the sharp ten-step factor **295/77** exactly (sat at 295/77, unsat 10⁻⁹ below). S5 isolated predicate 9(4), with T_I built from **my own ρ** (not the signed words) at ε = 10⁻⁴ and 10⁻⁶: unsat/unsat. |
| Lemma 10 (integer recurrence) | Exhaustive C, cross-checked against Python on a sub-box | Second assertion: 50 triples × ℓ ∈ [−27,3000] × w ∈ [−3000,2] (4.5·10⁸ pairs), 0 counterexamples. **Composite claim (gain ≥ 1 within 19 guarded returns from p ≤ 2, or guard exit)** on a 5-D box of 2.16·10⁹ states: 0 counterexamples, max 16 returns. The first assertion's box: 419,626 cases, 0 violations. 10⁸ random large-coordinate states in the application regime: 0. |
| Lemma 12 (lower hits) | 2.1·10⁸ random physical X (R from 121 to ~10⁷) plus 63k greedy local searches | 0 failures. **Max first-hit delay 30**, which matches the paper's max tree depth 30 < 32. |

---

## 3. Fresh-base tests of the headline laws (all [me])

Exact word subtraction on uniform nonconstant ordered words (`sim.c`). σ_B is found directly. The first cyclic difference index μ and period λ come from Brent's algorithm. T_B is settled by **comparing actual words** w_μ and w_{μ+λ}. Bases are 3534–3543, 12340–12349, 123450–123459, 1234560–1234569, 12345670–12345679 and 123456780–84, so each decade covers every residue mod 10 (both parities, all classes mod 5) and every class mod 6 and mod 7. Of these, 3534–3539 also appear in an authors' check script (`eleven_centered_clock_checks.py`). The other **49 bases are fresh**.

| Scale | n | σ/B median (range over bases) | max σ/B | max T/B | max B·‖Z−v/6‖∞ |
|---|---|---|---|---|---|
| 3.5·10³ | 200,000 | 0.088–0.091 | 0.161 | 0.208 | 6.00 |
| 1.2·10⁴ | 200,000 | 0.071–0.072 | 0.148 | 0.176 | 6.00 |
| 1.2·10⁵ | 50,000 | 0.065–0.066 | 0.143 | 0.167 | 5.33 |
| 1.2·10⁶ | 10,000 | 0.064–0.068 | 0.142 | 0.166 | 5.00 |
| 1.2·10⁷ | 1,000 | 0.059–0.068 | 0.141 | 0.165 | 5.00 |
| 1.2·10⁸ | 60 | 0.045–0.072 | 0.127 | 0.151 | 4.00 |

Finite outcomes, all 461,060 samples:
- σ_B ≤ 0.162B in every sample, against the 4509B cutoff. σ_B = ∞ in 0 samples, and 0 runs hit the cap. The distribution of σ/B is stable across scales, as expected from tightness and linear time. The 4509 constant is very conservative.
- T_B/B ≤ 0.21. (T_B − σ_B)/B ≤ 0.055. T_B = T_diff + 1 in every sample, which fits T_diff ≤ T_B ≤ T_diff+1.
- B·‖Z_B − (5,4,3,2,1)/6‖∞ ≤ 6 (proved bound 420). Z_B → v/6 at rate O(1/B).
- Along every post-S_B trajectory: F_B = P_B with q₅ > 0 (permanent pairing) at every step. One-step R increase ≤ 6 (allowance 6; recorded 12). Every lower hit came within 12 steps (bound 32). Envelope R ≤ R_σ + 384 held, and R ≤ 126 ≤ 504 after the first visit to R ≤ 120. Cycle states have R ≤ 504. Word and difference iterations were consistent. **0 failures.**
- Terminal periods depend only on B mod 6: {1,4,12}, {24,46,58}, {4,5,12}, {4,8,24}, {1,3,8,50}, {8,20,24} for B ≡ 0…5. This fits the claim that terminal cycles are cycles of the base-independent centered map Ψ in each physical class.
- The exceptional cycles exist in base classes 0, 1 and 2 mod 5 (§2). For the four-cycle in class 0 the bound on its basin mass is 391·11!B⁵/(B¹¹−B), about 1.6·10⁻¹⁴ at B = 10⁴. The simulations met none of the listed cycles.

---

## 4. Findings

**Blocking:** none.

**Non-blocking (exposition only):**
1. **[me] Corrupted macro.** paper.tex lines 210, 214 and 217 contain `N_{` + line break + `m bad}`: the `\r` of `\rm` became a newline. The PDF prints the subscript as italic "mbad" in all three places (confirmed in the PDF content streams). Fix: `N_{\rm bad}`.
2. **[inf] Grammar, §5.2 line 210.** "…and starts in an open set U admit blocks of at most J steps…" is missing a connective. Suggest "…per step, and suppose that starts in an open set U admit blocks…".
3. **[me] Ambiguous equality claim, Appendix B line 321.** "the outside vertex where equality holds for the extended triangle formula" does not say *which* constant. My computation shows that the extended triangle and edge formulas give W = −2/7 at V₁ (h‖y‖ = 1/2), which is **equality in the 4/7 full-closure test**. The actual boundary orbit gives W = −1/7, which is equality for 2/7. Suggest: "where the extended triangle formula attains equality in the 4/7 test".
4. **[inf] Lemma 10 first assertion (line 416).** The case analysis is very compressed ("Thus 1 ≤ z ≤ 4, …"). A referee may ask for it in full. The assertion itself holds in every exhaustive and random test (§2).

Informational (not issues): the threshold m ≥ 20 in Lemma 13 is conservative (the basin is 235 from m = 18). The 20/B energy allowance is generous (4/B suffices). The maximum recorded σ_B/B was about 0.1613, compared with the sufficient cutoff of 4509. The paper already presents 4509B as sufficient only.

---

## 5. Cross-check with the authors' certificate checker (run after my checks were frozen)

I ran `independent/exhaust.py --frozen-certificates` and `audit.py` from a private copy under the compute lock (54 s). Exit codes 0/0. Output: **122,278 rejected systems validated with 116,556 distinct Motzkin certificates**; 791/251/165 faces; 906/172 oriented words; 268 cells (192/6/70); 191 in the ball; 185 ordinary in the ball; sphere suprema 32 open and 1 attained; constants ordinary 2/7, in-ball 4/7. All five regenerated records are **byte-identical** to the frozen ones. Log: [certificate checker](logs/authors_independent_checker.log). These agree with the paper's numbers and with my lattice census.

---

## 6. What I could not check

- **[not] Full exhaustion of the catalogue (Lemma 7)** was outside the independent check. It is covered only by my lattice census (denominators ≤ 140) and by the authors' checker run.
- **[not] The complete Lemma 8 vertex optimization** (748 vertices, 189 cells). I spot-checked it from the definition of w(p) at about 120 periodic orbits, with finite R ≤ 6, one large base per residue class, and limits in R and B approximated rather than taken.
- **[not] The complete 97,465-node centered tree (Lemma 12)**, and the **unbounded** quadrants of Lemma 10 in symbolic form. Both were tested only in finite boxes and by random or adversarial sampling.
- **[not] The six seam predicates of Lemma 11 by SMT.** They were tested numerically on 360k exact states at B ≈ 10⁶.
- **[not] Z3**, since no CLI was available. cvc5 and Yices were used instead.
- **[not] The authors' full `reproduce.sh`.** Only the independent certificate checker was run.
- **[inf]** The analytic steps (occupation measures, Helly selection, BV chain rule, the compression and recession-slope argument of Appendix B, and the proof of Lemma 2 from G4) were checked by reading, not mechanically.
- The simulations give finite-B evidence only. They cannot prove the limit statements.

## 7. Files

- Report: this file.
- Code: `kap.h` (exact F_B, P_B), `selftest.c`, `basin.c`, `cycles.py`, `lattice.c`, `leading.py`, `work.c`, `paired.py`, `smt_checks.py`, `radii.py`, `lemma10.c`, `lowerhit.c`, `sim.c`, `campaign.py`, `summarize.py`
- The raw per-sample archive is omitted. The campaign code and seeds reproduce those records.
- Logs: `logs/` (`sim_summary.md`, `basin_m*.log`, `work_*.log`, `paired.log`, `smt/` with the SMT-LIB files and results, `lemma10_mode*.log`, `lowerhit.log`, `leading.log`, `lattice_*`, `radii.log`, `cycles.log`, `selftest.log`, `authors_independent_checker.log`)
