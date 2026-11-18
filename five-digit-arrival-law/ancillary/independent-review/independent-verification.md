# Independent verification of P2: "A common arrival-time law for five-digit Kaprekar dynamics"

- **Scope:** Independent re-derivation of the reductions from the text, with separately written encodings.
- **Date:** 17 November 2026.
- **Object:** `docs/drafts/p2-five-digit/paper.md` (776 lines, dated 17 November 2026).
- **Method.** I reconstructed every computer-assisted step from the paper's text alone. I wrote and ran my own Z3 encodings and brute-force programs. Only then did I open the project's encodings (`ancillary/source/five_excursions.py`, `five_clock_audit.py`, `five_joint_audit.py`, `five_quantitative_audit.py`, and the integration routine in `five_limit_law.py`) to compare them. I did not modify the paper or any existing file. All artefacts are in this directory. Every computation longer than a few seconds ran under the exclusive `fcntl` lock through `lockrun.py`.
- **Versions.** Z3 5.4.1 (the `z3-solver` package in `kaprekar-investigation/.venv-research`), Python 3.9.6, numpy 2.0.2, and Apple clang 21.0.0 (`-O2`).

## Summary

- **Verdict: PASS WITH ISSUES.** I found no mathematical error, no false statement and no gap in any reduction. Every quantitative claim I tested holds, at every base I tested.
- **The computer checks are sound.** Each computer-assisted step, as the paper states it, really gives what the theorems use. My own encodings use a different transition model (a total sorting network instead of 26 guarded branches). All 43 claims return UNSAT for every integer B ≥ 30, with no upper bound on B. Probes show that the encodings are not vacuous.
- **The predictions hold at fresh bases.** The exact finite-base laws at B = 997, 1500, 2401, 2999, 3000, 5003, 10007 and 10009 converge to the paper's F at the rate O(1/B). The 1,000,000-sample Monte Carlo runs agree with the exact values within 2 standard errors. The cycle weights 2/3 and 1/3 hold to within 10⁻⁵ at B = 997, and to within 10⁻⁷ at B = 10009. The maximum stays inside the paper's bounds on every base 30–3000 and at the large bases. The two-cycle basin formula (7) is exact on every even base 8–36 by enumerating all strings, and on every even base 30–3000 by the pair graph. The overshoot law (12) matches at B up to 2,560,009.
- **Issues (none of them changes a theorem):**
  - **I1 (required before release).** The "Verification scope" section still says this verification is pending.
  - **I2 (minor).** The project's audits assume, without a solver check, that the 26 guarded branches partition the state domain. If they did not, two of the audits could be UNSAT for no reason. I checked that the partition holds.
  - **I3 (editorial).** Several symbols carry two meanings.

---

## 1. Computer-assisted steps, restated precisely

Notation follows the paper. *Valid pair:* integers 1 ≤ u ≤ B−1, 0 ≤ v ≤ u. f_B is the pair map: sort the digits of (B⁴−1)u + (B³−B)v and take (d₁−d₅, d₂−d₄). The margins are a = u+v−B and b = 2u−v−B. ρ is the base-free margin rule (20):

ρ(a,b) = (b, a−3) if a ≥ 1; (b, −2a) if a ∈ {−1, 0}; (b+1, 3) if a = −2.

The following sets are as in the paper: H_B (19), J_B (30), T_{2,B} = {(B/2, B/2−1), (B/2+1, 0)} for even B (empty for odd B), and R_B = T_{2,B} ∪ {(B/2, 0)}.

| # | Where in the paper | Precise statement (all quantifiers over integers) |
|---|---|---|
| CA1 | §Core, (20), "substitution shows each image remains in H_B"; App. A item 1 | ∀B ≥ 30 ∀ valid (u,v) ∈ H_B: f_B(u,v) is valid, f_B(u,v) ∈ H_B, and (a′,b′) = ρ(a,b). In (u,v) form, the image is (u−1, u−v+1), (B−v, 2u−B) or (B−v, u−v−1). The paper also uses: on H_B, 2u ≥ B, v > 0 and max(a,b) ≤ (B−3)/2. |
| CA2 | §Global entry, first paragraph; App. A item 1 | ∀B ≥ 30 ∀ valid (u,v) with v ≥ 1: (u′,v′) = f_B(u,v) is valid, 2u′ ≥ B−1 and 2(u′−v′) ≤ B+2. |
| CA3 | Finite affine entry lemma; App. A item 2 | ∀B ≥ 30 ∀ valid q₀ ∈ J_B ∃j ∈ {0,…,6}: f_B^j(q₀) ∈ H_B ∪ T_{2,B}. |
| CA4 | (21); App. A item 3 | ∀a,b ≥ −2 with max(a,b) ≥ 10: max(ρ²(a,b)) = max(a,b) − 3 + 𝟏{min(a,b) = −2}. |
| CA5 | (22) and the cycle list; App. A "finite terminal calculation" | For each of the 144 starts (a,b) ∈ [−2,9]²: ρ²(a,b) ∈ [−2,9]²; every coordinate on the orbit is ≤ 10; the first-repetition index h satisfies max h = 10; the cycles are exactly the four listed; and min and max of h − (2/3)max(a,b) are −3 and 16/3. |
| CA6 | §Exact two-cycle basin; App. A item 4 | ∀ even B ≥ 30 ∀ valid q ∉ R_B: f_B(q) ∉ R_B. (The forward closure r₁ → r₂ → r₁ and r₃ → r₁ is proved by hand in the text.) |
| CA7 | (31); App. A item 4 | ∀B ≥ 30 with B ≡ 1 (mod 3), ∀ valid q ∈ H_B: u+v ≡ 0 ⇔ u′+v′ ≡ 0 (mod 3). |
| CA8 | (33); App. A (enumeration) | Each of the four linear parts of (13) has determinant ±1 and permutes (ℤ/3)². |
| CA9 | (36); App. A item 5 | ∀B ≥ 30 ∀ valid q₀ ∈ J_B with 4u₀ ≥ 3B, ∀j ∈ {0,…,6}: if q_i ∉ H_B for all i < j and q_j ∈ H_B, then 2max(a_j,b_j) ≤ 2u₀ − B + 40. |
| CA10 | (2), (5); App. A "Exact integration" | F(1/25) = 713903/39062500. With cutoff N = 300, the exact rational enclosure of E[L], built from (25)–(29), lies strictly inside (0.1827011465, 0.1827011669). |

The appendix also lists finite-base trajectory comparisons, maximum witnesses and an enumeration of input weights at bases 2–8. These support the reductions but no theorem rests on them. The paper says so correctly, so I treat them as corroboration and not as steps of a proof.

## 2. Reduction verdicts

For each step, I checked whether the check, if true, gives exactly what the theorem uses. I looked for missing cases, quantifier order, small bases, boundary cases and off-by-one errors in the clocks.

| # | Used for | Verdict | Notes |
|---|---|---|---|
| CA1 | (20) on H_B, hence (23) | **Complete** | The sorting argument in the text is correct; I re-derived the three cases by hand. The theorem needs only invariance plus ρ on H_B, because orbits started in H_B never leave it. The base-free box enumeration is therefore valid for H_B starts, even though some box points are not in H_B at B = 30 (for example (9, −2) needs B ≥ 47). Covering extra box points does no harm. The threshold B ≥ 30 is sufficient; my probe shows that CA1 already holds for B ≥ 24 and fails at B = 20. |
| CA2 + high-region lemmas | Strip reach in O(√B), Theorems 2 and 4 | **Complete** | The text proves some of these analytically: the high branches (u−1, u+v−B+1) and (u, 2u−B); "a diagonal arises only from (B−1, B−2)"; the bound v ≤ B−1−k(k−1)/2; and the zero-coordinate rule v′ = min(u−1, B−u). I re-derived them and also checked them in Z3 (A2″, A3–A3⁗). The step count ⌈√B⌉+2 off-diagonal + 2 diagonal + 3 preparatory + 6 strip steps = ⌈√B⌉+13 is correct (the off-diagonal count is in fact ≤ ⌈√B⌉+1). |
| CA3 | (3), (9) | **Complete, and tight** | "Within six steps" covers j = 0,…,6; both encodings model it this way. T_{2,B} is needed: without it the check fails, for example at B = 30, q₀ = (15, 14). Six steps are needed: five fail, for example at B = 31, q₀ = (16, 0). States visited before H_B ∪ T₂ cannot be cyclic, because H_B and T₂ are forward invariant. So the time to reach H_B ∪ T₂ is a valid prefix of the clock. |
| CA4 + the "at most once" argument | (23) | **Complete** | Combining CA4 with CA5 gives h − (2/3)M₀ ∈ [−3, 16/3 + 2/3] = [−3, 6], exactly (23). The text proves "the exceptional block occurs at most once" by hand. I encoded it as an induction: the block-start minimum lies in {0, 3} after an exceptional block, and this is preserved while M ≥ 10 (A5″, A5‴). I also checked the condition behind "cannot skip first capture": the odd-time state of a block with M ≥ 10 has max > 3, so it is not cyclic (A5⁗). |
| CA5 | (22), cycle labels | **Complete** | The paper leaves one point implicit: every listed cycle state lies in H_B for all B ≥ 30. I checked it (A10). |
| CA6 | (7)–(8) | **Complete** | With the hand-proved forward closure (A6′), the basin is exactly {r₁, r₂, r₃}. The histogram (7) then follows from (18); I confirmed it by enumerating all strings. |
| CA7 | (31), Theorem 3 | **Complete** | The text also gives a hand proof via ρ. |
| CA8 | (33), the 1/3 : 2/3 split | **Complete** | Determinant ±1 means q ↦ Aq + (Bb+c) is a bijection mod 3, so exactly 3 of the 9 residue classes satisfy (31) at the endpoint. With (32), each class carries 1/9 of the mass. |
| CA9 | Upper bound in (9) | **Complete, very slack** | The use is correct. Core entry from {4u₀ ≥ 3B} ∩ J_B cannot meet T₂, because R_B has u ≤ B/2+1. My probe shows the claim even holds with +0 in place of +20, and first entry occurs only at j ≤ 3. |
| CA10 | (2), (5) | **Complete** | The termination fact (28) is correct; I re-derived it, including A + C = 1 − 3δ. The omitted-mass polynomial (25) is correct; I re-derived it by inclusion–exclusion on Dirichlet(2,2,2). The width formula for the first moment is correct. |

**The analytic reductions around the checks.** I re-derived the rest of each theorem's argument by hand. This covers the borrowed digits (17) (also checked in Z3); the weights (18) (also checked by full enumeration for B = 4–36); the relation T − h ∈ {0, 1} with the cycle-input correction; the leading map (13); the strip and range invariance (24); the residue transport (31)–(33); the step accounting in the proof of Theorem 4 (constants 92/3, 51 and the witness bound −4); and, for Theorem 5, (40)–(46), including e_B ∈ [−3, 7] for even B and [−10/3, 22/3] for odd B, inside the claimed |e_B| ≤ 8. I found no errors. The rate bound (10) is analytic and deliberately loose. I checked its constant chain (27000, 2112, 4224 against 10⁵) and found it consistent.

## 3. Encodings, results, and comparison with the project

### 3.1 Reviewer files (in this directory)

| File | What it does | Result |
|---|---|---|
| `z3_independent.py` → `z3_independent.log` | 43 checks: the identity (17) and its digit ranges; CA1 (four parts); CA2 (three parts); the high-region lemmas (five); CA3; CA4 plus the "at most once" induction and the non-cyclic odd-time state; CA6 plus forward closure plus R_B ∩ H_B = ∅; CA7 (both directions); CA9 (j = 0…6); the overshoot exit states (43) and (B−R, H−ν) ∈ H_B; CA5 (enumeration); CA8; and "every cycle state lies in H_B for B ≥ 30". **Transition model:** the five borrowed digits pass through a bubble-sort comparator network (a total function, checked by the 0–1 principle). The image is (s₁−s₅, s₂−s₄). I validated it against plain numerical subtraction at 1,275 concrete points. No upper bound on B. | **43/43 pass**, total run time 2.6 s |
| `z3_independent.py --probe` → `z3_probes.log` | Non-vacuity and mutation probes. | Every hypothesis set is satisfiable. Weakened claims are SAT, so the encodings do constrain: 5 steps fail; dropping T₂ fails; {r₁, r₂} alone is not predecessor-closed. The threshold for CA1 is between B = 20 (SAT) and B = 24 (UNSAT). (36) holds with +0. First entry from {4u₀ ≥ 3B} happens only at j ≤ 3. (The probe labels the j-feasibility lines "UNEXPECTED"; that is a cosmetic flaw in the probe, not a failure. Those lines are informational.) |
| `strings_exhaustive.c` → `strings_exhaustive.log` | All B⁵ strings for B = 4…36, by the definition only (desc − asc), on the full functional graph. Checks (18) and lists every cycle with its arrival histogram. | (18) is exact for all B = 4…36. Formula (7) matches exactly for **every even B = 8…36**, including B = 10, 12, 16, 20, 24, 30. It fails at B = 4 and 6, which the paper does not claim. No period-2 cycle occurs in any odd base. |
| `pairgraph.c`, `sweep_pairgraph.py` → `sweep_30_1200.log`, `sweep_1201_3000.log` | Exact pair graph. The pair map is computed from the integer (B⁴−1)u + (B³−B)v, not from (17). The program computes h, the cycles and their labels, and the weighted T distribution. On **every base B = 30…3000** it tests (3), (9), (22)–(23) on H_B, the six-step lemma on J_B, the global count ⌈√B⌉+13, (36), the witness (37), the labels by B mod 3, and (7). | **No failure on any of the 2,971 bases.** 3h − 2M on H_B lies in [−9, 16] (claim [−9, 18]). J_B entry takes at most 6 steps (attained). Global entry has slack ≥ 9. The (36) excess is ≤ −1 (claim ≤ 40). max T − B/3 − √B/3 lies in [−0.83, 7.48] (claim [−5, 60]). |
| `mc_fresh.c` → `mc_fresh.log` | Monte Carlo with 10⁶ samples per base, at B = 997, 1500, 2401, 2999 and 3000. It uses **only the numerical Kaprekar map**: Floyd's algorithm gives the exact tail length (= T_B) and the cycle reached. No pair coordinates are used. | See §4 |
| `limit_exact.py` → `limit_exact.log` | My own exact rational evaluation of F and E[L] from (13)–(16) and (25)–(29): half-plane clipping and triangle-fan integration with ∫ sⁱtʲ = i! j!/(i+j+2)!. It also checks that the retained mass equals 1 − R(ε) exactly. | **F(1/25) = 713903/39062500, exactly as in (2).** F(1/10) = 51653/400000. F(1/5) = 31645397267/54000000000 ≈ 0.5860258753. F(3/10) ≈ 0.9683938275. E[L] ∈ [0.182701146507, 0.182701166802] at N = 300, **strictly inside (5)**. |
| `tau_mc.py` → `tau_mc.log` | Plain Monte Carlo of τ, 2·10⁷ samples. | Agrees with `limit_exact.py` within 1.6 standard errors at all four t values. |
| `overshoot_mc.c`, `overshoot_fast.c`, `overshoot_limit.py` → `overshoot.log`, `overshoot_large.log` | Theorem 5 at fresh bases (see §4.5). The fast variant gives byte-identical output to the slow one at B = 10007, 40009 and 160001. | (12) and (44) hold |
| `compare_branches.py` → `compare_branches.log` | Comparison with the project's transition encoding (run after all of the above). | See §3.2 |

From the bundle root, run `python3 independent-review/z3_independent.py [--probe]`. Compile the C files with `cc -O2 -o <bin> <file>.c -lm`. The logs record run parameters (the overshoot logs omit the seeds), and the drivers are `sweep_pairgraph.py`, `summarize_fresh.py` and `lockrun.py`.

### 3.2 Comparison with the project's encodings

**What is the same.** The project encodes the same ten statements as CA1–CA10, with the same quantifiers and the same domain B ≥ 30, unbounded. The six-step and (36) formulas match mine clause for clause. The T₂ encoding is equivalent: (2u = B ∧ 2v = B−2) ∨ (2u = B+2 ∧ v = 0) forces B to be even. The terminal enumeration and the determinant check are the same computation.

**What differs:**

1. **Transition model.** The project uses 26 guarded affine branches, each requiring "guards ⇒ image". I use one total comparator network. `compare_branches.py` shows that the two models are the same relation for every B ≥ 2: the guards cover every valid pair, the 325 pairs of branches are pairwise disjoint, and each branch's image equals the network's image on that branch (UNSAT on all 26).
2. **No coverage check in the project's audits.** The project's six-step and (36) audits model one step as a disjunction over branches, and its other audits check one branch at a time. None of them asserts that the branches cover the domain. If a state were covered by no branch, the disjunction would drop it, and no single-branch check would ever look at it. The checks could then be UNSAT for no reason. My total model has no such hole, and the partition check above closes it for the project's model. See issue I2.
3. **Extra checks in my encoding.** I also encoded in Z3 several lemmas the paper proves by hand: the high-region branches, the diagonal lemma, the zero-coordinate rule, max(a,b) ≤ (B−3)/2 on H_B, "at most once", the cycle states lying in H_B, (17), (43), and the exit state lying in H_B. All hold.
4. **Exact integration.** The project integrates with Green's theorem on polygon edges; I use triangle fans, an independent implementation. Both follow the cell decomposition stated in the paper's text. The exact values agree.

## 4. Fresh-base results

### 4.1 Distribution: exact finite-base value, Monte Carlo, and the limit F (paper's definition (13)–(16))

The limit values: F(1/25) = 0.0182759, F(1/10) = 0.1291325, F(1/5) = 0.5860259, F(3/10) = 0.9683938. Columns: "exact" is P(T_B ≤ tB) over all B⁵−B inputs (pair graph plus (18)). "MC" uses 10⁶ samples (numerical map only), with its standard error.

| B | B mod 3 | t | exact | exact − F | B·(exact − F) | MC ± se | (MC − exact)/se |
|---|---|---|---|---|---|---|---|
| 997 | 1 | 1/25 | 0.0159213 | −2.35e−3 | −2.35 | 0.016122 ± 0.000126 | +1.59 |
| 997 | 1 | 1/10 | 0.1224529 | −6.68e−3 | −6.66 | 0.122763 ± 0.000328 | +0.95 |
| 997 | 1 | 1/5 | 0.5747297 | −1.13e−2 | −11.26 | 0.575070 ± 0.000494 | +0.69 |
| 997 | 1 | 0.3 | 0.9641901 | −4.20e−3 | −4.19 | 0.964457 ± 0.000185 | +1.44 |
| 1500 | 0 | 1/25 | 0.0167746 | −1.50e−3 | −2.25 | 0.016809 ± 0.000129 | +0.27 |
| 1500 | 0 | 1/10 | 0.1248375 | −4.29e−3 | −6.44 | 0.124825 ± 0.000331 | −0.04 |
| 1500 | 0 | 1/5 | 0.5775774 | −8.45e−3 | −12.67 | 0.577844 ± 0.000494 | +0.54 |
| 1500 | 0 | 0.3 | 0.9646748 | −3.72e−3 | −5.58 | 0.964693 ± 0.000185 | +0.10 |
| 2401 | 1 | 1/25 | 0.0175302 | −7.46e−4 | −1.79 | 0.017669 ± 0.000132 | +1.05 |
| 2401 | 1 | 1/10 | 0.1268435 | −2.29e−3 | −5.50 | 0.127129 ± 0.000333 | +0.86 |
| 2401 | 1 | 1/5 | 0.5816440 | −4.38e−3 | −10.52 | 0.582539 ± 0.000493 | +1.82 |
| 2401 | 1 | 0.3 | 0.9664543 | −1.94e−3 | −4.66 | 0.966665 ± 0.000180 | +1.17 |
| 2999 | 2 | 1/25 | 0.0179396 | −3.36e−4 | −1.01 | 0.017874 ± 0.000132 | −0.50 |
| 2999 | 2 | 1/10 | 0.1282356 | −8.97e−4 | −2.69 | 0.128488 ± 0.000335 | +0.75 |
| 2999 | 2 | 1/5 | 0.5848104 | −1.22e−3 | −3.65 | 0.584343 ± 0.000493 | −0.95 |
| 2999 | 2 | 0.3 | 0.9677676 | −6.26e−4 | −1.88 | 0.967899 ± 0.000176 | +0.75 |
| 3000 | 0 | 1/25 | 0.0175261 | −7.50e−4 | −2.25 | 0.017483 ± 0.000131 | −0.33 |
| 3000 | 0 | 1/10 | 0.1269836 | −2.15e−3 | −6.45 | 0.127032 ± 0.000333 | +0.15 |
| 3000 | 0 | 1/5 | 0.5818382 | −4.19e−3 | −12.56 | 0.582048 ± 0.000493 | +0.43 |
| 3000 | 0 | 0.3 | 0.9665797 | −1.81e−3 | −5.44 | 0.966421 ± 0.000180 | −0.88 |
| 5003 | 2 | 1/25 | 0.0183080 | +3.21e−5 | +0.16 | – | – |
| 5003 | 2 | 1/5 | 0.5856711 | −3.55e−4 | −1.77 | – | – |
| 10007 | 2 | 1/25 | 0.0182788 | +2.9e−6 | +0.03 | – | – |
| 10007 | 2 | 1/10 | 0.1290297 | −1.03e−4 | −1.03 | – | – |
| 10007 | 2 | 1/5 | 0.5858770 | −1.49e−4 | −1.49 | – | – |
| 10007 | 2 | 0.3 | 0.9683326 | −6.1e−5 | −0.61 | – | – |
| 10009 | 1 | 1/25 | 0.0180598 | −2.16e−4 | −2.16 | – | – |
| 10009 | 1 | 1/10 | 0.1283550 | −7.77e−4 | −7.78 | – | – |
| 10009 | 1 | 1/5 | 0.5846863 | −1.34e−3 | −13.41 | – | – |
| 10009 | 1 | 0.3 | 0.9678515 | −5.42e−4 | −5.43 | – | – |

**Reading.**
- The Monte Carlo runs agree with the exact values; the largest gap is 1.8 standard errors.
- The exact error exact − F is of order 1/B. Within one class of B mod 3 the scaled error is stable: for t = 1/5 it is about −12.6 at B = 1500 and 3000 (class 0), and between −10.5 and −13.4 at B = 997, 2401 and 10009 (class 1). This is consistent with an O(1) clock offset that depends on the terminal cycle. It is not a discrepancy in F, and the theorem asserts only the limit.
- The mean E[T/B] at B = 10007 is 0.182740, against E[L] = 0.1827012.

The full table, including the extra t values, is in `fresh_table.md`.

### 4.2 Maximum

| B | exact max T | T of witness (37) | B/3 + √B/3 | max − (B/3 + √B/3) | paper (9): [B/3+√B/3−5, B/3+√B/3+60] | MC max (10⁶) |
|---|---|---|---|---|---|---|
| 997 | 345 | 345 | 342.86 | +2.14 | holds | 343 |
| 1500 | 518 | 517 | 512.91 | +5.09 | holds | 514 |
| 2401 | 819 | 818 | 816.67 | +2.33 | holds | 816 |
| 2999 | 1018 | 1017 | 1017.92 | +0.08 | holds | 1010 |
| 3000 | 1024 | 1017 | 1018.26 | +5.74 | holds | 1017 |
| 5003 | 1691 | 1689 | 1691.24 | −0.24 | holds | – |
| 10007 | 3369 | 3368 | 3369.01 | −0.01 | holds | – |
| 10009 | 3372 | 3371 | 3369.68 | +2.32 | holds | – |

On all bases 30–3000, the residual max T − B/3 − √B/3 lies in [−0.83, 7.48]. Every witness meets the paper's lower bound B/3 + √B/3 − 4, and comes within 0–7 of the true maximum. The uniform bound (3) holds on every base 30–3000. The data agree with a second-order coefficient of exactly 1/3. The stated constants −5 and +60 are correct but conservative.

### 4.3 Cycle weights and the joint law

| B | B mod 3 | exact masses | MC (10⁶) | conditional P(T ≤ B/5 \| label), exact |
|---|---|---|---|---|
| 997 | 1 | T₃ 0.6666567, T₄ 0.3333433 | period 3: 0.666813 ± 0.00047; period 4: 0.333187 ± 0.00047 | T₃ 0.56809, T₄ 0.58801 |
| 2401 | 1 | T₃ 0.6666654, T₄ 0.3333346 | 0.666145 ± 0.00047; 0.333855 ± 0.00047 | T₃ 0.57866, T₄ 0.58762 |
| 10009 | 1 | T₃ 0.66666659, T₄ 0.33333341 | – | T₃ 0.58402, T₄ 0.58602 (F = 0.58603) |
| 1500, 3000 | 0 | T₁ 1 − 1.18e−8 and 1 − 1.48e−9; the rest is T₂ | period 1 only | – |
| 2999, 5003, 10007 | 2 | T₅ 1 | period 4 only | – |

- The weights 2/3 and 1/3 are confirmed.
- The gap between the conditional laws of T₃ and T₄ shrinks like 1/B: 0.020 at B = 997, 0.009 at 2401, 0.002 at 10009. This confirms asymptotic independence.
- The T₂ masses at B = 1500 and 3000 equal (8) exactly. For example, at B = 3000 the count is 359,759,990 = 40B² − 80B − 10.
- On every base 30–3000, the set of cycle labels is exactly the one predicted by B mod 3, plus T₂ if and only if B is even.

### 4.4 Rare two-cycle basin (7), by exhaustive enumeration

I enumerated all strings (`strings_exhaustive.log`). The observed counts (N₀, N₁, N₂) equal (2, 35B²−75B−12, 5B²−5B), with no other arrival times, for every even B = 8, 10, …, 36. Examples: B = 10 gives (2, 2738, 450); B = 12 gives (2, 4128, 660); B = 16 gives (2, 7748, 1200); B = 20 gives (2, 12488, 1900); B = 24 gives (2, 18348, 2760); B = 30 gives (2, 29238, 4350). The formula fails at B = 4 and 6, where the basin is larger; this is outside the claimed range B ≥ 30. The pair-graph sweep confirms that the basin is exactly 3 pairs, and confirms (7), on every even B = 30…3000.

### 4.5 Overshoot law (Theorem 5)

I sampled starts weighted by (18) in two bands: E₁ = {λ ∈ [0.5, 1.5], z ∈ [0.6, 0.9]} with M = 6.48, and E₂ = {λ ∈ [1, 2.5], z ∈ [0.55, 0.95]} with M = 21.945. The limit values come from (12) by quadrature (`overshoot_limit.py`).

| Band | B (samples) | B·P(E_B) (limit M) | E[Y \| E] ± se (limit) | P(Y ≤ y₁ \| E) ± se (limit) | P(Y ≤ y₂ \| E) ± se (limit) | max \|e_B\| |
|---|---|---|---|---|---|---|
| E₁, y₁ = −1, y₂ = −0.75 | 10007 (2·10⁴) | 6.2124 (6.48) | −0.87296 ± 0.00190 (−0.86916) | 0.33360 ± 0.00333 (0.33439) | 0.66105 ± 0.00335 (0.65969) | 5 |
| | 40009 (2·10⁴) | 6.3479 | −0.86064 ± 0.00191 | 0.32345 ± 0.00331 | 0.65060 ± 0.00337 | 3.67 |
| | 160001 (10⁵) | 6.4148 | −0.87131 ± 0.00085 | 0.33622 ± 0.00149 | 0.66238 ± 0.00150 | 5 |
| | 640007 (10⁵) | 6.4474 | −0.86950 ± 0.00085 | 0.33388 ± 0.00149 | 0.65982 ± 0.00150 | 5 |
| | 2560009 (5·10⁴) | 6.4637 | −0.86735 ± 0.00121 | 0.33094 ± 0.00210 | 0.65934 ± 0.00212 | 3.67 |
| E₂, y₁ = −1.25, y₂ = −1 | 10007 (2·10⁴) | 20.233 (21.945) | −1.53225 ± 0.00286 (−1.52838) | 0.72750 ± 0.00315 (0.72756) | 0.89115 ± 0.00220 (0.88871) | 5 |
| | 40009 (2·10⁴) | 21.093 | −1.51920 ± 0.00285 | 0.72220 ± 0.00317 | 0.88500 ± 0.00226 | 3.67 |
| | 160001 (10⁵) | 21.521 | −1.52966 ± 0.00127 | 0.72921 ± 0.00141 | 0.88893 ± 0.00099 | 5 |
| | 640007 (10⁵) | 21.733 | −1.52757 ± 0.00127 | 0.72620 ± 0.00141 | 0.89086 ± 0.00099 | 5 |
| | 2560009 (5·10⁴) | 21.839 | −1.52655 ± 0.00181 | 0.72508 ± 0.00200 | 0.88772 ± 0.00141 | 3.67 |

The statistics converge to (12) at the stated rate O(B^{−1/2}):
- At B = 40009 the largest deviation is 0.0085 in E[Y | E₁], which is 1.7/√B.
- At B = 640007 and 2,560,009, all 12 entries are within 2.2 standard errors of the limit, and 10 of them are within 1.5.
- A 2·10⁴-sample run at B = 160001 had shown a 2.8-standard-error low value of P(Y ≤ −1.25 | E₂) (in `overshoot.log`). The 10⁵-sample rerun shows +1.2 standard errors, so that value was a sampling fluctuation.

- **Uniform formula (44).** Over every sample, max |e_B| = 5 (claim ≤ 8). There are no violations of (41) and no high-branch mismatches.
- **Band probability.** At B = 10007, 40009 and 160001, the scaled error √B·(B·P(E_B) − M) is −26.8, −26.4, −26.1 for E₁ and −171.3, −170.5, −169.5 for E₂. So the correction is O(B^{−1/2}), which is Pr(E_B) = M/B + O(B^{−3/2}) as stated. Most of it comes from the factor u − v = B − r − v in (47).

## 5. Overall verdict and issues

**PASS WITH ISSUES.** No theorem, constant, formula or reduction is wrong. Every computer-assisted step is logically complete and reproduces with independent encodings. Every headline prediction holds at fresh bases: the limit law, the moment bounds, the maximum asymptotics, the joint law and its cycle weights, the rare basin, and the overshoot law. The issues below are about the record and readability, not the mathematics.

**I1 — The verification section is out of date (required edit before release; blocks release only as a process step).**
- *Location:* §"Verification scope", last paragraph, lines 629–632.
- *What is wrong:* the paragraph says that the independent re-derivation and the fresh-base validation "remain pending". That section should record the completed checks.
- *Fix:* replace the paragraph with a short summary of this report. Suggested content: "An independent review re-derived each reduction from the text and wrote separate Z3 encodings using a total sorting-network transition model. All 43 statements were UNSAT for every B ≥ 30. It checked that the 26 branch guards used here partition the domain. It confirmed (2) and (5) by an independent exact integration. It tested the predictions at fresh bases B = 997–2,560,009: the CDF, the maximum, the cycle weights, (7) for even B = 8–36 by full enumeration, and the overshoot law. No discrepancy was found." Point to `independent-review/independent-verification.md`.
- *Clarity point:* the same section could say explicitly that Theorems 1 and 3 already rely on App. A items 1, 3 and 4 (the residue part) and on the terminal enumeration. Each of these also has a hand proof in the text. At present only the theorems that need the six-step and (36) implications are tied to the solver.

**I2 — The branch partition is assumed but not checked in the solver (minor; no effect on any result).**
- *Location:* App. A, "Exact transition encoding", the sentence "These guards … partition the positive-coordinate domain"; and `five_clock_audit.py` and `five_quantitative_audit.py`.
- *What is wrong:* the six-step and (36) audits model one step as a disjunction over the 26 guarded branches. A state covered by no branch would leave the orbit, and the formula would become UNSAT for no reason. The audits never assert coverage or disjointness; they rely on the combinatorial remark.
- *Severity:* low. I proved coverage and disjointness for all B ≥ 2, and equality with a total sorting map (`compare_branches.log`). So the published UNSAT results stand.
- *Fix:* add two solver assertions to `five_clock_audit.py`: "valid ∧ ¬(∨ guards)" is UNSAT, and the guards are pairwise disjoint. Then change the appendix sentence to "…partition the domain (checked by the solver)".

**I3 — Symbols with two meanings (editorial, low).**
- In (33), "q_end = Aq + Bb + c" uses A for a matrix, while A is the margin in (15) and (24); B there means the base times a vector b. Suggest writing P q + B·β + γ.
- M means max(a,b) in (21)–(23) and the band mass in Theorem 5.
- δ means min(y, 1−x) in (14) and the 0-or-½ offset in "H = B/2 + δ" in §"The conditional overshoot law".
- R means R(ε) in (25), r + n in (41), and the set R_B in App. A.
- H_B is the core, while H = ⌊(B+1)/2⌋ in Theorem 5.
- *Fix:* rename locally. No mathematical change.

**Observations, not issues.**
1. Several constants are much looser than needed: (36) holds with +0; the clock error in (23) observed on H_B is [−3, 16/3]; and the residual of the maximum lies in [−0.83, 7.48] on bases 30–3000. Tightening is optional.
2. The claims for B ≥ 30 also hold somewhat below that threshold: (7) holds for every even B ≥ 8, and CA1 holds for B ≥ 24. The paper does not claim these, so nothing needs to change.
3. I did not check the bibliography, Prichett's labels ("These agree with Prichett's labels [1, Tables 4.1–4.3]") or the disclosure section.
