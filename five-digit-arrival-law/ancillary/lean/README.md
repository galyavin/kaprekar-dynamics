# Five-digit Kaprekar dynamics in Lean

The theorems below formalize the finite lemmas in *A common arrival-time law for five-digit Kaprekar dynamics*. Each pair-map theorem quantifies over integer bases with no upper bound. Five-slot words allow leading zeroes.

| Paper statement | Lean theorem in namespace `FiveDigit` |
|---|---|
| Borrowed digits (17) and the difference-pair transition in Appendix A | `sort_ordered` and `sort_multiset` prove that the sorting network orders all five digits and preserves each digit's multiplicity. `subtraction_value` identifies the positional result with descending minus ascending, and `subtraction_valid` proves that every output digit lies in range. `borrowed_pair` identifies the next pair with sorting the borrowed digits. `constant_iff` characterizes zero outer difference. For every nonconstant valid word and every iterate in base at least two, `iterated_subtraction` equates its difference pair with the corresponding pair-map iterate. |
| Appendix A, claim 1. Core (19), rule (20), first-step strip inequalities | `core_pair_rule` gives all three coordinate rules. `core_rule` proves forward invariance and the margin recurrence for every base at least 30. `strip_bounds` gives validity and both first-step inequalities whenever the second coordinate is positive. |
| Appendix A, claim 2. Six-step entry from (30) | `six_step_entry` gives an index from zero through six at which the pair lies in the core or on the additional even-base two-cycle. The two-cycle predicate uses doubled coordinates, so it is empty in odd bases. |
| Appendix A, claim 3. Two-step descent (21) | `two_step_descent` proves the exact change in the maximum of two integer margins at least −2 when their maximum is at least 10. |
| Appendix A, claim 4. Rare predecessors (34) and residue rule used in (31) | `rare_predecessors` proves backward closure for every even base at least 30. `rare_images` gives all three forward images. `core_residue` proves that divisibility of the coordinate sum by three is preserved in both directions on the core for bases congruent to one modulo three. |
| Appendix A, claim 5. First-entry bound (36) | `first_entry_bound` bounds both doubled margins by `2*u₀ − B + 40` for every first core-entry index from zero through six, when the starting pair is in the strip and `4*u₀ ≥ 3*B`. `high_first_entry` also proves that such a first entry exists. |
| Terminal enumeration (22), all 144 starts | `box_invariance` gives two-step invariance of [−2,9]² and an intermediate maximum at most 10. `listed_cycles` checks all transitions on the four displayed cycles and that their vertices are distinct. `box_bounds` bounds the first visit time and its scaled error. `box_transient` proves that this time is the first visit to any periodic state. `box_extrema` attains time 10 and the upper error 16/3 at (5,7), and the lower error −3 at (9,4). |

The six-step proof splits the sorting orders that can occur on legal inputs, follows each affine image, and closes the remaining inequalities with Lean's integer-arithmetic tactic. It checks the entry-time and margin conclusions on the same paths. The box proof enumerates two indices in `Fin 12` by kernel reduction. The global probability laws, asymptotic estimates, and basin multiplicities are proved separately in the paper.

## Reproduction

Run these commands in this directory with Lean 4.28.0 and its bundled Lake. The project uses only the bundled standard library.

```sh
lake build
lake env leanchecker --fresh FiveDigit
```

The build prints the axiom dependencies of the main theorems; the second command replays the compiled declarations in a fresh kernel environment.

On an Apple M5, a clean copy built all nine jobs without warnings in 20.36 seconds. Fresh kernel replay completed in 66.21 seconds. The 26 axiom listings used only standard Lean axioms. [Verification record](verification.txt)
