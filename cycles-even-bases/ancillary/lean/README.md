# Theorem C in Lean 4

Contact: rustem@galyavin.com

Theorem C classifies the nonconstant terminal cycles of four-digit Kaprekar subtraction in every even base. Its corollary determines exactly which even bases have one such cycle and gives its least period. Four-digit words may begin with zero. The proof and its documentation are covered by the bundle's [LICENSE.md](../LICENSE.md).

## Build and check the axioms

Use elan with Lean 4.28.0. From this directory, run the following commands. Elan reads the version from `lean-toolchain`. The proof uses only the libraries distributed with Lean.

```sh
lake build
lake env lean TheoremC.lean
```

On 17 November 2026, a fresh copy with no Lake build products and only elan's pinned Lean toolchain completed the build in 16.09 s on macOS arm64. Lake's final line was

```text
Build completed successfully (16 jobs).
```

The second command took 0.74 s and printed the axiom dependencies of both entry points.

```text
'KaprekarBoundary.every_even_base' depends on axioms: [propext, Classical.choice, Quot.sound]
'KaprekarBoundary.theoremC' depends on axioms: [propext, Classical.choice, Quot.sound]
```

These are Lean's propositional extensionality, classical choice and quotient soundness. The [per-declaration axiom output](formal_verification.txt) covers the intermediate results as well.

## Correspondence with the paper

In [TheoremC.lean](TheoremC.lean), `every_even_base` takes any even integer base `B ≥ 2` and supplies its factorization `B = m·2ⁿ`, with positive odd `m` and `n ≥ 1`. The factorized entry point `theoremC` proves `TheoremCStatement m n`. A cycle is the orbit of a periodic word, and its length is the least positive number of steps that returns that word to itself.

| Field of `TheoremCStatement` | Statement in the paper |
|---|---|
| `sorting` | Sorting places the four digits in descending order and preserves the multiplicity of every digit. |
| `subtraction` | One step gives a valid four-digit word whose value is the descending arrangement minus the ascending arrangement. |
| `odd_core` | For odd `m ≥ 5`, the lifted odd-core cycles retain their periods. There is one additional cycle of length `n+1` if `2ⁿ ≡ 1 (mod m)`, or one of length `2(n+1)` if `m = 5` and `n ≡ 2 (mod 4)`. The conditions are disjoint, each additional cycle is distinct from every lift, and these cases exhaust the periodic words. |
| `core_three` | For `m = 3`, odd `n` gives one cycle of length `3(n+1)`. Even `n` gives two cycles, of lengths `n+1` and `2(n+1)`. |
| `binary`, `binary_fixed` | In base two the only nonconstant periodic words are `0111` and `1001`. Both are fixed. |
| `pure_power` | For `m = 1` and odd `n ≥ 3`, there are `(n−1)/2` cycles of each length `n` and `n+2`. Even `n ≥ 2` gives `(n−2)/2` cycles of each of those lengths, plus one of length `n/2` and one of length `n/2+1`. |
| `uniqueness` | Exactly one nonconstant cycle occurs for `m = 3` or `5` with odd `n`, or for `m = 7` when three does not divide `n`. |
| `unique_lengths` | In the core-three and core-five cases, write `n = 2k+1`. The unique cycles then have lengths `6(k+1)` and `1`, respectively. In the core-seven case the length is `3`. |

Theorems A, B and D and the fixed-point basin limits are outside this formalization.
