# Kaprekar dynamics

Research on Kaprekar's routine (repeatedly subtracting the ascending arrangement of a number's digits from the descending one) in arbitrary integer bases.

## Papers

**Cycles of the four-digit Kaprekar map in even bases** (2026), in [`cycles-even-bases/`](cycles-even-bases/) as [PDF](cycles-even-bases/kaprekar-cycles-even-bases.pdf) and [LaTeX source](cycles-even-bases/kaprekar-cycles-even-bases.tex). Archived at Zenodo: [doi:10.5281/zenodo.24325110](https://doi.org/10.5281/zenodo.24325110).

The paper determines every cycle of the four-digit Kaprekar map in every even base B = m·2ⁿ (m odd). For m ≥ 5 the cycles are the lifts of the odd-base cycles, together with one cycle of length n+1 when 2ⁿ ≡ 1 (mod m), and one of length 2(n+1) when m = 5 and n ≡ 2 (mod 4). The case m = 5 recovers Hasse and Prichett's classification. The cases m = 1 and m = 3 are given completely. The even bases with exactly one cycle are 3·2^(2k+1), 5·2^(2k+1) and 7·2ⁿ with 3 ∤ n. The classification is formally verified in Lean 4.

## Reproducing the results

```sh
cd cycles-even-bases/ancillary
bash reproduce.sh            # numerical checks (Python 3 and a C++17 compiler)
cd lean && lake build        # Lean 4.28.0 proof of the classification
lake env lean TheoremC.lean  # prints the axioms used
```

See [`cycles-even-bases/ancillary/README.md`](cycles-even-bases/ancillary/README.md) for expected outputs and versions.

## Priority record

[`COMMITMENTS.txt`](COMMITMENTS.txt) lists SHA-256 digests of research materials that are not yet released.

## License

Text, figures and data: CC BY 4.0. Code: MIT. See [`LICENSE.md`](LICENSE.md).

## Contact

rustem@galyavin.com
