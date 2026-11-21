# Kaprekar dynamics

Research on Kaprekar's routine (repeatedly subtracting the ascending arrangement of a number's digits from the descending one) in arbitrary integer bases.

Author's website: [galyavin.com](https://galyavin.com/).

## Papers

**Cycles of the four-digit Kaprekar map in even bases** (2026), in [`cycles-even-bases/`](cycles-even-bases/) as [PDF](cycles-even-bases/kaprekar-cycles-even-bases.pdf) and [LaTeX source](cycles-even-bases/kaprekar-cycles-even-bases.tex). Archived at Zenodo: [doi:10.5281/zenodo.24325110](https://doi.org/10.5281/zenodo.24325110); revised 18 November 2026 to credit earlier work (Zenodo [doi:10.5281/zenodo.24343078](https://doi.org/10.5281/zenodo.24343078)). Paper record: [doi:10.5281/zenodo.24358739](https://doi.org/10.5281/zenodo.24358739).

The paper determines every cycle of the four-digit Kaprekar map in every even base B = m·2ⁿ (m odd). For m ≥ 5 the cycles are the lifts of the odd-base cycles, together with one cycle of length n+1 when 2ⁿ ≡ 1 (mod m), and one of length 2(n+1) when m = 5 and n ≡ 2 (mod 4). The case m = 5 recovers Hasse and Prichett's classification. The cases m = 1 and m = 3 are given completely. The even bases with exactly one cycle are 3·2^(2k+1), 5·2^(2k+1) and 7·2ⁿ with 3 ∤ n. The classification is formally verified in Lean 4.

**Arrival times for the four-digit Kaprekar map in even bases** (2026), in [`arrival-times-even-bases/`](arrival-times-even-bases/) as [PDF](arrival-times-even-bases/kaprekar-arrival-times-even-bases.pdf) and [LaTeX source](arrival-times-even-bases/kaprekar-arrival-times-even-bases.tex). Archived at Zenodo: [doi:10.5281/zenodo.24329321](https://doi.org/10.5281/zenodo.24329321). Paper record: [doi:10.5281/zenodo.24358744](https://doi.org/10.5281/zenodo.24358744).

For every even base B = m·2ⁿ with m ≥ 5 odd, the paper gives the exact distribution of the number of steps needed to reach any union of the cycles inherited from base m. As n grows, the number of steps divided by n converges to a law with at most five atoms. In the bases 5·2ⁿ the proportion of inputs reaching the fixed point tends to 1033/1125 for n ≡ 0 and to 94/125 for n ≡ 2 (mod 4), against the values 3/4 and 15/16 conjectured by Devlin and Zeng. As m grows, the number of steps to reach any cycle obeys a total-variation law that is uniform in n.

**A common arrival-time law for five-digit Kaprekar dynamics** (2026), in [`five-digit-arrival-law/`](five-digit-arrival-law/) as [PDF](five-digit-arrival-law/kaprekar-five-digit-arrival-law.pdf) and [LaTeX source](five-digit-arrival-law/kaprekar-five-digit-arrival-law.tex). Archived at Zenodo: [doi:10.5281/zenodo.24330159](https://doi.org/10.5281/zenodo.24330159). Paper record: [doi:10.5281/zenodo.24358746](https://doi.org/10.5281/zenodo.24358746).

With five digit slots and a growing base B, the number of steps needed to reach a cycle, divided by B, converges through all integer bases to one explicit continuous law. Along each residue class of B modulo 3 the limiting time is independent of the terminal cycle. The paper also determines the exact basin of the even-base two-cycle, proves that the maximum arrival time is B/3 + √B/3 + O(1), and identifies the overshoot behind the square-root term. Its finite integer-affine checks are verified by SMT encodings and in Lean 4.

**Majorization in odd-digit Kaprekar dynamics and the nine-digit arrival law** (2026), in [`odd-digit-majorization/`](odd-digit-majorization/) as [PDF](odd-digit-majorization/kaprekar-odd-digit-majorization.pdf) and [LaTeX source](odd-digit-majorization/kaprekar-odd-digit-majorization.tex). Archived at Zenodo: [doi:10.5281/zenodo.24342270](https://doi.org/10.5281/zenodo.24342270). Paper record: [doi:10.5281/zenodo.24358748](https://doi.org/10.5281/zenodo.24358748).

For every odd number of digits, the continuous large-base map decreases the vector of successive gaps in majorization order, so every real orbit approaches a periodic orbit. For nine digits, periodicity is reached in finite time, and the initial differences, arrival time and terminal differences, all divided by the base, have five pairwise distinct joint limits indexed by the base modulo 5. A covariance principle under change of base explains the modulus. The nine-digit proofs use computer-assisted affine inequalities and exact lattice counts.

**Seven-digit Kaprekar dynamics: global law, two clocks, parity** (2026), in [`seven-digit-dynamics/`](seven-digit-dynamics/) as [PDF](seven-digit-dynamics/kaprekar-seven-digit-dynamics.pdf) and [LaTeX source](seven-digit-dynamics/kaprekar-seven-digit-dynamics.tex). Archived at Zenodo: [doi:10.5281/zenodo.24346853](https://doi.org/10.5281/zenodo.24346853). Paper record: [doi:10.5281/zenodo.24358752](https://doi.org/10.5281/zenodo.24358752).

With seven digit slots and a growing base B, the normalized starting differences and the arrival time divided by B have a joint limit along the bases of each parity. At almost every starting point the limiting time takes one of two values; the faster one has probability 3/8 in odd bases and 1/2 or 3/4 in even bases. Every arrival time is at most 54B + 614132 for B ≥ 400, and at first arrival the normalized differences lie within 34/B of (3/4, 1/2, 1/4). The proofs use computer-assisted affine inequalities and exact rational calculations.

**Linear-time access for eleven-digit Kaprekar dynamics** (2026), in [`eleven-digit-access/`](eleven-digit-access/) as [PDF](eleven-digit-access/kaprekar-eleven-digit-access.pdf) and [LaTeX source](eleven-digit-access/kaprekar-eleven-digit-access.tex). Archived at Zenodo: [doi:10.5281/zenodo.24358476](https://doi.org/10.5281/zenodo.24358476). Paper record: [doi:10.5281/zenodo.24358756](https://doi.org/10.5281/zenodo.24358756).

For a uniformly random eleven-digit string in a growing base B, Kaprekar subtraction reaches a centered regime, where all six gaps between the symmetric digit differences are at least ⌈(B+1)/7⌉ + 64, within 4509B steps with probability tending to one, through all bases. Consequently the time to the first cycle is of order B, and the first cyclic word's normalized differences converge in probability to (5,4,3,2,1)/6. Exact cycles that never reach the centered regime exist for bases congruent to 0, 1 or 2 modulo 5, with vanishing probability. The proof combines almost-everywhere entry for the continuous leading map with a stopped energy estimate; its computer-assisted steps include an independent exhaustion with exact rational infeasibility certificates.

**Arithmetic selection of arrival times in eleven-digit Kaprekar dynamics** (2026), in [`eleven-digit-selection/`](eleven-digit-selection/) as [PDF](eleven-digit-selection/kaprekar-eleven-digit-selection.pdf) and [LaTeX source](eleven-digit-selection/kaprekar-eleven-digit-selection.tex). Archived at Zenodo: [doi:10.5281/zenodo.24371697](https://doi.org/10.5281/zenodo.24371697). Paper record: [doi:10.5281/zenodo.24371745](https://doi.org/10.5281/zenodo.24371745).

For a uniformly random eleven-digit word in a growing base B, the joint law of the normalized input and the arrival time divided by B depends on the arithmetic of B. A fixed test separates the bases congruent to 0 and 1 modulo 6, so the joint law has no limit through all bases. On a small rectangle of initial gaps the conditional law converges along each class of B modulo 6 to six distinct explicit mixtures of affine clocks, so the joint law converges along no progression whose modulus is not divisible by 6. On a second box the conditional law converges along the classes B ≡ 5 and B ≡ 11 (mod 12) to two distinct explicit laws, whose largest difference on [0,1]-valued tests is 17/7128; hence a system of limits along the classes modulo N would require 12 to divide N. Whether the full joint law converges along the classes modulo 12 is open. The proofs use computer-assisted certificates of affine routes and their lattice weights.

**Eleven-digit Kaprekar arrival laws modulo twelve** (2026), in [`eleven-digit-law/`](eleven-digit-law/) as [PDF](eleven-digit-law/kaprekar-eleven-digit-law.pdf) and [LaTeX source](eleven-digit-law/kaprekar-eleven-digit-law.tex). Archived at Zenodo: [doi:10.5281/zenodo.24396353](https://doi.org/10.5281/zenodo.24396353). Paper record: [doi:10.5281/zenodo.24396437](https://doi.org/10.5281/zenodo.24396437).

For a uniformly random eleven-digit word in a growing base B, couple all bases through eleven uniform real digits. On an entrance event of probability at least 92.9921494990%, the normalized initial differences, the arrival time divided by B and the normalized differences of the first cyclic word have a joint limit along every class of B modulo 12; these laws converge along every class modulo M if and only if 12 divides M. A clock kernel describes the limit: each repeated affine return is charged its limiting duration, and each admissible exit residue has probability equal to the reciprocal of its preimage index in the inherited input lattice. Some rank-one arrivals repeat only with period 66 in B, but the extra factor 11 disappears from the joint law. Convergence on the rest of the inputs is open. The proof uses exact polyhedral and lattice computations.

**Off-centre Kaprekar cycles with positive basin density** (2026), in [`off-centre-cycles/`](off-centre-cycles/) as [PDF](off-centre-cycles/kaprekar-off-centre-cycles.pdf) and [LaTeX source](off-centre-cycles/kaprekar-off-centre-cycles.tex). Archived at Zenodo: [doi:10.5281/zenodo.24399274](https://doi.org/10.5281/zenodo.24399274). Paper record: [doi:10.5281/zenodo.24399342](https://doi.org/10.5281/zenodo.24399342).

For uniformly random words of fixed odd length n in a growing base B, at all 47 strict lengths of folded order 3 to 13 (the least r with 2^r ≡ ±1 mod n), a share of inputs bounded below for all large bases ends on cycles whose peak normalized distance from the centred profile exceeds a fixed bound; at 15 and 17 digits these cycles have period eight. At 45 and 91 digits, cycles of least period 24 and of least period 132 each attract a positive share. At thirteen digits, in bases divisible by 7, a positive share reaches one period-177 cycle near the centre. The proved shares are tiny lower bounds and say nothing about the typical input. The proofs use exact certificates for guarded affine routes, lifted to actual subtraction with exact input weights.

**Macroscopic cycles in odd-length Kaprekar dynamics: a chamber criterion and a residue-transfer rule** (2026), in [`macroscopic-cycles/`](macroscopic-cycles/) as [PDF](macroscopic-cycles/kaprekar-macroscopic-cycles.pdf) and [LaTeX source](macroscopic-cycles/kaprekar-macroscopic-cycles.tex). Archived at Zenodo: [doi:10.5281/zenodo.24408684](https://doi.org/10.5281/zenodo.24408684). Paper record: [doi:10.5281/zenodo.24408723](https://doi.org/10.5281/zenodo.24408723).

For every odd word length at least seven, an exact criterion decides when macroscopic Kaprekar cycles can exist inside the chamber of positive pairing: a nonzero gap profile invariant under folded doubling and satisfying weak inequalities. When folded doubling has one orbit on the nonzero slots, a macroscopic cycle must lose pairing or breach a secondary comparison by an amount proportional to the base within two folded turns. A residue-transfer rule, valid at every odd length, gives exact limiting weights for exit residues of certified routes by an inherited lattice index. The length 1069 lies outside the natural three-way classification of lengths. No population share beyond earlier work is claimed.

**The arithmetic spectrum of Kaprekar lengths** (2026), in [`arithmetic-spectrum/`](arithmetic-spectrum/) as [PDF](arithmetic-spectrum/kaprekar-arithmetic-spectrum.pdf) and [LaTeX source](arithmetic-spectrum/kaprekar-arithmetic-spectrum.tex). Archived at Zenodo: [doi:10.5281/zenodo.24419764](https://doi.org/10.5281/zenodo.24419764). Paper record: [doi:10.5281/zenodo.24419878](https://doi.org/10.5281/zenodo.24419878).

Which odd word lengths can host macroscopic Kaprekar cycles inside the chamber of positive pairing is decided by the arithmetic of the length. At a prime length the question reduces to whether an alternating walk of incomplete character sums on the cosets of ⟨−1, 2⟩ surrounds the origin; at index two it is a Legendre walk ending at 1 − h(−4p). Unconditionally, at each fixed index, prime lengths admitting such cycles have density zero; under the Riemann hypothesis for the relevant Kummer fields, almost every prime length admits none. Every odd composite length admits them, with explicit fixed points; among composites with two or three folded orbits only 9, 15 and 21 are strict. Exact certificates classify every odd length below 10⁴. Basin mass and cycles leaving the chamber are not decided.

**Capture times, radius bounds and periodic roots in odd-length Kaprekar dynamics** (2026), in [`capture-times/`](capture-times/) as [PDF](capture-times/kaprekar-capture-times.pdf) and [LaTeX source](capture-times/kaprekar-capture-times.tex). Archived at Zenodo: [doi:10.5281/zenodo.24421912](https://doi.org/10.5281/zenodo.24421912). Paper record: [doi:10.5281/zenodo.24422109](https://doi.org/10.5281/zenodo.24422109).

Inside the chamber of positive pairing, the gaps on the folded orbit of slot 1 evolve by an integer recurrence of order r (the folded order), and this recurrence is the clock: the probability of staying in the chamber and being captured by gap time T is O(((T+1)/B)^(k−d)), with k − d = r at strict lengths and a matching lower bound at the 47 strict lengths of orders 3–13. Every wholly positively paired cycle has peak radius at most 2/k + 6/B. The loss-free periodic points of the recurrence with period dividing 2r are classified by a fixed finite graph, and an exact full-rank construction extends positive off-centre basin density to all twelve lengths of folded order 15.

## Reproducing the results

```sh
cd cycles-even-bases/ancillary
bash reproduce.sh            # numerical checks (Python 3 and a C++17 compiler)
cd lean && lake build        # Lean 4.28.0 proof of the classification
lake env lean TheoremC.lean  # prints the axioms used
```

```sh
cd arrival-times-even-bases/ancillary
bash reproduce.sh            # numerical checks (Python 3 and a C++17 compiler)
```

```sh
cd five-digit-arrival-law/ancillary
bash reproduce.sh            # numerical checks, SMT audits and the Lean build
```

```sh
cd odd-digit-majorization/ancillary
bash reproduce.sh            # solver audits and exact calculations (stages in ancillary/README.md)
```

```sh
cd seven-digit-dynamics/ancillary
bash reproduce.sh            # solver audits, exact calculations and graph replays (about five minutes)
```

```sh
cd eleven-digit-access/ancillary
bash reproduce.sh            # catalogue, work, local and centered-clock checks (about 20–35 minutes)
```

```sh
cd eleven-digit-selection/ancillary
bash reproduce.sh            # exclusion tree, route verifiers, box matching, residual routes and guarded coverage (about 13 minutes)
```

```sh
cd eleven-digit-law/ancillary
bash reproduce.sh            # exhaustive geometric and arithmetic replays and frozen panels (about two hours)
```

```sh
cd off-centre-cycles/ancillary
bash reproduce.sh            # primary exact checks and the independent reconstruction (about two minutes)
```

```sh
cd macroscopic-cycles/ancillary
bash reproduce.sh            # exact certificates, independent checks and evidence totals (about four minutes)
```

```sh
cd arithmetic-spectrum/ancillary
bash reproduce.sh            # census certificates, independent checks and both index-two scans below 10^9
```

```sh
cd capture-times/ancillary
bash reproduce.sh            # exact replays of the capture, entrance, root-map and order-15 certificates (about a minute)
```

Each `ancillary/README.md` gives the expected outputs and versions.

## Priority record

[`COMMITMENTS.txt`](COMMITMENTS.txt) lists SHA-256 digests of research materials that are not yet released.

## License

Text, figures and data: CC BY 4.0. Code: MIT. See [`LICENSE.md`](LICENSE.md).

## Contact

rustem@galyavin.com
