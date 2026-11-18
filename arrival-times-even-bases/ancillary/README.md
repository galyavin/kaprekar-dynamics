# Arrival times for the four-digit Kaprekar map in even bases

Contact: rustem@galyavin.com

This bundle is licensed as described in [LICENSE.md](LICENSE.md).

Run the following command in this directory with Python 3.9 or later and a C++17 compiler. Python uses only its standard library. The environment variables `PYTHON` and `CXX` can select the executables.

```sh
bash reproduce.sh
```

The script compiles both C++ programs in a temporary directory, recomputes the exact results, and compares them with the two JSON files in [results](results/). It also checks the SHA-256 manifest of the source, documentation, and results. A missing result or manifest is generated. A mismatch stops the run. On an Apple silicon Mac with Python 3.9.6 and Apple clang 21.0.0, allow about 10 seconds including compilation. The final output reports elapsed time and software versions.

## Section 8 computations

[The results](results/section8.json) include the following comparisons.

| Computation | Expected result |
|---|---|
| Formula (5) and weighted predecessor layers | All coefficients agree in 300 target/exponent cases, odd parts 5 through 31, exponents 1 through 3. |
| Literal digit subtraction and weighted predecessor layers | All coefficients agree for 146,904 nonconstant inputs counted across seven target cases, with bases 5, 10, 7, 14, 13, 13, 14. |
| Full difference-pair graph at base 5120 | All 13,109,759 pairs enumerated. Every coefficient and both moment sums agree with (5). |
| Largest finite time for the fixed point at odd part 5 | All 40 maxima agree with Devlin and Zeng, Theorem 1. |
| Table 2 and conditional moment coefficients | Four exact rows. Their sums are 1033/1125, 1, 94/125, 1 in phase order. |
| Exact total variation distance to the uniform-law model | All 20 distances for odd parts 5, 9, 31, 53 and exponents 0 through 4 satisfy the stated bound. |

[run_checks.py](run_checks.py) runs the programs under [source](source/). The polynomial implementation in [symbolic_arrivals.py](source/symbolic_arrivals.py) generates the powers of four by successive multiplication. The large-base C++ enumeration works directly from the digits of subtraction outputs.

## Independent computations

[independent.py](independent.py) evaluates (5) by separating residue and valuation classes. It sums their quadratic weights by finite differences. An independent numerical calculation sorts and subtracts each four-digit input and follows its orbit. The six comparisons include a proper union of two cycles at base 26 and the empty target at base 14. In base ten the coefficients at times 0 through 7 are 1, 383, 576, 2400, 1272, 1518, 1656, 2184.

[graph_check.cpp](graph_check.cpp) constructs the forward graph by subtracting the positional words with digits (u,v,0,0) in descending and ascending order. It identifies the cycles, entrance distances and maximal times for every basin. For odd parts 5 through 41 and exponents 1 through 6, the arrival times of the critical states E, U, Z(2^t) (including S at t=n) and A(B/2+2^t) agree with Lemma 7 for each cycle and their union. The largest finite input time is 31. For odd parts 27, 29, 31 and 53 and exponents 0 through 4, every scaled total variation distance is below 3.94, against the bound 26.

[The independent results](results/independent.json) give both coefficient vectors, the boundary comparisons, the four limit rows, exact basin proportions through exponent 40, and all 20 rational distances. In the boundary records, cycles are indexed starting at zero in increasing order of their least pair. The entries in `H` follow E, U, the states Z(2^t) in increasing t including S at t=n, and then A(B/2+2^t) in increasing t. A null value denotes infinity. Each target’s maximum is the largest finite numerical arrival time. For basin proportions, `raw` uses all digit strings as denominator and `eligible` uses the nonconstant strings.
