# Fresh-copy verification, 17 November 2026

The encodings, certificates and kernel table were copied to a new temporary directory outside the project. The scripts ran from that copy against a newly created Python environment. Generated JSON, intermediate solver queries and compiler outputs stayed there. Historical progress messages were removed from some bundled checker summaries after replay; the mathematical checks are identical. Each substantive check held the project compute lock separately. The stage names below correspond to the public driver.

The system was macOS on arm64 with Python 3.9.6, SymPy 1.14.0, mpmath 1.3.0, z3-solver 5.4.1.0 (Z3 5.4.1), cvc5 1.4.2, Yices 2.7.0 and Apple Clang 21.0.0. Each selected Python checker ran as `python -B -u NAME.py` from the copied `source` directory. The public `./reproduce.sh fresh` entry point was also run from a second clean bundle copy. No stored solver verdict was used to certify a universal implication.

| Calculation | Fresh result |
|---|---|
| Core geometry and returns | All 13 route scripts passed. Reconstructed 45 core guard implications and both release inventories. Checked all 33 routes, 16 inverse branches and 5 residues in 2,640 compressed captures; 33 full difference trajectories took 4,524,703 updates and two full number trajectories took 274,777. |
| Four real entry lemmas | All four negations UNSAT in both Z3 and cvc5. Two false strengthenings SAT in both. Three deliberately corrupted results rejected. |
| Global prefixes and terminal phase | All 320 compressed captures matched a unique atom with the first cyclic number’s differences independently reconstructed by whole-integer subtraction. The tested entry lengths were 2, 5, 8 and 631. There were 716 finite-base guard rejections before increasing the base. Two full pair paths used 320,648 updates. |
| Core lattice and probabilities | All 165 kernels reconstructed and normalized. Four corrupt tables rejected. All 15,840 finite compressed trials matched a unique joint atom. The maximum coordinate discrepancy in these trials was 29,651 / 720,000,000,000. |
| Covariance | All six nine-slot sorting orders and 78 guard directions passed. The four full matched captures used 1,029,184 raw steps across the source and translated systems. |
| First zero, return cover and one-step inequalities | The odd-base first-access negations passed in Z3 and cvc5, including the shared seven-slot control. All 21 return-cover and one-step negations were UNSAT in both engines. The separate C++ replay checked full functional graphs at bases 128 and 129. |
| Uniform descent | 240 root cases, 29,534 coverage nodes, no missing child. Every node was UNSAT in both Z3 and Yices. The full audit included 59,075 solver calls, two SAT controls and four rejected corruptions. Independent arithmetic replay checked 27,965 branch witnesses, the eight-step counterfamily, and 120 full number captures. |
| Fresh finite bases | All 20,000 direct trajectories matched a unique clock and terminal atom. The public driver reproduced the same counts. Results below. |

## Five fresh bases

The C++ program reconstructs the resulting number's nine digits by ordinary column borrowing, sorts them, and finds the first repeated difference phase with Floyd's algorithm. It advances one step to read the differences of the first cyclic number. For every base it takes 4,000 uniform integer points in a box of radius 210 about the original difference point (43,39,31,28)/50 times the base (center coordinates rounded down). The sixteen box vertices pass every defining inequality of the separating route in exact arithmetic. Every observation was compared against all four joint atoms. This is a finite test on one interior box per base class. Exact differences between the five probabilities follow from the image-lattice count.

| Base | 7/2 | 9/2 | 50/9 | 37/6 | Maximum normalized clock error | Maximum normalized terminal error |
|---:|---:|---:|---:|---:|---:|---:|
| 300035 | 954 | 407 | 1866 | 773 | 0.000133688 | 0.000016665 |
| 300036 | 916 | 386 | 1950 | 748 | 0.000120726 | 0.000018332 |
| 300037 | 770 | 411 | 2155 | 664 | 0.000125170 | 0.000014999 |
| 300038 | 767 | 432 | 2156 | 645 | 0.000111097 | 0.000014999 |
| 300039 | 799 | 433 | 1992 | 776 | 0.000136279 | 0.000014999 |

The largest absolute difference of an observed proportion from its rational prediction was 341/28,800. An initial test design at base 75,035 crossed a route boundary (one sampled point had minimum defining slack −3364/75035). Those data do not test the specified conditional law. The final program checks containment of the entire sampling box before measuring any trajectory.

## First cyclic number

For the original differences obtained by rounding down base 2,000,003 times (43,39,31,28)/50, the evaluator gives difference time 251,121 and number time 251,122. The corresponding difference vectors are (1600004,1200003,800004,400002) and (1600003,1200003,800001,399999). Whole-integer subtraction confirms the latter as the successor. The new global evaluator returns both vectors explicitly.

## File identity and builds

- Descent tree SHA-256: `ceb0ca4f02f7de12538fd32fc3e10654aa3cc5520cc8fa16d3502977e769c429`.
- Kernel table SHA-256: `b8432d441ea05d3ba91cb7c2285d244738025e47060a2627163fe34ea798732d`.
- Paper compiled in independent clean build directories with pdfLaTeX (pdfTeX 1.40.29, TeX Live 2026), BibTeX 0.99e, and Tectonic 0.17.0. Both produced 14 pages. The supplied PDF is the pdfLaTeX build.

An independent review re-derived the reductions and tested them with separately written encodings. Its report records the results, the statistical checks and their limits in [independent-review/](independent-review/independent-verification.md).
