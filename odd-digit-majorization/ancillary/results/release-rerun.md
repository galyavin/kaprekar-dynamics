# Release-copy verification, 18 November 2026

All computations ran in new directories outside the repository, using a private copy of the pinned Python environment. The seven original stages used individual invocations of `bash reproduce.sh STAGE`. Running `descent` separately includes its zero-preparation prerequisites a second time. Each independent substage ran from the supplied scripts after compiling `kap9.c` with `cc -O2`, without any pre-existing generated data. The independent replay comparisons were applied to the fresh output. No stored solver result supplied a proof verdict. Elapsed times below exclude waits for access to the machine.

Environment: macOS 27.2.1 on arm64, Python 3.9.6, SymPy 1.14.0, mpmath 1.3.0, Z3 5.4.1 (z3-solver 5.4.1.0), cvc5 1.4.2, Yices 2.7.0, NumPy 2.0.2, SciPy 1.13.1 and Apple Clang 21.0.0 (clang-2100.4.1.1).

| Stage | Wall time (seconds) | Fresh observation |
| --- | ---: | --- |
| [entry](entry.log) | 295.3 | Four universal entry negations UNSAT in both engines. Two false strengthenings SAT and three corruptions rejected. All 320 compressed captures matched. |
| [routes](routes.log) | 89.6 | All 13 scripts passed. All 33 routes, 16 inverse branches, five residues and 2,640 compressed captures checked. |
| [weights](weights.log) | 129.7 | All 165 kernels rebuilt, four corruptions rejected and all 15,840 compressed trials matched. |
| [covariance](covariance.log) | 4.1 | Six sorting orders and 78 guard directions passed. Four matched captures used 1,029,184 subtraction steps. |
| [zeros](zeros.log) | 74.8 | First access and all 21 return/one-step queries passed in both engines. The complete finite graphs at bases 128 and 129 were regenerated. |
| [descent](descent.log) | 361.5 | All 29,534 coverage nodes on 240 roots UNSAT in Z3 and Yices. Two SAT controls and four rejected corruptions. The arithmetic replay checked 27,965 branching witnesses and 120 full captures. |
| [fresh](fresh.log) | 41.7 | 20,000 of 20,000 direct trajectories matched a unique atom. The five count rows reproduce the table in RERUN.md. |
| independent | 2571.7 | The independent model, solvers, finite enumerations, fixed-seed simulations, prediction integration and grouped comparisons reproduced the archived evidence. |

Sum of the observed stage times: 3,568 s (0.99 h), including the repeated zero prerequisites. The default invocation shares those prerequisites.

| Independent substage | Wall time (seconds) |
| --- | ---: |
| Independent C build | 0.2 |
| Independent model, exact arithmetic and SMT (four concurrent checkers) | 202.7 |
| Independent exhaustive and stratified descent | 81.5 |
| Separating route, uniform strings and finite maximum clocks | 927.5 |
| 48,000 prediction inputs (eight workers) | 84.0 |
| All 165 kernels at two families of bases | 1028.2 |
| Targeted replications | 246.7 |
| Final statistical comparisons | 1.0 |

The 65 exhaustive descent totals sum to 811,157,145 qualifying sources. The stratified and slow-family tests add 3,132,549 sampled sources. All have zero descent violations. The two route families and the targeted replications match 1,086,000 observations to their prescribed groups. The separating route adds 1,500,000 orbits, and the uniform experiment adds 1,600,000 strings. The two all-route families give pooled χ² = 673.4 on 663 degrees of freedom. The targeted checks give p = .42 and .80.

All twelve independently generated JSON summaries and all four regenerated compressed route/replication sample archives matched the historical files byte-for-byte. The archived [review](../independent-review/independent-verification.md) describes the reach of the computations. The [reading notes](../independent-review/README.md) identify the SAT controls, vacuous cases and superseded log.
