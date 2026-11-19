# Verification for Linear-time access for eleven-digit Kaprekar dynamics

Contact: Rustem Galyavin <rustem@galyavin.com>

The text and data are licensed CC BY 4.0 and the code is licensed MIT. See [LICENSE.md](LICENSE.md).

## Requirements and command

The recorded environment used Bash, Python 3.9.6, Z3 5.4.1 and Yices 2.7.0. Activate a virtual environment outside this directory and install `requirements.txt` (z3-solver 5.4.1.0, SymPy 1.14.0 and mpmath 1.3.0). Add the Yices SMT-LIB executable `yices-smt2` to PATH. The Z3 wheel supplies both the Python package and the command-line executable. Set `PYTHON` or `Z3` to an absolute executable path if needed.

From a fresh copy of this directory, run

```bash
bash reproduce.sh
```

The program verifies the supplied file inventory and SHA-256 digests, then prints a private output directory for all generated files. Use `--workers N` to set the number of centered-coverage workers, which defaults to four. Allow about 30 minutes for standalone verification with four workers.

The checks run in this order.

1. Verify the input manifest and pinned software versions.
2. Check the local radius predicates, adjacent-level sorting and returns, affine cycles and their inverse basin, and guarded-return quadrants. Check the exact sharp ray factor and independently reconstruct the four-cycle basin at five bases.
3. Reconstruct and audit the independent catalogue, compare it with the archived catalogue, and replay both archived sorting trees in both solvers.
4. Regenerate both sorting trees and the exact energy and work tests, then compare again with the independent catalogue.
5. Check the centered lower hits and the coverage of every node in the centered tree.

## Expected results

- The generator builds both leading sorting trees. The positive-outer-gap tree has 38,937 queries and 906 oriented return words. The zero-outer-gap tree has 11,937 queries and 172 oriented return words. Both solvers reproduce every archived decision.
- The independent search and exact-fraction checker reconstruct 791 sorting faces, 251 zero-loss faces, 165 stationary-outer-gap faces and 268 cyclic cells. They validate 122,278 rejected systems using 116,556 distinct reduced rational Motzkin certificates. They replay all 1,078 oriented witnesses with exact sorting.
- The independent arithmetic audit checks strict cell membership, energy extrema, and 748 work vertices in 189 cells. It reproduces the global work constant 2/7 and the full-closure vertex-test constant 4/7 for cells meeting the energy ball. The comparison with the archived and newly generated catalogues reports no differences.
- Both solvers refute the six radius counterexample predicates and the six adjacent-level sorting and return predicates. A second encoding checks the radius predicates with Yices, and a pairwise-rank encoding checks the sorting and return predicates with both solvers.
- Both solvers refute the fifty unbounded guarded-return quadrant counterexamples and the included reduction counterexamples. The deliberately false two-observation bound has a replayed witness.
- The strict ten-step ray query is UNSAT in both solvers at the bound 295/77. The point (−53, 12, −15, 77)/77 attains the bound with image (295, 256, 57, −295)/77.
- Independent reverse search at m = 20, 21, 50, 1000 and 20001 gives 235 states with the affine layer sizes, including 45 zero-innermost states. At m = 20 it matches the supplied exhaustive list state for state and depth for depth.
- The inherited centered lower-hit checker accepts all 97,465 nodes. It verifies branch coverage with fresh solver calls and splits any unsettled Yices query into all six physical residue classes.
- Exact affine subtraction confirms periods 4, 8, 8, 10 and 8 in base classes 0, 1, 1, 1 and 2 modulo 5, respectively. The complete reverse search in class 0 gives the layer sizes 4, 27, 31, 38, 57, 36, 25, 17 and 0.

## Files

`original/` holds the incremental tree transcripts, status lists, feasibility witnesses and exact arithmetic results. `generator/` reconstructs the tree and its work and energy tests. `independent/` specifies a different exhaustive search with an exact rational checker. Its frozen records include all feasibility witnesses, energy and work vertices, and the rational certificates. It also contains the strict ray query, the exact ray linear program, the full-cell affine return check and the [inverse-basin reconstruction](independent/basin_reconstruction.md), with the state list at m = 20. `local/` contains radius predicates, adjacent-level sorting predicates, and the affine cycle and inverse-basin checks. `centered/` contains the guarded-return and centered lower-hit verifiers and the finite centered tree. `results/` holds archived decision and comparison records and the SHA-256 manifest.

## Independent verification

[independent-review/report.md](independent-review/report.md) describes the separate checks of an earlier draft and their scope. The inputs, SMT-LIB encodings and summary logs are in `independent-review/`. The SMT checks use cvc5 1.4.2 through its Python package and Yices 2.7.0 through `yices-smt2`. The word-subtraction experiments use a C compiler with optimization and an unsigned 128-bit integer type.

From the bundle root, run these two standalone checks with Python's standard library.

```bash
python3 -B independent/ray_expansion.py
python3 -B independent/affine_returns.py
```

- The ray program uses exact rational linear programming on the eight facets of the four-dimensional unit box, with the guard `c >= d`. There are 79 feasible ten-step branches and the minimum output norm is 295/77. It takes several minutes.
- The affine reduction reconstructs the sorting and exact return identities throughout all 268 cells in `original/energy.json.gz`. Grouping each cell by the shortest return identity valid throughout it gives 16 at length 1, 14 at 2, 80 at 3, 153 at 5, 2 at 6, 2 at 7 and 1 at 12. It takes about one second.

To regenerate the [forward log at m = 19 and m = 18](independent-review/logs/basin_m19_m18.log), build and run the exhaustive forward program.

```bash
run_dir=$(mktemp -d)
cc -O3 independent-review/basin.c -o "$run_dir/basin"
{ "$run_dir/basin" 19 && "$run_dir/basin" 18; } > "$run_dir/basin_m19_m18.log"
cmp "$run_dir/basin_m19_m18.log" independent-review/logs/basin_m19_m18.log
```

The forward runs take a few seconds each. At each base, the four-cycle closes and the basin has 235 states with layer sizes 4, 27, 31, 38, 57, 36, 25 and 17. For minimum gaps 0, 1, 2 and 3, the multiplicities are 141, 48, 30 and 16. The sum of `minG + 1` over the basin is 391.

To repeat the SMT checks and the simulation campaign, use a fresh copy of the bundle and the following commands. The archived campaign ran for 107 seconds with nine concurrent jobs and produced 461,060 sample records at 55 bases.

```bash
python3 -m pip install cvc5==1.4.2
run_dir=$(mktemp -d)
python3 -B independent-review/smt_checks.py "$run_dir/smt"
cc -O3 independent-review/sim.c -o "$run_dir/sim"
python3 -B independent-review/campaign.py "$run_dir/sim" "$run_dir/samples"
python3 -B independent-review/summarize.py "$run_dir/samples"
```

The archive `independent-review/logs/sim_logs.tar.gz` (5,882,492 bytes) is omitted. It contained the raw per-sample logs that the campaign regenerates. The per-base summaries and job timings are included.
