# Reproducing the five-digit paper

Contact: rustem@galyavin.com

This bundle is licensed as described in [LICENSE.md](LICENSE.md).

From inside this folder, run

```sh
bash reproduce.sh
```

The default `all` run executes the three solver audits, exact integration,
complete direct-subtraction comparisons, maximum witnesses, the five-digit
overshoot checks, and the manuscript-number checks. It compiles the separate
C++17 reference implementation and regenerates its comparison histograms.
It then runs the independent 43-check review, a second exact CDF calculation and complete weighted graphs at bases 30–200, and the Lean build and fresh kernel replay. The default command passed from a standalone copy in 132.15 seconds on an Apple M5 (18 November 2026 UTC).

The runner copies the bundled snapshot in `source/` to `results/run/`, where it regenerates the numerical outputs. Logs, software versions, source SHA-256 hashes and the status of each check are in `results/`.

`results/artifact_sha256.json` records the SHA-256 of every distributed source and result, using paths relative to this folder. Check these hashes before running the driver. A new run records its own software versions and elapsed times. From a fresh distribution, verify the listed bytes with

```sh
python3 - <<'CHECK'
import hashlib, json
from pathlib import Path
manifest = json.loads(Path('results/artifact_sha256.json').read_text())
for name, digest in manifest.items():
    if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
        raise SystemExit(f'Checksum mismatch: {name}')
print(f'{len(manifest)} files verified')
CHECK
```

## Requirements and exact versions

The recorded run used Python 3.9.6, Z3 5.4.1 (`z3-solver==5.4.1.0`) and
Apple Clang 21.0.0 with C++17. Python rational arithmetic provides all
integration enclosures. The default run also uses Lean 4.28.0 and its bundled Lake (5.0.0-src+7e01a1b), on the executable search path. The Lean project imports the bundled standard library. The long independent-review replay additionally uses NumPy 2.0.2 and a C compiler with 128-bit integers. `results/versions.json` records the actual executable versions.
No solver environment is silently installed or changed.

The shell script uses `python3` by default. To prepare a local environment with the recorded solver version, run the following commands inside this folder before running the checks.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install z3-solver==5.4.1.0
```

Set `PYTHON` to select another interpreter with Z3 installed. Set `CXX` to select a C++17
compiler. Run with assertions enabled. The numerical driver rejects `-O`, `-OO` and `PYTHONOPTIMIZE`. The audits require UNSAT for every counterexample query on an asserted theorem and SAT for the missing-branch control. An unexpected solver result, a failed assertion or a nonzero subprocess exit fails the run.

## Statement-to-command map

Each stage uses `bash reproduce.sh STAGE`.
The default `all` runs the numerical stages below, then the independent Z3 check, `headline` and `lean`. Theorem constants are proved
in the paper; numerical reproductions do not replace those arguments.

| Stage | Encodings / calculations | Paper statements and outputs |
|---|---|---|
| `audit` | `five_clock_audit.py` | Coverage and disjointness of the 26 branch guards for every integer base ≥ 30, with a missing-branch control; all-base core formula and invariance; first-step bounds; six-step strip entry; two-step descent; all 144 terminal starts, maximum tail 10 and errors -3 and 16/3. Writes `five_clock_audit_results.json`. |
| `audit` | `five_joint_audit.py` | The 26 rare-predecessor and residue checks; four residue bijections; terminal labels and periods. Writes `five_joint_audit_results.json`. |
| `audit` | `five_quantitative_audit.py` | All seven first-core-entry margin implications, for every integer base at least 30. Writes `five_quantitative_audit_results.json`. |
| `limit` | `five_limit_checks.py`, `five_limit_law.py`, `five_holdout.cpp` | Exact F(1/25), normalized polygons, missing-mass polynomial via retained mass checks, moments at cutoff 300; 121,024 direct pairs, 10,080 core states, 4,560 rational leading points. Recreates complete references at bases 313–318 and 1553–1558. Writes `five_limit_results.json` and individual histograms. |
| `joint` | `five_joint_checks.py`, `five_joint.py`, `five_holdout.cpp` | 223,344 direct pairs, 18,606 core labels; exact additional basin weights and actual-time histogram; joint CDF comparisons. Recreates complete references at bases 431–436 and 1973–1978. Writes `five_joint_results.json` and individual histograms. |
| `maximum` | `five_quantitative_checks.py`, `five_quantitative.py` | 241 constructed witnesses: bases 30–250, 1000–1010, 1973–1978, 100003–100005. Checks their actual clocks against the sharp bounds, and existing joint-reference maxima. Evaluates conservative CDF bounds at 2^40, 2^80, 2^120. Writes `five_quantitative_results.json`. Run `joint` first for the reference-maxima portion. |
| `overshoot` | `p2_overshoot_checks.py`, `arithmetic_scaling_checks.py`, `five_boundary_layer.py` | Only the five-digit source routine runs: 9,208 clocks and 3,971 direct high steps; verifies actual center error at most 8. Weighted exit blocks at bases 10000, 40000, 160000. Writes `p2_overshoot_results.json`. |
| `numbers` | `p2_numbers.py` | Exact fraction 713903/39062500; tail 38348597/39062500 = 0.9817240832; strict decimal bounds 0.1827011465 and 0.1827011669 for the mean; independent ordered-input weight enumeration (61,740 strings, bases 2–8); arithmetic check of the conservative CDF constant. Writes `p2_numbers_results.json`. |

`source/five_excursions.py` generates the common guarded sort branches;
`source/five_digit.py` supplies elementary digit operations and weights;
`source/five_clock.py` supplies the stated core and limiting clock formulas.
The four-digit modules in the snapshot are import dependencies of the
unchanged mixed-digit check file; its four-digit campaign is not run here.
The bundle modifies the research copy of `source/five_clock_audit.py` to check the partition before using it in the clock proof. It asks whether a valid state satisfies zero guards or at least two, using exactly the guards from `five_excursions.py`. A third query omits branch `Z0` and requires an uncovered state. The remaining files copied from the research directory are byte-identical. The two additional scripts are `p2_overshoot_checks.py` and `p2_numbers.py`.

The paper's CDF error constant and boundary-layer threshold are conservative
analytic estimates. They are not solver outputs. Exact integer checks and
finite trajectory comparisons are listed separately in the appendix.

## Independent headline checks and Lean

`bash reproduce.sh headline` uses exact triangle cubature for (2) and positional subtraction for the pair graphs. It must print 713903/39062500 and its complementary tail 38348597/39062500. All 171 bases must pass (9) and (37), including the full histogram and three-pair basin at each of the 86 even bases. The per-base records are in `results/headline_results.json`.

`bash reproduce.sh lean` runs the following commands. The expected build reports nine completed jobs. Across the 26 axiom listings the union is `propext`, `Classical.choice`, and `Quot.sound`. The fresh replay exits with status zero without additional output. A previously recorded clean build took about 20 seconds and the replay about 66 seconds on this machine.

```sh
cd lean
lake build
lake env leanchecker --fresh FiveDigit
```

## Paper production record

The earlier Markdown draft produced 17 pages with Pandoc and Tectonic. Its numerical verification completed in 28.32 seconds before the Lean and additional review stages were added.

Two clean pdfLaTeX passes with TeX Live 2026 produced 16 pages. Every page was rendered and inspected. The second pass reported one underfull line in Appendix A.2 (badness 1033), zero overfull boxes and zero unresolved references.

## Reference verification, 17 November 2026

The entries below were checked against the primary source, including the
location actually cited. URLs are also in the paper. Source inspection is
bibliographic checking, not a claim of exhaustive historical-priority review.

| Ref. | Verified bibliographic details and cited location |
|---|---|
| 1 | Gordon D. Prichett, *Terminating cycles for iterated difference values of five digit integers*, J. reine angew. Math. 303/304 (1978), 379–388. Göttingen digitization: opening page, Section 1, Section 4 and Tables 4.1–4.3 (pp. 387–388). The table's T1–T5 labels agree with the paper's definitions. |
| 2 | Pat Devlin and Tony Zeng, *Maximum distances in the four-digit Kaprekar process*, Integers 21 (2021), A97. Publisher PDF p. 1 supplies metadata; Theorems 1–2 (pp. 3–4) give the cited maximum/fixed-point-basin scope. |
| 3 | Evan Chen, Ken Ono, Richard E. Schwartz and Dinesh Thakur, *Four-Digit Kaprekar Dynamics in Odd Bases*, Journal of Integer Sequences 29 (2026), Article 26.4.7. Publisher PDF title page lists all four authors; Theorems 1 and 24 give the structural and cycle-count statements. arXiv:2606.20439v3, revised 11 August 2026, agrees. The publisher's HTML abstract page header lists only Chen and Ono; the PDF and volume index establish the full authorship. |
| 4 | Christoph D. Dahl, *Coarse-Grained Drift Fields and Attractor-Basin Entropy in Kaprekar's Routine*, Entropy 28(1) (2026), 92, DOI 10.3390/e28010092; published 12 January 2026. Publisher article/PMC metadata, Sections 2–3 and Appendix A were checked for the finite-base statistical scope. |
| 5 | Christoph D. Dahl, *Information funnels and multiscale gap-space dynamics in Kaprekar's routine*, arXiv:2512.05124 (2025). The arXiv record confirms the distinct preprint title. It is cited separately so that the journal title is not assigned to the preprint. |

## Independent verification

The reviewer independently reconstructed the proof reductions using a sorting-network transition. The review includes 43 checks and fresh-base tests through 2,560,009. Run the independent Z3 checks from this folder with

```sh
"${PYTHON:-python3}" independent-review/z3_independent.py
```

On 17 November 2026, `bash reproduce.sh audit` passed all three audits in 0.73 seconds in a standalone copy outside the repository. Z3 returned UNSAT for both an uncovered state and a state satisfying multiple guards across the full integer domain `B ≥ 30, 1 ≤ u < B, 0 ≤ v ≤ u`. Removing `Z0` returned SAT with the uncovered state `(B,u,v) = (31,16,0)`. The independent command passed all 43 checks in 2.61 seconds in the same copy. Both times are elapsed wall time with Python 3.9.6 and Z3 5.4.1.

The full historical calculations, including the two complete base sweeps, can be rerun with `bash reproduce.sh review-full`. Install NumPy in the selected Python environment for this command (`python -m pip install numpy==2.0.2` for the recorded version). Results go to `results/review-replay/`. The full review replay took 962.74 seconds (about 16 minutes) on the same machine. Both sweeps had zero failing bases. The fresh overshoot run checked 620,000 samples across the two bands, with maximum absolute error 5 in (44) and zero high-branch or exit-interval failures. The two overshoot bands use seed 42 because the historical logs omit their seeds. Both orbit implementations must agree at each directly simulated base. The comparison to archived numerical logs uses the recorded toolchain and its printed rounding.

The [independent review](independent-review/README.md) includes the original report, logs and individual commands. Its report concerns the 17 November draft. Its recommended update on completed verification appears in Section 11 of the release manuscript. The 43 successful checks in its summary include 39 UNSAT queries and four numerical or enumerative checks.

The historical pair-graph logs print large individual cycle counts through binary64; their last integer digits can be rounded. Internal weights and total-weight checks use 128-bit integers.
