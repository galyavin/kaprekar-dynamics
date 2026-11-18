# Independent verification of the nine-digit results

The [review](independent-verification.md) describes which implications were re-encoded, which route weights were tested statistically, and where the original descent certificate supplies the all-base proof. Its line references refer to the 17 November manuscript. The programs and historical logs accompany it.

Run `bash reproduce.sh independent` from the parent directory. The driver builds `kap9.c` with `cc -O2`, runs the independent arithmetic and Z3 checks, enumerates the 65 descent bases, simulates the route and uniform populations, recomputes predictions, and compares the pooled statistics. It uses a temporary working directory with an empty set of generated results. The scripts derive predictions from the distributed kernel table. The C simulator obtains observations by digit subtraction. The numerical libraries and solver versions are pinned in the parent directory.

## Reading the logs

- The 39 queries per solver seed comprise 38 theorem negations (UNSAT) and the L2d boundary example (SAT). The earliest `logs_z_entry.txt` tried a stronger boundary claim, found a counterexample, and retains its original failure marker. The distributed encoding labels L2d as an informative SAT case. Both seeded reruns record the corrected expectation in `logs_z_reseed.txt`.
- The review's §2 sentence about non-vacuity has exceptions. In `logs_z_vacuity.txt`, the first return at step four is impossible at both tested ranks. One attempted mutation (reducing the minimum-gap coefficient from 100 to 2) still gives UNSAT. Three further mutations yield SAT witnesses in `logs_z_vacuity_iv.txt`. The historical logs preserve these outcomes.
- The paper's Z coordinate is already normalized by the base. In the numerical discussion, expressions written as “Z_B/B” denote the raw terminal differences divided by the base once. Simulation scripts store raw differences and form this quotient once.
- The two simulation families can agree within sampling uncertainty. Such agreement supports their compatibility with the law. A finite-sample homogeneity test cannot establish exact equality of limiting distributions. The first-column difference between classes 2 and 4 follows from the exact kernel. The reported third-column separation also appears in the simulations.
- The historical `logs_allroutes_150k_ungrouped_superseded.txt` records a discarded first pass. The current analysis groups atoms within the stated finite-base resolution. The main pooled fit has 663 degrees of freedom. The four groups selected for replication receive independent samples, as described in §5.3.
- `logs_*_lock.txt` record elapsed waits in the historical runs. They contain no solver results.

## Files omitted from the distribution

The two all-route samples (`allroutes_raw_300k.npz`, 9,343,578 bytes, and `allroutes_raw_450k.npz`, 9,485,419 bytes) and two replication samples (`targeted_sameB.npz`, 962,066 bytes, and `targeted_bigB.npz`, 1,015,197 bytes) are regenerated during the independent stage. Their grouped counts remain in `allroutes_analysis_families.json` and `targeted_analysis.json`.

The original ten `uniform_B*.npy` arrays were discarded after analysis. Their regenerated copies occupy 115,201,280 bytes in total. The fixed seeds in `sim_uniform.py` regenerate them at bases 60010–60014 and 60035–60039. The platform-specific executable `kap9` is built from the included C source. The directory retains every historical text log, including the superseded diagnostic.
