# Independent verification of the five-digit arrival law

From the bundle root, with Python and Z3 available as described in the bundle README, run

```sh
"${PYTHON:-python3}" independent-review/z3_independent.py
```

The program derives the transition from a sorting network and reports the number of checks passed out of 43. A failed check gives a nonzero exit status. The reviewer's findings and proof reductions are in [independent-verification.md](independent-verification.md).

| Files | What they check |
|---|---|
| `z3_independent.py`, `z3_independent.log`, `z3_probes.log` | Integer-base proof checks and finite terminal enumeration. The optional `--probe` argument also tests weakened claims and hypothesis feasibility. |
| `compare_branches.py`, `compare_branches.log` | Coverage and disjointness of the 26 guards for every integer base ≥ 2, and agreement of each affine image with the sorting network. |
| `limit_exact.py`, `limit_exact.log` | Exact rational integration for the limiting CDF, retained mass and mean. |
| `tau_mc.py`, `tau_mc.log` | Monte Carlo comparison for the limiting CDF and mean. |
| `strings_exhaustive.c`, `strings_exhaustive.log` | Full five-digit string graphs, ordered-input weights, cycles and arrival histograms at bases 4–36. |
| `pairgraph.c`, `pairgraph_fresh.log`, `pairgraph_10009.log` | Exact weighted pair graphs at selected bases through 10009, including maximum arrival time and cycle weights. |
| `sweep_pairgraph.py`, `sweep_30_1200.log`, `sweep_1201_3000.log` | Finite-base reduction checks and cycle labels at each base from 30 through 3000. |
| `mc_fresh.c`, `mc_fresh.log` | One million independent input draws per base using direct numerical subtraction. |
| `summarize_fresh.py`, `fresh_table.md` | Comparison tables from the pair-graph and direct-subtraction logs. |
| `overshoot_mc.c`, `overshoot_fast.c`, `overshoot_limit.py`, `overshoot.log`, `overshoot_large.log` | Weighted sampling in two overshoot bands, with direct and accelerated orbit evaluation and numerical limit predictions. |
| `lockrun.py` | Wrapper used to serialize the original runs on the reviewer's computer. Set `COMPUTE_LOCK` to share a specific lock file. |

## Optional computations

The commands below run from `independent-review/`. The two sampling and quadrature scripts `tau_mc.py` and `overshoot_limit.py` require NumPy. The C programs require a compiler with 128-bit integer support, such as Clang or GCC.

```sh
python3 compare_branches.py ../source
python3 limit_exact.py 300
python3 tau_mc.py 20000000
python3 summarize_fresh.py

for name in strings_exhaustive pairgraph mc_fresh overshoot_mc overshoot_fast; do
    cc -O2 -o "$name" "$name.c" -lm
done
for base in $(seq 4 36); do ./strings_exhaustive "$base"; done
python3 sweep_pairgraph.py ./pairgraph 30 1200
python3 sweep_pairgraph.py ./pairgraph 1201 3000
for base in 997 1500 2401 2999 3000 5003 10007 10009; do ./pairgraph "$base" checks; done
for base in 997 1500 2401 2999 3000; do ./mc_fresh "$base" 1000000 "$((7919 * base + 17))"; done
```

For example, these commands compare the two overshoot implementations on the same sample and compute the predicted limit. The two logs list the bases, bands and sample counts for the larger campaign. Their headers omit the random seeds, so the commands below use an explicit example seed.

```sh
./overshoot_mc 10007 0.5 1.5 0.6 0.9 20000 42
./overshoot_fast 10007 0.5 1.5 0.6 0.9 20000 42
python3 overshoot_limit.py 0.5 1.5 0.6 0.9
```
