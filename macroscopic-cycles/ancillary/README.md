# P10 reproduction

Run the following from a fresh copy of the paper folder.

```sh
bash ancillary/reproduce.sh
```

From inside the ancillary folder, run `bash reproduce.sh`. The command needs Bash, `shasum` and Python 3.9 or newer. It uses the Python standard library and requires no external solver, additional data, network access, or write access to the parent repository. Assertions are enabled. A successful run ends with `PASS: P10 exact certificates, independent checks and evidence totals.`

## Checks

| Stage | Result reproduced |
|---|---|
| `rc11_reconstruction/cones.py` | Complete zero-loss homogeneous circuit enumeration at thirteen digits, exact feasibility and work tests, and the exceptional period-two record |
| `rc11_reconstruction/physical.py` | Paired update, entrance and cycle controls, the 1069 orbit and row data, exact prime-square fixed points, and the nine-digit insertion certificate with physical checks on all five cosets |
| `grc/check_instances.py` | Nine-digit moving weights, eleven-digit exit indices, and isolated-wall sections using the original finite route implementations |
| `grc/check_corner.py` | Guards, exact drift and literal first failures for the seventeen-digit relative corner at two bases |
| `grc_reconstruction/instances.py` | Independent branch compiler and lattice replay, including tied queries, [G4] and [G8] exact exits, and the distinct centred endpoints |
| `grc_reconstruction/physical.py` | Independent literal isolated-wall and all-residue seventeen-digit corner replay at fresh bases |
| `recompute.py` | Release fidelity checks of the 1069 rows and the prime-square lift, the [G4] scalar exit index, all corner residues at base 130243, and a nine-digit chamber fixed point |
| `check_rows.py` | Orbit coefficients of all secondary rows at the prime and its square, with the zero indices 534, 1070 and 1603 at the square |
| `check_evidence.py` | Seventeen-digit finite sample totals from the archived panels and held-out count table |

`SHA256SUMS` pins the executable inputs and the finite evidence. The replay totals the held-out count cells, with 2,048 originals per base. `KER17-PANEL-DIGESTS.tsv` pins the omitted per-base trajectory streams by SHA-256. The command checks this digest table and reads the archived counts. All required inputs are below 2 MB.

The scripts reproduce the finite checks used by this paper. The published results cited as inputs, including the nine-digit physical limiting support law, have their own proof supplements.
