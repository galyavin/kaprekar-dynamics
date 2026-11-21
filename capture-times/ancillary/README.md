# Reproduce the exact computations for P12

From a fresh copy of this directory, run:

```sh
bash reproduce.sh
```

The script requires Bash and CPython 3.9 or newer, with only the Python
standard library. Set `PYTHON` to choose an interpreter. The script enables
Python assertions, checks the bundled SHA-256 manifest and compares the full
output of each checker against the frozen output. The order-15 replay checks
over 11 million transitions. The checks write temporary output outside this
directory and remove it at exit.

| Directory | Inputs and computation |
| --- | --- |
| `recon14_physical` | The released off-centre cycles paper and its ancillary ZIP supply the strict profiles and retirement records for orders 3–13. The checker reads three JSONL entries from the ZIP and executes none of the archived scripts. It also checks literal small-base controls for the paired update, inverse bound, breach and digit replication. |
| `recon15` | Exact finite controls for rank allowances, side-word counts, one-step reflected siblings, entrance and decay inequalities, clock-tail geometry, read recurrence, and physical guard translations. |
| `recon_carriers` | Transfer counts, exact graph and ternary word enumerations on finite ranges, alternative entry checks, six-block checks and explicit cycles of periods outside the graph. The arbitrary-order graph and count follow from the paper's proof. |
| `recon_root15` | Two frozen JSONL files encode four reference sources, their insertion intervals and their proposed read-stage data. The checker independently replays their paths and reconstructs the determinants, all-time coefficients, and exact folded-order and strictness certificates for the twelve lengths. |

The two files in `recon14_physical` named after the off-centre cycles paper
are byte-for-byte copies of the released P9 source and ancillary archive.
The two JSONL files in `recon_root15` are the frozen order-15 inputs. The
manifest pins these inputs, the scripts and the expected outputs.
[P9 release](https://doi.org/10.5281/zenodo.24399342)

All inputs needed by the four commands are bundled and below 2 MB. The
commands run offline using exact integer and rational arithmetic. The paper
gives the analytic arguments for arbitrary bases and lengths.
