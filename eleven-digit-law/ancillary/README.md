# Eleven-digit Kaprekar arrival laws modulo twelve

Run this command from a fresh copy of the directory.

```sh
bash reproduce.sh
```

Use Python 3.9 or newer with SymPy, a C++17 compiler available as both
`clang++` and `c++`, Yices (`yices-smt2`), and `shasum`. The verification
environment used Python 3.9.6, SymPy 1.14.0, Apple Clang 21.0.0 and Yices
2.7.0. The inherited results cited in the paper, including the centered
clock bound, have their own certificates in the cited releases.

The command checks the SHA-256 of the supplied source and reference files
in `inputs/` and works on a private copy. The 9.7 MB box leaf stream is
generated from its pinned inputs and verified against `LARGE.sha256`.
Generated files and executables live in temporary directories. The default
run includes exhaustive geometric and arithmetic replays, as well as the
frozen physical panels. Several checks run six or eight independent roots
or sample shards in parallel. The component replays together took roughly two hours on the
verification machine, apart from waits for the shared compute lock.

Set `PYTHON` to select an interpreter, for example a virtual environment
with SymPy. Set `P8_COMPUTE_LOCK` to an absolute existing lock path when
sharing compute with another campaign. The runner otherwise makes a
private lock.

| Paper result | Replay |
| --- | --- |
| Exact finite route and original-lattice controls, paired cube and box | `bash reproduce.sh --part local` |
| Basin exclusions, entrance through step four, folded boxes, both anchor and three-step integrations | `bash reproduce.sh --part mass` |
| Primary and independent six-root rank and collar coverage, four-lift first-coordinate bound | `bash reproduce.sh --part rank` |
| Original and independently reconstructed arithmetic inventories and terminal period-66 sections | `bash reproduce.sh --part period` |
| Frozen physical diagnostics and prospective evidence | `bash reproduce.sh --part physical` |

`--part quick` runs the local route checks and the smaller basin and
arithmetic calculations, while reading the pinned three-step integration
ledger. `--part mass` regenerates that ledger by both integrations.
For the small cube and box the symbolic checks use exact rational weights.
The original-word sample reports use floating-point arithmetic for
estimation and intervals. For the prospective panels, the raw shard hashes
and integer counts must agree exactly. Floating-point summary values are
compared with an absolute tolerance of \(10^{-13}\) to allow for differences
between Python's summation implementations. The empirical checks have no role in the proof
of the joint limit or the mass lower bound.
