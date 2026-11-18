# Reproducing the nine-digit Kaprekar law

Contact: Rustem Galyavin, rustem@galyavin.com. The paper and data are CC BY 4.0 and the code is MIT. See [LICENSE.md](LICENSE.md).

Use Python 3.9 with the pinned packages below, a C compiler named `cc`, a C++17 compiler named `c++`, and Yices 2.7.0 with `yices-smt2` on the search path. Z3 and cvc5 use their Python bindings. Install Yices for your operating system from [SRI](https://yices.csl.sri.com/). On a machine with those prerequisites, run these commands from this directory.

```sh
python3.9 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
bash reproduce.sh
```

The default runs all eight stages. It runs the original encodings and the independent review in disposable working directories. An assertion failure, solver timeout, missing child, unrecognized status or rejected control makes the relevant stage fail. For a permanent terminal log, redirect stdout and stderr and check the command's exit status. The clock audit uses four worker processes. The full run includes long exact route recompilations and millions of individual subtraction steps.

| Command argument | Calculation | Evidence for |
|---|---|---|
| `entry` | Reconstruct the four real implications with Z3 5.4.1 and cvc5 1.4.2, then compare complete integer trajectories and raw real folds | Finite periodicity and almost-everywhere strict entry |
| `routes` | Rebuild each sorted word from subtraction, verify complete phase sets, reconstruct every guard, compare first-cycle states against intermediate iterates | Bulk, faces, corners, releases, endpoints and first capture |
| `weights` | Recompile the 33 routes in all five base classes from exact integer lattice images | All 165 rational kernels, including the index-eighteen line and the displayed separation |
| `covariance` | Reconstruct the complete paired-order inventory, check preservation under translation, and follow full matched orbits | Applicability of the all-odd analytic identity on the nine-digit core |
| `zeros` | Check the positive first-zero cutoff, sixteen zero-return implications and five one-step inequalities, using Z3 and cvc5; compare complete finite graphs | Linear zero preparation and the dispersion increment |
| `descent` | Rebuild the 240 root cases and all 29,534 affine nodes, requiring UNSAT in both Z3 and Yices. Replay primitive columns, independent cycles and the obstruction at eight steps | Uniform linear capture time. This argument first runs the zero and finite-graph prerequisites when selected individually |
| `fresh` | Subtract digit columns and find the first cycle by direct iteration, then classify the realized time and terminal state against the four separating-route predictions | Independent finite test of the five conditional laws |
| `independent` | Compile the separate C model, rerun exact checks and Z3 queries, enumerate descent states, simulate the conditional and uniform populations and reconstruct their comparisons | Independent verification of the reductions and finite-base predictions |

For example, `./reproduce.sh entry` runs the real-entry stage alone. The driver uses the local `.venv` if available. Set `PYTHON` to select another interpreter.

## Reading the kernel

`source/nine_kernel_results.json` orders the rows by route and then residue. A route names the successive decreasing facets and the final channel. The 5-entry clock row and each row of the 4-by-5 terminal matrix act on the column consisting of the four normalized core differences followed by one. The weight is a rational number. To evaluate the law at a generic original input, use the rational evaluator:

```sh
cd source
../.venv/bin/python - <<'PY'
from fractions import Fraction as F
from nine_global import limit_kernel
x = tuple(F(i, 50) for i in (43, 39, 31, 28))
for beta in range(5):
    result = limit_kernel(x, beta)
    print(beta, result['route'], [(a['time'], a['weight']) for a in result['atoms']])
PY
```

Finite-cap evaluators raise an error on ties, unresolved entry, or a failed finite-base guard. The proof takes compact subsets of strict entry cells and then exhausts their limiting input mass. The fresh numerical check uses a box strictly inside one such cell. The universal clock applies to every input, including ones outside those cells.

The finite core routines report the first cyclic difference vector. For a transient initial difference state the first cyclic *number* is one update later, and the global evaluator returns its differences separately as `first_cyclic_number_differences`. The entire limiting terminal segment is pointwise fixed under the leading map, as proved in the paper.

The fresh runs reconstructed all 165 kernels, passed every descent node in both solvers, and matched every direct observation at bases 300035–300039 to a unique joint atom. See [the fresh-copy record](RERUN.md) for the environment, check counts and observed frequencies. The analytic proofs of reflection, convergence and the conditional period restriction are in the paper.

## Environment and runtime

The reference machine uses macOS 27.2.1 (arm64), Python 3.9.6, SymPy 1.14.0, mpmath 1.3.0, Z3 5.4.1 (Python package 5.4.1.0), cvc5 1.4.2, Yices 2.7.0, NumPy 2.0.2, SciPy 1.13.1 and Apple Clang 21.0.0 for C and C++. See [the measured stage times and fresh results](results/release-rerun.md).

The independent sample comprises 2,586,000 conditional orbits and 1,600,000 uniform strings. It also checks 811,157,145 qualifying states exhaustively at 65 bases, and 3,132,549 sampled descent states. The expected boundary and mutation diagnostics and the list of omitted raw files are in [independent-review/README.md](independent-review/README.md).

The [SHA-256 manifest](results/artifact_sha256.json) lists every distributed file except the manifest itself. From this directory, verify the original files before running the computations:

```sh
python3 - <<'PYCODE'
from hashlib import sha256
from pathlib import Path
import json
for name, digest in json.loads(Path("results/artifact_sha256.json").read_text()).items():
    if sha256(Path(name).read_bytes()).hexdigest() != digest:
        raise SystemExit(f"Checksum mismatch for {name}")
print("Checksums match.")
PYCODE
```

Measured wall times (seconds) on the reference machine: entry 295, routes 90, weights 130, covariance 4, zeros 75, descent 362, fresh 42, independent 2572. The stage-by-stage total was 0.99 h.
