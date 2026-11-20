# Reproducing the off-centre Kaprekar certificates

From this directory in a fresh copy, run:

    bash reproduce.sh

The replay verifies SHA-256 for every input and the two replay programs in
[SHA256SUMS](SHA256SUMS). It copies the inputs into a private working directory,
runs with assertions enabled, and removes that directory on exit. No network
access, external solver or nonstandard Python module is needed. The tested
configurations used Python 3.13.15 and 3.9.6 with Apple clang 21.0.0 in C++17 mode.
Set `PYTHON` or `CXX` if those executables have different names locally.
For a replay inside the research checkout, set `P9_COMPUTE_LOCK` to the
absolute path of its shared compute lock. A fresh standalone copy uses a
private lock.

| Paper result | Replayed inputs under the archived paths |
|---|---|
| Period eight at 15 and 17 | Tasks 155 and 157 sign-word verifiers, physical profile checks and literal word orbits |
| Strict lengths and retirement | Tasks 164, 168 and 170 root certificates, order-12 itinerary, divisor screen, rational lifts and triangular arrays |
| Periods 24 and 132 | The ordinary order-12 record and the shifted record `period132-root.jsonl`, with the independent lifts at both lengths |
| Full-orbit and forced-tie examples | Literal fixed-point controls at length 23 and Task 164's ten-step feeder and controls at length 25 |
| Thirteen digits | Task 158 six-ray certificate and earlier route premise, Task 167 second source box and outer key, exact polynomial mass, and a comparison of the rescaled terminal cycles |
| Order-eight remarks | Exact component arrays and rational volumes, both physical lifts at 51 digits, and replay of the finite cube evidence |
| Independent reconstruction | Task 177's four replay stages. Their combined output must match the pinned reconstruction log byte for byte |

The highlighted claims use finite exact predicates over affine parameter
families. Literal whole-word trajectories are controls for the reduction
from sorting and subtraction to gaps. The finite parity and mixture
experiments are empirical; the classifiers and decisions are replayed
exactly. The independent reconstruction includes the wall-feeder stage to
reproduce the complete archived run, although the paper omits the standalone
wall-feeder proposition.

All inputs required by this command are supplied. No required input exceeds
2 MB. The bundle omits unrelated population-panel archives. The only archived
outer-section input for the Task 167 primary proof is the compressed
`encounters.tsv.gz`; its referenced key is also independently reconstructed
from the source route. [Fresh-copy replay output](verification.log)
