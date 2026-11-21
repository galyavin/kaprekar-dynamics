# Reproducing the arithmetic census

Run `bash ancillary/reproduce.sh` from the paper folder or `bash reproduce.sh` from this folder. The bundle needs Python 3.9 or newer, the standard shell tools (including `shasum`), a C++17 compiler, and a C compiler supporting `unsigned __int128` (Clang or GCC). The independent scan uses about 1 GB of RAM. The script prints the interpreter and compiler versions and runs in a disposable directory for generated files.

The run checks SHA-256 pins, replays the census and five supplemental certificates, verifies the four census totals and the overlap with the earlier census, and reruns the independent prime, composite and exact prefix checks. It compiles and runs both exhaustive index-two scans below 10⁹. It compares deterministic output against the frozen logs, omitting two timing lines from the census log and one from the independent scan log. No solver, network retrieval or external data file is required for replay. All required inputs are smaller than 2 MB. The manifest pins every file in the bundle except the manifest itself.

In a shared research checkout, set `P11_COMPUTE_LOCK` to the absolute path of the existing project compute lock before running the command. The replay takes that lock for the duration of the run.

| Folder | Source and purpose |
| --- | --- |
| `arith/` | Task 207 primary census certificates, supplemental planar certificates, and verifier. The bundle uses Python integer modular elimination at the same prime, 65521, in place of the NumPy implementation of that step. The certificate data and arithmetic conditions are unchanged. |
| `arith2/` | Task 209 independent index-three and index-four census and prefix check. |
| `recon14_arith/` | Task 205 independent prime walk, index-two identities and finite prefix checks, plus the exact order-of-4 and Legendre-row scan. |
| `recon14_lowt/` | Task 208 independent composite and physical-subtraction checks. |
| `recon14_sieve/` | Task 212 independent exact index-three and index-four feasibility and prefix replay. |
| `independent/` | Fourteenth-compass frame check, including the earlier census and the independent order-of-2 and literal-slot scan. |

The frozen certificate replay is separate from heuristic certificate discovery. The copied census source also retains the original candidate generator, which was run with NumPy and SciPy (HiGHS). The one-command replay checks each proposed certificate exactly and does not invoke the generator. Fixed-prefix counts and seeded samples in the logs serve as finite controls for the analytic argument.
