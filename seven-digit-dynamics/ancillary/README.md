# Seven-digit Kaprekar dynamics

Reproduction files for Rustem Galyavin, *Seven-digit Kaprekar dynamics: global law, two clocks, parity*. Contact: rustem@galyavin.com.

Run from this directory in a fresh copy.

```sh
bash reproduce.sh
```

The command needs Python 3.9 or later, Bash, a C++17 compiler named `c++`, and Yices with `yices-smt2` on the executable path. It uses an installed Z3, cvc5 and SymPy if available. Otherwise it creates a temporary Python virtual environment and installs the pinned packages in `requirements.txt` using pip. Set `PYTHON` to select an interpreter, including one with the pinned packages already installed for an offline run. Keep Python assertions enabled. The checked system used Python 3.9.6, Z3 5.4.1 (Python package 5.4.1.0), cvc5 1.4.2, Yices 2.7.0, SymPy 1.14.0 and Apple clang 21.0.0. Wall time was about five minutes on the checked machine. The descent uses four worker processes.

The script verifies the SHA-256 manifest, copies the checkers to a temporary directory, regenerates the solver queries and checks the results. It prints every stage and exits unsuccessfully on any failure. Intermediate solver inputs and timing files are temporary. To retain the console record, run `bash reproduce.sh > replay.log 2>&1`. The supplied `results/reference.log` records the original ten-stage, 268-second run. The release command also runs the independent transition-model self-check. The final fast residue table is checked in full against the computed rates.

| Check | Expected result |
|---|---|
| Leading-map entry | Eight real-linear counterexample queries UNSAT in both Z3 and cvc5 |
| First access to zero | A three-step positive prefix ending at zero is UNSAT in both engines for every integer base at least 18 |
| Returns from zero | Seven integer formulas UNSAT in both engines, six explicit clearance controls |
| Potential | Seventeen integer queries UNSAT in Z3 and Yices. Two false strengthenings SAT, with replayed paths |
| Edge | 884 even and 845 odd reachable offset states, two unbounded guard checks in Z3, 140 lifted terminal families; error in the interval [−43/2, 343/8] |
| Facet | 338 restricted-source and 490 two-sided rational guard identities. Every printed column of both phase tables agrees. Fast fractions for positive/negative sources are 1/2 and 3/4 in even bases, 3/8 in both odd-base sources |
| Direct numerical replay | 952 bulk exits, 288 rotated/predecessor checks, 12 complete leading-prefix clocks over 1,643,929 digit steps; complete populations in bases 96, 97, 192 and 193 |
| Independent transition model | Direct integer subtraction agrees with the borrow model on all nonconstant difference states for bases 2 through 40 and on 20,000 random states at larger base bounds |
| Independent C++ graphs | 22,315,002 difference states in bases 129, 130, 400 and 401. All nine-step prefixes tested in the first two bases |

The final line reads `All seven-digit checks passed.` The finite populations and sample integrals check implementations. The universal descent, entry inequalities and lifted edge paths, together with the arguments in the paper, establish the asymptotic assertions. See the appendix for each finite claim and its encoding.

### Independent verification files

The [independent report](independent-review/independent-verification.md) identifies the domains and outcomes of the separately written checks, including the partial Yices result on nine-step descent. The `independent-review/` directory contains the Python and C encodings, exported SMT-LIB inputs, clock-box and exhaustive logs, and the numerical summary `logs_analyze_fresh.txt`. The independent clock boxes used bases 200000–200007, 1000003 and 1000004. The exhaustive cyclic-radius measurements used bases 18–300, 400, 401, 500–507, 1000 and 1001.

The one-command replay checks the shipped hashes, runs the ten original checker programs in the order printed, and ends with the independent transition-model test. The independent solver and population campaigns have separate recorded logs. To rerun the independent numerical campaign, copy `independent-review/` to a working directory, then run the following commands there. Use a C compiler with 128-bit integer support and Python with SymPy 1.14.0 and NumPy 2.0.2. The recorded campaign took 252 seconds on the review machine. The optional Yices partition drivers also require `timeout`.

```sh
cc -O2 exhaust.c -o exhaust -lm
cc -O2 traj.c -o traj -lm
bash campaign_fresh.sh
"${PYTHON:-python3}" analyze_fresh.py
```

Six one-arrival-per-word sample files are omitted: `logs/pop_20011.txt`, `logs/pop_20012.txt`, `logs/pop_20013.txt`, `logs/pop_20014.txt` (100,000 rows each), `logs/pop_100005.txt` and `logs/pop_100006.txt` (20,000 each). The campaign regenerates them from their recorded seeds. Their original summary is in `logs_analyze_fresh.txt`. The two platform-specific executables are supplied as C source. The machine-local scheduling wrapper is omitted.

Text, figures and data are licensed under CC BY 4.0. Code is licensed under the MIT License. See [LICENSE.md](LICENSE.md).
