# Cycles of the four-digit Kaprekar map in even bases

Contact: rustem@galyavin.com

Four digits, including leading zeroes, are sorted in descending and ascending order, then the smaller number is subtracted from the larger. The checks enumerate the cycles reached from nonconstant inputs. Their state is the pair of outer and inner differences of the sorted digits.

The text, figures and data are licensed under CC BY 4.0 and the code is licensed under the MIT License. [License terms](LICENSE.md)

## Reproduce the numerical results

From this directory run

```sh
bash reproduce.sh
```

Requires Bash, Python 3.7 or newer, and a C++17 compiler. The script uses the Python standard library. Set `PYTHON` or `CXX` to choose executables. It verifies the SHA-256 manifest, builds the C++ program in a temporary directory and compares both regenerated JSON files byte for byte with the distributed results.

| Check | Expected result |
|---|---|
| Digit examples | 6174 has pair (6,2) and is fixed. In base 20 the listed fixed point and six-cycle close in order. In base 56 the listed three-cycle and four-cycle close in order. |
| Boundary table | Every timed first return agrees in 120 bases with odd part 3 through 41 and exponent 1 through 6. Of 138,240 boundary states, 18,540 return. The paper's range through odd part 31 and exponent 4 accounts for 18,945 tested starts. |
| Critical states and Figure 1 | The set of observed targets equals the critical set at each base. Of 1,200 critical starts, 147 return to S, E or U after exactly n+1 steps. The other 1,053 enter a boundary-free periodic orbit. All 360 decisions for S, E and U agree with the figure. |
| Even bases through 1000 | All 500 classifications agree on 83,708,250 difference pairs. For odd parts at least five, every directed lift, every additional cycle representative and its least period are checked. For odd parts one and three, all period multiplicities are checked. The uniqueness criterion agrees throughout. |
| Powers of two | The multiplicities of Table 1 agree at every exponent from 1 through 12. At base 4096 there is one cycle of each length 6 and 7, and five of each length 12 and 14. |

Both [the Python boundary check](rules.py) and [the C++ cycle enumeration](cycles.cpp) perform four-place subtraction with borrows. A boundary-free repetition proves absence of a return. The numerical results for each base are in [boundary.json](results/boundary.json) and [cycles.json](results/cycles.json). The [SHA-256 manifest](results/artifact_sha256.json) covers all distributed files except itself.

## Build the formal proof

With Lean 4.28.0 and Lake available, run

```sh
cd lean
lake build
lake env lean TheoremC.lean
```

The source uses the libraries distributed with Lean. A successful build finishes 16 jobs. The second command reports `[propext, Classical.choice, Quot.sound]` for both `KaprekarBoundary.every_even_base` and `KaprekarBoundary.theoremC`. The [Lean guide](lean/README.md) maps the formal statements to the paper.

On 17 November 2026, a fresh copy on macOS 27.2.1 (arm64) completed the numerical command in 1.61 seconds using Python 3.9.6 and Apple Clang 21.0.0. The 26 file hashes and both regenerated results matched. With Lean 4.28.0 and Lake 5.0.0, a clean build took 13.97 seconds and the axiom check took 0.72 seconds.
