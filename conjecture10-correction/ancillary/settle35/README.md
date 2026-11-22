# SETTLE-35

An exact 12-stage root certificate and guarded transport give, for every integer \(B\ge232892391600\),
\[
\Pr_{35,B}\{ Z_B>1/396\}\ge \frac{35!}{(232892391600)^{17}38^{18}}.
\]
Every constructed terminal lies in the original paired chamber. The probability concerns uniformly sampled nonconstant ordered words. Length 35 is non-strict because \(d_2=0\) on all invariant root-equal profiles. This result refutes the necessity and every-cut decay claims in P10's Conjecture 10 (Strictness and survival) at length 35. [Root certificate and theorem](FAMILY.md) · [Closed update and transport](TRANSPORT.md)

The twelve independent return directions give \(L^{12}\) roots for every integer \(L\ge1\). Lifetime spread is at most \(1334114L\), total export at most \(2940551L\), and the clock through one root period at most \(3785731L\). The terminal relative period is 11. The physical period is a multiple of 11 and divides 66. Literal subtraction verified the complete histories and every phase of three period-66 lifted controls at bases \(232892391600\), \(232892391601\) and \(349338587417\). Their least entries were \(4549628\), \(4853696\) and \(6843992\). [Exact root-stage record](family-certificate.json) · [Literal lift records](literal-lifts.tsv)

## Frozen finite samples

The finite diagnostics used separate independent-root and physical-box denominators, fixed cuts, forward caps and a finite decline criterion. [Frozen design](FORECAST.md) · [Exact tallies](finite-summary.txt)

| Root cube half-width \(L\) | Draws | Repeat censored | Positive recurrent export | Budgeted zero-export repeats |
|---:|---:|---:|---:|---:|
| 64 | 96 | 94 | 1 | 1 |
| 128 | 96 | 95 | 0 | 1 |
| 256 | 96 | 94 | 0 | 2 |
| 512 | 96 | 78 | 11 | 7 |
| 1024 | 96 | 67 | 21 | 8 |
| 2048 | 96 | 0 | 90 | 6 |
| 4096 | 96 | 0 | 86 | 10 |

The low pool's share floor is \(4/288\) and the high pool's share ceiling is \(91/288\) under the specified censor allocations. The finite decline criterion was not met.

| Physical-box base | \(L\) | Draws | Literal cycle censored | Cleared both cuts in \(C\) | Enlarged-chamber leaves / returns |
|---:|---:|---:|---:|---:|---:|
| 200000 | 20 | 24 | 24 | 0 confirmed | 5 / 5 |
| 2000003 | 200 | 24 | 24 | 0 confirmed | 1 / 1 |
| 20000000 | 2000 | 24 | 22 | 2 | 0 / 0 |

Every enlarged-chamber exit in this panel came from a complementary row. All projected updates and loss labels agreed with literal subtraction up to the first exit or the recorded horizon, for 25,352,729 transitions. Pairing and height held at every observed phase. Both hits have least period 66, zero-export root period 11 and no enlarged-chamber exit. Their exact peak numerators over \(18B\) are 9584668 and 9633448, with entry times 76225 and 50416. The per-label exports to initial positions \(5,7,10,14,15\) are \((21170,0,21480,0,22240)\) and \((16066,0,16587,0,16844)\). The other 22 projected roots at the largest scale showed positive recurrent export. Their literal terminal cycles remain censored.

The whole source-box original weights were close to \(4.23\cdot10^{-46}\), \(2.90\cdot10^{-46}\) and \(2.79\cdot10^{-46}\) at the tested bases. The hit fractions are conditional uniform gap-sample counts. The archived \(L=20000\) results form an external comparison, and the all-base theorem uses the explicit root family.

## Lengths 97 and 119

The entire sampled boxes have positive full-height pairing and positive margins on every row with a nonzero centre margin. The rotated complementary row ranks are 0 and 23, respectively. At 97 the enlarged guard exempts the root-only \(d_2\). At 119 the diagnostic exempts the pinned \(d_8\) while tracking every other guard. Both exact centres returned by literal subtraction as off-centre fixed points in \(C\) and passed \(1/1000\). [Frozen design](COMPARISON-FREEZE.md) · [Exact geometry and original source weights](comparison_geometry.json)

| Length | Random draws | Literal cycles resolved | Censored | Detected \(C\) hits above \(1/1000\) | Other-row leaves / returns | Entire box original weight |
|---:|---:|---:|---:|---:|---:|---:|
| 97 | 24 | 3 | 21 | 0 | 22 / 22 | \((4.4180\text{–}4.4182)\cdot10^{-205}\) |
| 119 | 24 | 0 | 24 | 0 | 24 / 24 | \((7.0645\text{–}7.0648)\cdot10^{-247}\) |

The 97 box has 12 initial states outside \(C\) on \(d_2\), and the 119 box has 14 outside \(C\) on \(d_8\). The three resolved 97 cycles are exact central fixed points in \(C\). The 119 paths are all censored. The random panels yielded no finite replication hit. Original-population transfer remains open at both lengths. [97 literal records](comparison-97.tsv) · [119 literal records](comparison-119.tsv)

## Reproduction

Copy these files to a private directory and use standard-library Python and C++17 with assertions enabled. Each substantial run takes an exclusive flock on the project's absolute .task-monitor/compute.lock. Release it between parameter cells.

    c++ -O3 -std=c++17 -Wall -Wextra -Werror probe.cpp -o probe
    c++ -O3 -std=c++17 -Wall -Wextra -Werror retirement.cpp -o retirement
    c++ -O3 -std=c++17 -Wall -Wextra -Werror literal_extra.cpp -o literal_extra
    ./retirement 3500 16384 350245 4194304 5 | cmp - family-route.txt
    python3 -B verify_family.py family-route.txt | cmp - family-certificate.json
    python3 -B comparison_geometry.py | cmp - comparison_geometry.json
    python3 -B summarize.py | cmp - finite-summary.txt

Each cube output is produced by **./probe cube L**, each 35 box output by **./probe box35 B**, and the other comparison outputs by **./literal_extra 97** and **./literal_extra 119**. The lift run uses the following pipeline.

    python3 -B lift_inputs.py family-certificate.json | ./literal_extra lifts

Root columns give least entry and period, detection steps, maximum spread, export through entry or the censor horizon, recurrent export, maximum step export, budgeted hit and initial relative coordinates. Literal columns give least entry and period, detection steps, peak numerator over \(mB\), all-cycle \(C\) and widened-guard flags, first height, pairing and other-row failures, first widened-guard leave and return, first admitted row breach and return to \(C\), first \(C\) leave and return, last compared time, minimum pairing margin over the compared prefix, and row minima. At 35, five per-label exports precede the initial gaps. Elsewhere only the gaps follow. Unobserved times use \(-1\), and a zero period encodes a censored repeat. [Probe reader](probe.cpp) · [Lift check](literal_extra.cpp)
