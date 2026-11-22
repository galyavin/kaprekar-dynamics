# Length 35 retains macroscopic memory

For every integer \(B\ge232892391600\), a uniform nonconstant ordered word of 35 base-\(B\) digits satisfies

\[
\Pr_{35,B}\{\text{terminal in }C,\ Z_B>1/396\}
\ge \frac{35!}{(232892391600)^{17}38^{18}}
=2.170140224479\ldots\times10^{-182}.
\]

Here \(C\) is the original chamber in P10. Every phase of each constructed terminal has radius above the cut. These orbits preserve positive pairing and rows \(d_3,\ldots,d_9\) from the source, and incur transient breaches of \(d_2\). The exact replay and the physical proof below establish the all-base statement. The two direct-subtraction controls provide finite checks at fresh bases.

[Released P10, §§1–2, 6](../../../docs/papers/kaprekar-macroscopic-cycles.pdf) · [Task 245 certificate data](../settle35/family-certificate.json) · [Exact replay](symbolic.log) · [Literal controls](literal.log)

## Derivation from subtraction

Write \(n=2k+1\), \(m=k+1\), and let \(q_i\) be the difference between the \(i\)-th largest and \(i\)-th smallest digit. Use \(q_0=B,\ q_{k+1}=0\) so that \(g_j=q_j-q_{j+1}\) for \(0\le j\le k\). The signed columns, from least to most significant, are \(-q_1,\ldots,-q_k,0,q_k,\ldots,q_1\). If \(q_k>0\), the borrow survives through the central column and stops at \(q_k\). The central column produces \(B-1\). The other digits at level \(i\) are

\[
A_i=q_i-[i=k],\qquad D_i=B-q_{k+1-i}-[i<k].
\]

Positive pairing requires full height \(q_k>0\) and \(\min(A_i,D_i)\ge\max(A_{i+1},D_{i+1})\) at every adjacent level, where ties are allowed. The digit \(B-1\) is largest. The two same-column differences across a level boundary are a single gap plus an offset of magnitude at most 1. For \(A_i-D_{i+1}\) the gap contribution is the sum of the suffixes starting at \(i\) and \(k-i\), minus \(B\). For \(D_i-A_{i+1}\) the contribution is the sum of the prefixes of lengths \(k+1-i\) and \(i+1\), minus \(B\). Both sums use \(m+1\) gap terms, counted with multiplicity, and the offsets are at least \(-1\). Thus \((m+1)\min_j g_j\ge B+1\) implies the exact adjacent-level guard and full height.

The P10 row comparisons are \(d_i=\sum_{a=0}^{i-1}(g_{k-a}-g_a)=q_i+q_{k+1-i}-B\) for \(1\le i\le\lfloor m/2\rfloor\), with \(d_0=0\). Set \(h_i=(-1-d_i)_+\) and \(h_0=0\). In a provisional calculation use \(q_i,B-q_{k+1-i}-1\) at every level, and arrange their maximum \(H_i\) then their minimum \(L_i\) in the required level order. This provisional list only serves to compute differences at the prescribed positions. Its within-level spread is \(t_i=|1+q_i+q_{k+1-i}-B|\), with \(t_0=1\). Its ordered positions give \(q'_{2i-1}=L_{i-1}-L_{k+1-i}\) and \(q'_{2i}=H_i-H_{k+1-i}\), with the boundary convention \(L_0=B-1,\ L_{k+1}=-1\). Mirror levels have the same spread. Consequently the even gap \(g'_{2i-2}\) equals \(L_{k+1-i}-L_{k+2-i}\), and the odd gap \(g'_{2i-1}\) equals \(L_{i-1}-L_i\). Their sum is \(g_{i-1}+g_{m-i}\), by \(L_{k+1-i}-L_i=q_{k+1-i}-q_i\), and their even-minus-odd difference is \(t_{i-1}-t_i\), by \(L_{k+1-i}+L_i=B-1-t_i\). For odd \(m\) the last unmatched gap is \(g_{k/2}\). The identities also cover the ends of the list by the boundary convention. Inserting \(t_i=1+d_i+2h_i\) gives the provisional folded permutation \(P\), sending \(j\) to \(\min(2j,n-2j)\), plus

\[
(a_{2i-2},a_{2i-1})=(-h_i+h_{i-1},h_i-h_{i-1}).
\]

Restoring the last level changes its minimum by
\(\min(g_k-1,g_0)-\min(g_k,g_0-1)=e=\operatorname{sgn}(d_1)\) on integers. The maximum changes by \(-e\), as the level sum is unchanged. The first two symmetric differences change by \(-e,e\). Their gap correction is \(e(1,-2,1,0,\ldots)\). This yields the exact update on every positively paired, full-height physical state, including ties in the adjacent-level guard.

Let \(J\) be the orbit of 1 plus slot 0. The images under \(P\) of the support of row \(i\) are precisely slots \(0,\ldots,2i-1\). The rows supported wholly on \(J\) therefore form an initial block \(1,\ldots,I\). At a proper folded length, \(I<\lfloor m/2\rfloor\). Otherwise every slot is in \(J\), including the possible unmatched slot \(k\), which maps to 1. Slot \(2I+1\) is the first outside \(J\), since \(2I\) is the image of the root slot \(I\). If the complementary rows all satisfy \(d_i\ge-1\), their \(h_i\) vanish and each remaining \(h_i\) contributes

\[
h_i(-\mathbf e_{2i-2}+\mathbf e_{2i-1}+\mathbf e_{2i}-\mathbf e_{2i+1}).
\]

All correction slots lie in \(J\) except \(2I+1\). The root sum increases by \(h_I\). For \(f(j)=\min(2j,n-2j)\), a complementary label starting at \(j\) is at \(f^t(j)\) at time \(t\); its value there is \(g_j(0)-\sum_{0\le s<t,\ f^{s+1}(j)=2I+1}h_I(s)\). All root reads are unchanged by adding a constant to the root gaps, so subtracting the new anchor gives a closed relative map independent of \(B\) and the complement.

At 35 the fold order is \((1,2,4,8,16,3,6,12,11,13,9,17)\), the complementary orbits are \((5,10,15)\) and \((7,14)\), and \(I=2\). Each step charges the label arriving in slot 5, with first arrivals from the starting slots 15, 10, 5 in that order. The labels from 7 and 14 retain their values. For root offsets \(x_1,\ldots,x_{12}\) relative to slot 0, \(d_1=x_{12}\) and \(d_2=x_{12}+x_5-x_1\). The update right-rotates \(x\) and adds the following increments, where \(a=h_1-e\).

| Index | Increment |
|---|---:|
| 1 | \(2h_1-3e\) |
| 2 | \(2h_1-h_2\) |
| 3 | \(a+h_2\) |
| 6 | \(h_2-e\) |
| All others | \(a\) |

## Symbolic access, injectivity and the clock

The independently written [checker](reconstruct.py) uses the frozen certificate as numerical input. It traces the root in physical fold order using the kick vectors obtained above. In moving label order at phase \(t\bmod12\), assign fold index \(j\) (one-based) label \(t+12-j\bmod12\). The two phase reads are then \(f_t(b)=b_t\) and \(f_{12+t}(b)=b_t+b_{t+7}-b_{t-1}\), with label subscripts modulo 12.

At stage \(j\) a reference prefix of length \(p_j\) and a reference word of length \(w_j\) evaluate both reads in the current span \(V_j\) at their exact values. Reads outside \(V_j\) take the forced primary positive or secondary zero-export branch, respectively. Two inputs differing by a direction annihilated by \(V_j\) receive identical kicks and export throughout a projected word. The exact replay derives the prefix endpoint \(a_j\), word endpoint \(b_j\), return ray \(u_j=a_j-b_j\) in moving labels, and the next reference start \(s_{j+1}=b_j+v_j\), with certificate shift \(v_j\). It checks both \(u_j\) and \(v_j\) annihilate \(V_j\). All prefix and word lengths are multiples of 12.

One read form is inserted in the span at each stage using rational elimination. The checker derives the newly tracked set \(G_j\), which can contain several forms, and verifies that every \(f\in G_j\) has \(f(u_j)>0\). It also checks the tracked and untracked read values against the complete prefix and word minimum tables. The forms inserted in order are \(8,4,10,21,20,16,22,0,13,14,18,7\). The matrix with those forms as rows evaluated on the return rays as columns has zeros above the diagonal since \(u_j\) annihilates all earlier forms. Its positive diagonal product is \(5779085723146507793742893652480\). In particular, the rays are independent.

For integer \(L\ge1\) choose independently \(K_jL\le c_j<(K_j+1)L\) and start with \(s_0+\sum_j(v_j+c_ju_j)\) in moving labels at phase 0. Replay the stage-\(j\) prefix once and the word \(c_j+1\) times. In phase-aligned label coordinates the offset from the stage reference is

\[
v_j+b u_j+\sum_{\ell>j}(v_\ell+c_\ell u_\ell),
\]

with \(b=c_j\) in the prefix and \(b=c_j,c_j-1,\ldots,0\) in the word copies. At a completed word the offset relative to its start decreases by \(u_j\). The final copy reaches \(s_{j+1}+\sum_{\ell>j}(v_\ell+c_\ell u_\ell)\) at phase zero. Thus the stages compose, provided that every forced read is on its assigned branch.

The replay computes the exact minimum reference read for each form on each occupied segment and determines the integer floors \((0,0,0,0,0,15,0,1,1,11,1,1)\) from sufficient branch inequalities. It then verifies every forced comparison by independent endpoint substitution. A form installed in stage \(i\) annihilates all shifts and rays from stages after \(i\). For a prefix, its current repeat count ranges from \(K_jL\) to \((K_j+1)L-1\); for a word copy the remaining count ranges from 0 to \((K_j+1)L-1\). Counts in later stages have the prefix range. Each read on that segment equals its reference value plus a fixed linear form in these counts. For each coefficient the minimum occurs at the lower endpoint if nonnegative and at the upper endpoint if negative. The resulting least slack above 1 (untracked primary) or \(-1\) (untracked secondary) is \(A L+C\) with integers \(A,C\). All 343 occupied segment/form inequalities have \(A\ge0\) and \(A+C\ge0\). The equality of tracked reads follows at every phase from the annihilator condition, so this checks ties as well as strict branches on the full parameter range. The directions' independence gives exactly \(L^{12}\) distinct integer starts.

After the twelve stages the common release root is \((-1,0,2,2,0,11,-1,2,7,2,0,-1)\). Exact first-repeat detection finds entry after 30,428 steps, with export 26,526 before that entry. Its eleven cyclic states are distinct and all have \(h_2=0\). The replay also gives these bounds for each inflated orbit.

- The maximum over its lifetime of the root spread, including slot 0, is at most \(1334114L\). This is twice the bound given by the maximum absolute reference coordinate plus the summed sup-norm bounds for the inserted translations. It is at least \(1+77454L\), from primary read 7 in the first copy of the last word. The final copy has read at least 1 and the last return ray has read 77454.
- The total export lies between \(461994L\) and \(2940551L\). The first bound comes from at least \(L\) copies of the last word, each exporting 461,994. The upper bound sums each fixed prefix, at most \((K_j+2)L\) copies per word, and the release-to-entry export.
- The first zero-export recurrence is at or after \(596448L\), since the last word has length 596,448, has positive export, and at least \(L\) copies precede its last copy. One full relative period after release is complete by \(3785731L\), summing step lengths as for the upper export bound and including the eleven cyclic steps.

## Physical transport at every base

At 35 set \(\pi_j=3/55\) on \(J\) and \(\pi_j=16/275\) otherwise. This profile sums to one and is invariant under folded doubling. Its nine rows are \((0,0,1,2,2,1,1,1,1)/275\). The pairing slack above \(1/19\) is \(2/1045\). The suffix-sum radius of \(\pi\) is \(1/198\) and is attained at \(q_5\). For example, the profile deviations \(\pi_j-1/18\) in units of \(1/4950\) are \(-5\) on \(J\) and \(13\) outside it, and the largest absolute prefix sum is 25, at the first five slots. The normalized centred perturbation belongs to \(W_{35}\). Row \(d_2\) lies entirely on \(J\) and hence vanishes on every root-equal profile. Thus 35 is non-strict by P10's definition, and the other secondary rows can all be positive.

Take \(C_0=2940551\), \(E=18+2C_0+1=5881121\), \(D=116446195800=19800E\) and \(L=\lfloor B/D\rfloor\). For each of the \(L^{12}\) relative roots choose \(g_0\) in the interval of \(L\) consecutive integers starting at \(\lfloor B\pi_0-(\sum x_i)/13\rfloor\). The offsets determine the remaining root gaps. Choose the gaps at \(5,7,10,14\) from \(L\) consecutive integers each starting at \(\lfloor B(16/275)\rfloor\). Slot 15 enforces the exact sum \(B\). The independent choices give exactly \(L^5\) distinct integer lifts of every distinct relative root. The initial root mean deviates by at most \(L\), the four free complementary gaps by at most \(L\) each, and the dependent gap by at most \(17L\) by adding the mean error times 13 and the four free errors.

To certify transport, first generate a candidate path by the algebraic update with \(h_i=0\) for \(i\ge3\). It conserves total sum and follows the verified relative root. At every time the root mean increases from its initial value by cumulative export divided by 13, at most \(C_0L/13\). Every root coordinate differs from its current mean by at most the root spread \(C_0L\). Each complementary gap has cumulative loss at most \(C_0L\) and remains on its own profile level because \(\pi\) is invariant. Thus every root gap is within \((1+C_0+C_0/13)L\) of \(B\pi_j\), and every complementary gap is within \((17+C_0)L\). These are bounded by \(EL\le B/19800\), uniformly over the entire route and all choices.

This bound gives \(g_j\ge(3/55-1/19800)B=1079B/19800>B/19\). For \(B\ge2D\) it gives \(19g_j\ge B+1\) and every slot is positive. The bound for every row \(3\le i\le9\) is

\[
d_i(g)\ge B\bigl(d_i(\pi)-2i/19800\bigr)\ge 3B/1100>0.
\]

The sufficient pairing test proved above therefore verifies full height, every literal adjacent-level comparison and all remaining secondary guards at every candidate state. The integer positive gaps summing to \(B\) are physical by the digit construction below. Induction identifies the candidate path with actual subtraction. At each phase, \(q_5/B-13/18\ge1/198-13/19800>1/396\).

After export ceases, the root mean is constant, so the relative cycle repeats the absolute root every eleven steps. Complementary labels just permute with order dividing six, so the full gap period is a multiple of eleven dividing 66. A gap state determines its successor word by subtraction, while a word determines its own gap state. Thus the gap period and the word period agree on the cycle. On every terminal phase \(h_2=0\) gives \(d_2\ge-1\), while the seven other secondary guards already hold. This places the entire terminal in the original chamber \(C\), where \(h_1\) is allowed. Positive total export forces a breach of \(d_2\) somewhere on every transient route in this construction.

For general proper complement-strict lengths, suppose a relative-root family gives at least \(\alpha L^r\) distinct roots at all sufficiently large integer scales, with lifetime spread and total export both at most \(C_0L\), and eventual periodicity with zero recurrent \(h_I\). Mix an invariant root-equal profile with strict rows \(i>I\) with the uniform profile to secure \(\Delta=\min_j\pi_j-1/(m+1)>0\). Strictness of these rows also makes \(\gamma=\min_{i>I}d_i(\pi)>0\) and the profile radius \(A>0\). The same anchor and complementary floors give \(L^{k-r}\) lifts per root. The initial mean and free gap errors are at most \(L\), and the dependent gap error is at most \(mL\). The mean increases by export divided by \(r+1\), while each complement label loses at most the total export and stays on its own profile level. With \(E=m+2C_0+1\) every gap error is at most \(EL\). For \(D_n\ge4E/\min(\Delta,\gamma/m,A/m)\) and \(L=\lfloor B/D_n\rfloor\), pairing retains at least \(3B\Delta/4\) above the per-gap threshold, the complementary rows retain at least \(3B\gamma/4\), and the radius is at least \(3A/4\). For all sufficiently large \(B\), these margins ensure physicality and the literal pairing guard. The same disjoint word count gives positive all-base mass. The terminals lie in the enlarged paired chamber of rows \(i>I\). Membership in the original chamber additionally requires \(h_2=\cdots=h_{I-1}=0\) through all recurrent phases. At 35, \(I=2\), so this extra condition is empty.

## Ordered originals

For a particular source gap state, order the digits \(a_1<\cdots<a_{35}\) and write their 34 successive differences as 17 opposite pairs. For \(1\le j\le17\), the sum of adjacent differences in positions \(j\) and \(35-j\) is \(g_j\). The middle pair gives \(q_{17}=g_{17}\); the other pairs give \(q_j-q_{j+1}=g_j\). Each has \(g_j-1\) positive integer splits, and the gaps in all these pairs sum to \(B-g_0\). The smallest digit has \(g_0\) choices, namely 0 through \(g_0-1\), with the largest digit at most \(B-1\). Each resulting set of distinct digits has \(35!\) orderings. The exact number of distinct-digit ordered originals with this gap state is therefore \(35!g_0\prod_{j=1}^{17}(g_j-1)\). An ordered word has a unique gap state, and the root offsets and five free coordinates recover the source parameters uniquely. Different sources have disjoint sets of originals.

Since \(g_j>B/19\) and \(B\ge2D>38\), the count per source is at least \(35!(B/38)^{18}\). There are \(L^{17}\) sources, and \(L\ge B/(2D)\) for \(B\ge2D\). Divide the lower count by the total \(B^{35}-B\le B^{35}\) to obtain the theorem at the top of this file. The choice of integer lifts covers every base in the stated range.

## Literal controls at fresh bases

The checker chose bases \(2D+1601\) and \(3D+71393\) and fixed distinct choices of route counts and free integer gaps. For each of those lifts the [literal implementation](literal.cpp) reconstructs a sorted word of distinct digits and computes each new word by sorting and subtracting with the column borrow. It recovers all 18 gaps from each actual word, checks the true adjacent-level comparisons and complementary rows, and compares all gaps and all five labelled values against the reduced law. To determine entry it compares each word with the one 66 steps earlier, since the all-scale theorem proves that the period divides 66. The first equality identifies the first cyclic word 66 steps earlier. The first positive return within the recorded cycle gives the least period, and the gap entry is determined in the same way. The \(d_2\) and radius checks cover the entire cycle.

| Base | \(L\) | Least word entry | Least gap entry | Least period | Transient phases with \(d_2<-1\) | Terminal minimum \(d_2\) |
|---:|---:|---:|---:|---:|---:|---:|
| 232892393201 | 2 | 4826073 | 4826072 | 33 | 1020029 | −1 |
| 349338658793 | 3 | 6206889 | 6206888 | 66 | 1312578 | −1 |

The minimum terminal radius is \(21146051803/4192063077618\) at the first base and \(31724563543/6288095858274\) at the second. Both exceed \(1/396\) by integer cross-multiplication. The first lift assigns equal values to the 7 and 14 labels, so that orbit contributes no factor 2 to the value period. Full original words, margins and labelled charges are in [literal.log](literal.log).

## P10's conjecture and the scope of the refutation

> **Conjecture 10 (Strictness and survival).** For every odd \(n\ge7\), strictness is equivalent to the existence of \(\delta>0\) with \(\liminf_{B\to\infty}\mathbb P_B(Z_B>\delta)>0\), the limit inferior taken over all bases. At a non-strict length, \(\mathbb P_B(Z_B>\delta)\to0\) for every fixed \(\delta>0\).

[Released P10, §7, p. 12](../../../docs/papers/kaprekar-macroscopic-cycles.pdf)

The theorem disproves the necessity direction of the equivalence and the claimed decay at every positive cut. At the non-strict length 35, the event with \(\delta=1/396\) has a positive uniform lower bound over all large integer bases, even with the additional condition that every terminal phase is in \(C\). It concerns a population of original words with transients inside the enlarged paired chamber, and applies to length 35.

## Reproduction

Use Python 3.10 or later and a C++17 compiler. From this evidence directory, use a private temporary directory for generated files. The literal computation takes the project compute lock.

```sh
project=$(git rev-parse --show-toplevel)
work=$(mktemp -d)
python3 reconstruct.py ../settle35/family-certificate.json --controls "$work/controls.txt" > "$work/symbolic.log"
c++ -std=c++17 -O3 -Wall -Wextra -Werror literal.cpp -o "$work/literal"
python3 - "$work" "$project/.task-monitor/compute.lock" <<'PY'
import fcntl
import subprocess
import sys
from pathlib import Path

work, lockfile = map(Path, sys.argv[1:])
with lockfile.open('r+b') as lock, (work / 'controls.txt').open('rb') as inputs, (work / 'literal.log').open('wb') as output:
    fcntl.flock(lock, fcntl.LOCK_EX)
    subprocess.run([work / 'literal'], stdin=inputs, stdout=output, check=True)
PY
for result in symbolic.log controls.txt literal.log; do
  cmp "$result" "$work/$result" || exit 1
done
```
