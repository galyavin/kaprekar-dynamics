# Full-rank zero-export access at length 35

**Theorem.** For every integer \(B\ge232892391600\), uniformly sampled nonconstant ordered words of 35 base-\(B\) digits have probability at least
\[
\frac{35!}{(232892391600)^{17}\,38^{18}}
\]
of reaching an original-chamber cycle with \(Z_B>1/396\). Every cyclic phase in the constructed family clears the cut. This refutes the necessity and every-cut decay claims of P10's Conjecture 10 (Strictness and survival) at length 35. [P10, Strictness and survival](https://doi.org/10.5281/zenodo.24408723)

The proof uses a finite exact projected-return certificate and [guarded physical transport](TRANSPORT.md). The separate finite source-box probes play no role in the original-weight lower bound.

## Projected branch words

Number the rotating labels of the relative root coordinates from 0 to 11 in endpoint reading order from \((x_{12},x_{11},\ldots,x_1)\). At step \(\ell=t\bmod12\), the two observations in the recurrence are
\[
y_\ell(b)=b_\ell,\qquad z_\ell(b)=b_\ell+b_{\ell+7}-b_{\ell-1}
\]
with label subscripts modulo 12. Give them form indices \(\ell\) and \(12+\ell\), respectively. The formula follows because \(x_5\) then has label \(\ell+7\) and \(x_1\) label \(\ell-1\). In label coordinates a step is translation by the phase's assigned kicks. [Fold-order recurrence](TRANSPORT.md)

For a subspace \(V\) of read forms, evaluate observations in \(V\) at their actual values. Force any other \(y_\ell\) to use the positive branch and any other \(z_\ell\) to give zero export. The resulting projected word rule commutes with translation in any direction of \(\ker V\). The evaluated reads agree after translation and the other branches are prescribed. A projected word agrees with ordinary dynamics wherever the forced reads actually satisfy \(y_\ell\ge1\) or \(z_\ell\ge-1\), respectively.

The certificate gives twelve stages indexed \(j=0,\ldots,11\). Start with \(s_0=0\) in label order and \(V_0=\{0\}\). At stage \(j\) trace a prefix of length \(p_j\) and a word of length \(w_j>0\) under the rule for \(V_j\). Both lengths are multiples of 12. Denote the states just before and just after the word by \(a_j,b_j\). The return displacement is \(u_j=a_j-b_j\in\ker V_j\). The certificate supplies another translation \(v_j\in\ker V_j\) and sets \(s_{j+1}=b_j+v_j\). A form \(f_j\) increases the rank by one, so \(V_{j+1}=V_j+\mathbb Q f_j\). Let \(G_j\) consist of all read forms in \(V_{j+1}\) outside \(V_j\). Two verified properties hold for all \(f\in G_j\).

- \(f(u_j)>0\).
- In the complete word shifted by \(v_j\), every read of \(f\) meets its forced threshold.

The new form indices are \(8,4,10,21,20,16,22,0,13,14,18,7\). The matrix with entries \(f_i(u_j)\) is lower triangular with positive diagonal. Its determinant is \(5779085723146507793742893652480\) in the read-form coordinates, so the \(u_j\) are independent. The verifier uses rational rank and an independent exact implementation of the recurrence. [Stage inputs](family-route.txt) · [Verified stage data](family-certificate.json) · [Verifier](verify_family.py)

## Integer inflation and guards

For integer \(L\ge1\), take \(c_j\) independently over the \(L\) integers in \([K_jL,(K_j+1)L)\). Use
\[
s_0+\sum_j(v_j+c_j u_j)
\]
as the initial relative root in label order. At stage \(j\), take its one prefix and \(c_j+1\) copies of the word. Its start is \(s_j+\sum_{i\ge j}(v_i+c_i u_i)\). In a word copy the difference from the synthetic reference state at the same phase is
\[
v_j+b u_j+\sum_{i>j}(v_i+c_i u_i),\qquad b=c_j,c_j-1,\ldots,0.
\]
For the prefix use \(b=c_j\). All added directions vanish on \(V_j\), so the projected rule gives the reference kicks and exports. Each completed word reduces \(b\) by one. The last copy arrives at \(s_{j+1}+\sum_{i>j}(v_i+c_i u_i)\) with phase zero, as required for the next stage.

Fix \(f\in G_i\). After stage \(i\) it is read normally, and all additions from later stages vanish on it. In its own word, a read differs from the reference by \(f(v_i)+b f(u_i)\). The verified minima ensure the forced inequality at \(b=0\), and \(f(u_i)>0\) extends it to all copies. In each earlier stage \(j<i\), and in its own prefix with \(j=i\), its read is at least
\[
M_{j,f}+\sum_{t=j}^i f(v_t)
+ L\sum_{t=j}^{i-1}(K_t+1)\min\{f(u_t),0\}
+ K_i L f(u_i).
\]
Here \(M_{j,f}\) is the minimum reference read over the earlier prefix and word, or just the prefix at \(j=i\). A segment that never reads \(f\) imposes no inequality. The bound uses \((K_t+1)L\) on the remaining counts for earlier copies. Choose each \(K_i\ge0\) so the last two terms together are at least \(\max\{0,H_f-M_{j,f}-\sum_{t=j}^i f(v_t)\}L\) for every such pair \(j,f\), with \(H_f=1\) for primary forms and \(-1\) for secondary forms. The positive slope ensures this choice. The verified integer floors are
\[
(0,0,0,0,0,15,0,1,1,11,1,1).
\]
These inequalities account for every forced branch. Since the other branches are evaluated exactly, the entire itinerary follows \(T_{35}\). Independence of the twelve \(u_j\) gives exactly \(L^{12}\) distinct starts.

After the last stage every form is evaluated normally. The release root in fold order is \((-1,0,2,2,0,11,-1,2,7,2,0,-1)\). It first reaches a least-period-11 relative cycle after 30428 more steps. Its export to entry is 26526, and all eleven recurrent phases have \(h_2=0\). [Exact terminal phase table](family-certificate.json)

## Budgets and physical source

Let \(N\) be the largest absolute reference coordinate over the prefixes, one copy of every word and the common terminal trajectory through a period. Every absolute relative coordinate in the inflated path is bounded by
\[
NL+\sum_j\bigl(\|v_j\|_\infty+(K_j+1)\|u_j\|_\infty\bigr)L.
\]
The spread, including the anchor, is bounded by twice that value. Export sums the fixed prefixes, at most \((K_j+2)L\) copies of each word, and the shared terminal prefix. The clock uses the same count with step lengths. The resulting upper bounds are \(1334114L\) for spread, \(2940551L\) for export, and \(3785731L\) for time through one relative period.

The last stage has at least \(L\) extra copies, each of length 596448 and export 461994. Total export is at least \(461994L\), and a recurrent zero-export phase cannot start before \(596448L\) steps because the last copy still exports. The primary form numbered 7 has slope 77454 on the last return ray and has forced minimum at least 1 in its last copy. In the first copy the corresponding relative read is at least \(1+77454L\), which also lower-bounds lifetime spread.

Take the profile \(3/55\) on roots and \(16/275\) on complements. The margins \(d_3,\ldots,d_9\) are \((1,2,2,1,1,1,1)/275\). The sufficient pairing slack is \(2/1045\), and the homogeneous suffix radius is \(A=1/198\). With \(C_0=2940551\), physical transport gives \(D=116446195800\) and \(L=\lfloor B/D\rfloor\). The anchor and four independent complementary gaps contribute \(L^5\) integer lifts per root at every base \(B\ge2D\). Every complementary row stays positive and the source-to-cycle paths stay full-height and positively paired. Once export ceases, the complementary permutation has order six and the root has least period 11. The physical period is a multiple of 11 and divides 66. All terminal phases lie in \(C\) and have radius at least \(3A/4>A/2\). The \(L^{17}\) distinct interior gap sources give the stated original-word lower bound. [Physical guard and original multiplicity proof](TRANSPORT.md)
