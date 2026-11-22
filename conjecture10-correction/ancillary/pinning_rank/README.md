# Pinning rank and chamber capture

Task 247 works at a fixed odd length \(n=2k+1\ge7\). Put \(m=k+1\), \(I=\lfloor m/2\rfloor\), \(P(j)=\min(2j,n-2j)\) for \(j>0\), and \(P(0)=0\). Let \(r\) be the order of \(P\) and \(J=\{0\}\cup\operatorname{orb}_P(1)\). The secondary rows are

\[
 d_i(g)=\sum_{j=m-i}^{m-1}g_j-\sum_{j=0}^{i-1}g_j,\qquad 2\le i\le I.
\]

The ambient space \(D\) has zero gap sum and gaps equal on \(J\). Its independent coordinates are the differences \(g_j-g_0\) at slots outside \(J\), so \(\dim D=k-r\). In these coordinates the restriction of \(d_i\) is just its vector of coefficients on slots outside \(J\). P12's \(K\subset D\) requires \(d_i(P^t g)\ge0\) for every secondary \(i\) and \(0\le t<r\). Write \(d=\dim\operatorname{span}K\) and \(c=k-r-d\). See [P12, Gap coordinates and original measure, The capture exponent and the breach bound](../../../docs/drafts/p12-release/p12.tex).

## Codimension

Let \(U=D\cap\operatorname{Fix}(P)\) and \(\ell_i=d_i|_U\). Call \(i\) implicit if \(\ell_i\) vanishes on \(W=K\cap U\). For any \(u\in K\), its orbit average belongs to \(W\). Every term \(d_i(P^t u)\) is nonnegative, so an implicit row vanishes at every phase of every point in \(K\). If \(i\) is absent from the implicit set, there is a point of \(W\) on which \(\ell_i>0\). Summing one such witness per nonimplicit row gives an invariant point strictly positive on all nonimplicit rows at every phase. In the common kernel in \(D\) of the rotated implicit rows, a relative open ball about that point remains in \(K\). Thus this kernel is \(\operatorname{span}K\), and \(c\) is the rational rank of the **individual** rotated implicit rows.

A nonnegative relation \(\sum_i\lambda_i\ell_i=0\) certifies that its support is implicit. A relation with positive weight on every implicit row exists. For each implicit \(i\), the form \(-\ell_i\) belongs to the dual of \(W\) in \(U\) and is therefore a nonnegative combination of the \(\ell_j\). Add \(\ell_i\) and sum the resulting relations. Together with an invariant point positive off the support, such a relation certifies the entire implicit set. A strictly positive vector in the kernel of the implicit \(\ell_i\) also shows that the nonnegative relations span their entire coefficient kernel, by adding a sufficiently large multiple of that vector to any element of the kernel.

A row is **literal root-only** if its prefix and suffix slots all belong to \(J\). It then restricts to zero on \(D\). Any row with a complementary slot restricts to a nonzero form, since the two blocks are disjoint and the complement offsets in \(D\) are independent. Consequently \(c=0\) exactly when all implicit rows are literal root-only, equivalently when some invariant profile is strictly positive on every row with a complementary slot. This is complement-strictness. The case \(D=\{0\}\) is vacuous and is excluded from the proper-length count.

The seven requested ranks are certified as follows.

| \(n\) | \(r\) | \(k-r\) | \(d\) | \(c\) | implicit row indices |
|---:|---:|---:|---:|---:|:---|
| 25 | 10 | 2 | 1 | 1 | 2, 6 |
| 27 | 9 | 4 | 2 | 2 | 4 |
| 35 | 12 | 5 | 5 | 0 | 2 |
| 39 | 12 | 7 | 2 | 5 | 4 |
| 97 | 24 | 24 | 24 | 0 | 2 |
| 119 | 24 | 35 | 12 | 23 | 8 |
| 277 | 46 | 92 | 1 | 91 | 11, 12, 13, 14, 15, 19, 21, 26 |

At 277, even the rotations of all linear combinations from the full coefficient kernel of the invariant rows on the implicit set have rank 90 over \(\mathbb Q\). The individual rotated rows have rank 91. The loss persists after taking the span of all nonnegative Gordan relations.

## Original-word charge in the ordinary chamber

Here \(A_B(T)\) is P12's gap-time event. The original starts in \(C\), remains in \(C\), and is cyclic by integer time \(T\). In particular, the cycle is contained in \(C\). Assume \(B\ge2\). In \(C\), write \(z=g_k-g_0\), \(e=\operatorname{sgn}(z)\) and \(h=(-1-z)_+\). The kick in \(g^+=Pg+b(z)\) is

\[
 b(z)=e(1,-2,1,0,\ldots)+h(-1,1,1,-1,0,\ldots).
\]

The relative-root vector \(x=(g_j-g_0)_{j\in J\setminus\{0\}}\) follows the permutation on nonzero root slots plus the root coordinates of \(b(z)-b(z)_0\mathbf 1\). Both depend only on \(x\), including when slot 3 is in \(J\). The outgoing slot 1 receives \(z-3e+2h\ge-2\) in relative coordinates. At every other nonzero slot the relative kick is at least \(-1\). The label at slot \(k\) passes through slot 1 in the preceding \(r-1\) steps, so after \(r\) chamber updates \(z\ge-2-(r-1)\). All kicks from then on are \(O_n(1)\). This bound is independent of the primary losses in the first \(r\) steps.

By P12's cyclic tube lemma, every target phase is within \(O_n(1)\) of \(B\mathbf 1/m+K\), and in particular has relative-root norm \(O_n(1)\). The number of physical integer states in this tube is \(O_n(B^d)\). To see the lattice count, bound the physical region by a box of side \(O_n(B)\), include \(K\) in its \(d\)-dimensional span and thicken the normal directions by a fixed distance. Enlarge the resulting region by a bounded fundamental cell of the sum-\(B\) lattice. The volume is \(O_n(B^d)\), which bounds the number of disjoint cells.

Take a common target time \(t_*=\max(r,T)\). All eligible paths are cyclic at \(t_*\). From time \(r\) to \(t_*\) the bounded kicks and reverse permutation show that \(x(r)\) has norm \(O_n(T+1)\), hence \(O_n((T+1)^r)\) possible values. There are also \(O_n(B^r)\) possible physical relative roots. For two eligible states at time \(r\) with the same relative root, all subsequent kicks up to \(t_*\) agree. Their difference at \(t_*\) is \(P^{t_*-r}\) applied to their original difference, so the map to the target states is injective on the fibre even if the two orbits first become cyclic at different times. Each relative root therefore accounts for at most \(O_n(B^d)\) time-\(r\) states. P12's all-height predecessor bound

\[
 M_n=1+\sum_{\theta=1}^k3^\theta(2\theta)!
\]

costs at most \(M_n^r\) for the first \(r\) steps. Its original ordered-word probability per gap state is \(O_n(B^{-k})\). Thus

\[
 \mathbb P_{n,B}(A_B(T))
 \ll_n B^{-c}\left(\frac{\min(T+1,B)}{B}\right)^r.
\]

In particular this gives the requested bound for \(0\le T\le B\), and \(O_n(B^{-c})\) for uninterrupted residence with unlimited entry time. For the latter take a common time beyond all entry times at a fixed base, which exists because the physical state set is finite, and use the \(O_n(B^r)\) bound. A prescribed head of \(H\) arbitrary physical steps costs at most \(M_n^H\) when the rest is an uninterrupted ordinary chamber path with the clock measured from \(H\). A radius mark can only restrict the event. P12's separate static seeded-class order is \(B^{-(r+c)}\). Its marked two-sided dynamic estimate applies at the 47 certified strict lengths in its stated range and for its stated cuts and clocks. The fibre argument proves an upper estimate.

### Changed guards

The counting proof extends to a guarded path family with a common deterministic time after its reset, \(O_n((T+1)^r)\) possible relative roots there, and \(O_n(B^d)\) possible common-time targets, if the source-to-target map in each fibre has bounded multiplicity independent of \(B,T\). A sufficient condition for multiplicity one is that the suffix is \(g_{u+1}=Pg_u+b_u\) with the accumulated offset \(\sum_uP^{t_*-1-u}b_u\) fixed by the source relative root on that fibre. A uniformly bounded number of offset values gives the same bound with the corresponding factor. Pulling back a bounded head also requires its predecessor factor. For a family with pairing breaches or changed terminal guards, the reset, offsets and target count each require verification for that family.

For positively paired states the kick coefficients are explicit. The relative-root evolution can change when secondary breaches occur. Along a paired segment starting at \(H\) a label originally in slot \(j\) has net displacement

\[
 L_j=\sum_{u=H}^{t_*-1}
       (b_u)_{P^{u-H+1}j},\qquad
 (g_{t_*})_{P^{t_*-H}j}=(g_H)_j+L_j .
\]

At 27 the implicit row \(d_4\) restricts to \(e_{12}-e_3\) on the complement. Its three rotations on slots \(3,6,12\) give two independent differences. On each chamber-cycle phase these differences are \(O(1)\) by the cyclic tube. In the paired formula the kick values at destinations \(3,6,12\) are \(h_2-h_1,h_3-h_4,h_6-h_7\) respectively. The Task 241/243 hypothesis bounds oscillation of their net displacements on the three moving labels by \(V\), pulling back two \(O(1+V)\) width constraints. These are the two pins in the ordinary chamber charge. The registered head and oscillation estimate and its residual events retain their stated hypotheses. See [Task 241](../../../tasks/completed/241-forget-probe-27.md) and [Task 243](../../../tasks/completed/243-recon-l-forget27.md).

The unrestricted enlarged-chamber upper law already fails at \(n=23,T=0\). In the chamber that keeps positive pairing and drops literal root-only secondary guards, all secondary guards are dropped at this length. P9 supplies positively paired fixed gaps at base \(71a\), odd \(a\ge7\), with \(12q_1-11(71a)=-a\). At a common base \(B=71A\) for large odd \(A\), choose \(A/2\le a\le3A/4\) with \(a\equiv A\pmod{12}\) and translate all 12 gaps by \(u=71(A-a)/12\). This increases the sum by \(12u\). Every row \(d_i\) is fixed and each adjacent pairing margin gains \(u\), since each pair in level \(i\) gains \((12-i)u\). The paired kick is unchanged and \(P\mathbf1=\mathbf1\), so the new gaps are fixed. Their centred first difference numerators \(-a\) distinguish \(\Theta(B)\) targets and give \(Z_B\ge1/1704\). Their normalized gaps lie in a compact interior subset of the simplex, giving original weight \(\asymp B^{-11}\) each by P12. The resulting marked mass is \(\gg B^{-10}\), whereas the claimed bound with \(r=11,c=0,T=0\) would be \(O(B^{-11})\). This construction is from the fixed family in [P9, Two excluded lengths with macroscopic fixed points](../../../docs/drafts/p9-release/p9.tex).

## Prime criterion and exact finite count

Under \(P\), the prefix slots in \(d_i\) map to \(0,2,\ldots,2i-2\) and the suffix slots to \(1,3,\ldots,2i-1\). Membership in \(J\) is preserved. Delete exactly those \(i\) for which \(1,\ldots,2i-1\) all belong to \(J\). In particular \(d_2\) is literal root-only exactly when \(3\in J\), and if \(3\notin J\) every secondary row is retained. At a proper prime \(p\) this is the test \(1,\ldots,2i-1\subset H_p=\langle -1,2\rangle\). At least one row is retained, since the final prefix reaches all positive folded slots except possibly \(k\in J\).

Write \(G_p=\mathbb F_p^\times/H_p\). P11 identifies the invariant row with the walk vector \(V(2i-1)=\sum_{b=1}^{2i-1}(-1)^{b+1}\mathbf e_{[b]}\), with \(\mathbf e_{H_p}=0\). By the codimension criterion, a proper prime is complement-strict exactly when the **retained** walk vectors lie in a common strict open half-space through zero, equivalently their convex hull excludes zero. A retained zero vector is an obstruction. At folded index two the endpoint is retained and positive in the nonsquare offset orientation, since \(V(2I-1)=h(-4p)/2>0\) by P11. Thus every retained \(F(i)\) must be \(<1\), where \(V(2i-1)=(1-F(i))/2\). This threshold recovers ordinary strictness if \(3\notin J\). With \(3\in J\) it describes the additional non-strict prime class. See [P11, Prime walk and Index-two identities](../../../docs/drafts/p11-release/p11.tex).

The exact census tests every odd \(7\le n<10^4\) with \(k-r>0\) and \(3\in J\). These are precisely the candidates for proper complement-strictness without strictness. A positive certificate is an integer invariant orbit-offset vector positive on every retained row. A negative certificate is a nonzero, nonnegative integer relation on retained invariant rows. The result is 274 lengths, of which 94 are prime and 180 composite. The additional class includes 26 index-two primes in this range. For example, 97 deletes \(d_2\) and all its retained \(F(i)\) are \(<1\). At 17 and 41 there are no deleted secondary rows and all secondary \(F(i)\) are \(<1\). They are P11's ordinary strict index-two baseline in its stated finite range. The census also certifies 3343 and 9931, which account for the two additions to the preliminary count of 272. The infinitude of additional complement-strict primes remains open.

## Proper facets at positive codimension through 101

For each odd \(7\le n\le101\), an invariant interior point and relation certify the implicit indices and the rank. For the eleven lengths with \(c\ge1\), the replay lists every secondary row that vanishes on all of \(K\), defines a proper facet, or defines a proper face of higher codimension. Row indices are shown at phase zero, since \(P\) acts on \(K\) and preserves all face dimensions. Two rows can define the same facet. At 27, \(d_4\) vanishes throughout \(K\), \(d_2,d_3,d_7\) are proper facets, and \(d_5,d_6\) cut out only the apex.

For the certificates, restrict to \(S=\operatorname{span}K\) in the exact kernel coordinates. Orient each nonzero normal so it is nonnegative on \(K\), and divide out positive scalar factors when deduplicating. Given the normal \(a\) of interest, give a point in \(K\cap\ker a\). For the other distinct normals that vanish at that point, supply a positive combination equal to a positive multiple of \(a\). This forces them to vanish on the whole face. Every remaining normal is strictly positive at the point. By taking a small open ball in the common kernel of the active normals, the span of the face is exactly that common kernel. Codimension one in \(S\) gives a facet, and dimension zero gives the apex. On a ray the apex is a facet and is listed in the facet column. An implicit row defines the improper face \(K\).

The universal necessary-facet claim is false even for positive limiting weight. At length 9, \(r=3\), \(J=\{0,1,2,4\}\), and the only secondary comparison \(d_2\) is positive on the invariant offset that raises slot 3, after centring. P12's strict source transport uses a near-uniform profile on that ray with a positive secondary margin. Fix the source width and take \(T=\lfloor tB\rfloor\) for a sufficiently small fixed \(t>0\) in that proof. The transport and lifetime bounds keep the resulting family, whose original mass has a positive base liminf, at \(d_2\ge\eta B\) for some fixed \(\eta>0\) through capture and recurrence, with a positive radius cut. The distance from the facet \(\{0\}\) in centred normalized coordinates is bounded below by \(\eta/\|d_2\|_2\) at every phase. At 35, \(d_2\) is literal root-only and vanishes on \(D\), so that physical guard defines the improper face \(K\). Earlier wall visits at 35 are a separate dynamical question. See [P12, Marked lower bound](../../../docs/drafts/p12-release/p12.tex).

## Replay

Run `python3 -B verify.py` from this folder. It regenerates the rows and checks integer or rational identities and inequalities, then recomputes the ranks over \(\mathbb Q\). `short_certificates.json` has one record per tested length with its invariant integer interior witness and positive relation, indexed by increasing implicit row index. For positive-codimension lengths through 101 it also has the face witnesses in the reduced coordinates. `census_certificates.jsonl` records one positive invariant witness or one negative relation for every eligible census candidate. The discovery script uses SciPy to propose witnesses and exact arithmetic to recover or validate them. The verifier uses only the Python standard library and no feasibility status from SciPy.
