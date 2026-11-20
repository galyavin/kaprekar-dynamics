# Eleven-digit leading destination basins and entrance mass

**Task 138, outcome (b).** The proved original leading entrance mass is at least **16.2038159009%**. An exact calculation adds 2.8097941320 percentage points on boxes first entering after step four. The reconstruction confirms Astra's 54 finite-covector certificates and gives 17 additional ones. It proves that the full adjacent-return method has no finite normal covector on the other 90 rows. Three of those rows were already null by Lemma M, and the denominator-25 normal argument extends to one more whole cell. This leaves 86 unpaired families unresolved by the present ledger.

The population is the limiting law of ordered nonconstant eleven-digit inputs with denominator $B^{11}-B$. In five independent gap coordinates its density is $11!\prod_i g_i$. Write $\tau$ for first leading entry into the open regular chamber $\Pi$ and $D=\{\tau<\infty\}$. These results use the continuous class-independent source and leading map. [Population and chamber](../entry_pi/README.md#chamber-entry-ball-and-periodic-inventory) · [Task 130 entry ledger](../en11b/README.md)

## Paired families and the singular ray

Let $b=(1/6,\ldots,1/6)$ in gap coordinates and $h(g)=b+(g-b)/8$. Every coordinate of $h(g)$ is at least $7/48>1/7$. A strictly interlevel-paired step is in a closed affine sorting piece whose interlevel comparisons at $b$ are also strict. Its within-level comparisons vanish at $b$, so their signs, including ties, agree along the segment away from $b$. The affine leading step fixes $b$, giving $K(h(g))=h(K(g))$. An orbit converging by phases to the relative interior of a wholly strictly paired cell stays strictly interlevel-paired throughout a sufficiently late tail. Applying $h$ puts that entire tail in $S=\{\min_i g_i>1/7\}$ with a limit of the same paired word and ties. Tasks 098 and 105 prove nullity in $S$ for 27 of these cells. The inverse homothety and each finite leading preimage preserve null sets, giving global basin nullity. The checker verifies the contraction and the periodic itinerary at all 30 paired witnesses. [Task 098 paired boundary](../entry_pi/README.md) · [Task 105 paired destinations](../destinations/README.md)

On the fixed singular ray the raw radius vector is $(2/3+2f,1/3,-2f,-1/3,-2/3+2f)$ and its grouped magnitudes obey
$$
2/3+2f>2/3-2f>1/3>2f>0
\quad (0<f<1/6).
$$
Small radius offsets $(a,b',c,d,h')$ therefore have next offsets $(a,-d-h',\max(a+b',-c-d),\min(a+b',-c-d),-b'-c)$. In normals $(v,c,d,z)=(a+b',c,d,h'-a)$ the map is
$$
F(v,c,d,z)=(-d-z,\max(v,-c-d),\min(v,-c-d),-v-c).
$$
The independent exact linear query rules out $c\ge d$, input sup norm one and output sup norm at most two after ten steps. The one-step version is satisfiable as a control. By homogeneity a nonzero normal offset cannot tend to zero along a tail that stays in the phase neighborhood. The ray line is the common zero set of the normals for every reference point, so a converging tail must reach that line in finite time. This extends the null basin to the full open range of $f$. [Lemma M ray](../ball_entry/README.md#2-singular-ray) · [Independent query](ray_check.py)

## Finite normal covectors

For all 191 rational witnesses the checker regenerates phases by sorting the eleven entries $1,q_i,1-q_i$. It allows all tied magnitude orders and both signs at every zero, builds the Jacobians from the differences of sorted entries, and identifies the tangent by pulling the phase ties back to the initial phase and imposing the return equation. The tangent dimension agrees with the catalogue. A branch is the product of an adjacent-sum triangular matrix and a signed permutation, with absolute determinant one. Thus every finite preimage of a Lebesgue null set is null.

Let a nonzero normal covector have a finite phase-labelled pullback orbit up to sign. On a tail converging by phases to $z$ in the relative interior, the strict comparisons at $z$ persist. Each allowed affine branch sends the reference phase to the next reference phase. The covector's values at successive phase-zero sections are among finitely many $u(x-z)=u(x-q_*)$, where $x$ is the starting tail section, $q_*$ is the fixed witness, and each returned $u$ annihilates the tangent. These values tend to zero, placing $x$ on one of finitely many fixed proper hyperplanes independent of $z$. Countable union over prefix times and cyclic phases gives a null basin for the whole cell.

The independent search with the original caps returns the same 54 row IDs and counts $13,6,18,15,2$ in dimensions zero through four. All 54 archived seeds reproduce their phase-labelled orbit sizes and tangent conditions. [Independent implementation](catalogue_checks.py) · [Astra's certificate](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md#2-two-common-exclusions-for-the-leading-basins)

### Obstruction to enlarging the search

Choose a reference Jacobian at each phase, compose them into $R$, and let $P_i$ be the reference prefix up to phase $i$. Replacing the reference $A_i$ by $A_i'$ generates a deviation $U_i=P_i^{-1}A_i^{-1}A_i'P_i$. The complete set of one-period matrices generates the same group as $R$ and all $U_i$. A one-phase change gives $RU_i$, and a full period with changes is $R$ times their product in reverse phase order. A finite pullback orbit, even if recorded only up to sign, is finite under the generators and inverses.

For any integer word matrix $W$, a rational row $u$ with finite orbit under $W$ satisfies $u(W^{120}-I)=0$. Its cyclic span has dimension at most five and the minimal polynomial on that span divides $X^n-1$ for some $n$. Its cyclotomic factors have orders among $1,2,3,4,5,6,8,10,12$ whose least common multiple is 120. The checker imposes this condition on the generators and tested pairwise products, includes the tangent constraints, and takes the smallest column space containing these constraints and invariant under the generators. Every finite-orbit normal row must annihilate that space.

For 90 of the 107 rows the upper space has dimension zero. On each of the other 17 it has dimension one and its primitive generator has a closed phase-labelled orbit. Larger seed, coefficient and orbit-size bounds cannot certify any of the 90 using all adjacent returns. A finite orbit of nonzero positive semidefinite rational quadratic forms on the normal quotient would also force a finite normal covector. Its sum would give an invariant positive semidefinite form whose rational kernel contains the tangent. The induced action on the resulting quotient lattice would preserve a positive definite form, hence be finite. A nonzero covector on this quotient would then give the forbidden finite orbit.

The group permits independent choices of adjacent branches at all phases. Actual converging tails must satisfy the sorting guards together. The resulting obstruction concerns the method. A common argument based on the realizable return language or on the mass of survivors subject to those guards remains possible.

### Earlier null-basin arguments

The face words of rows 176, 177 and 179 match the three triangle cases proved in Lemma M §3. Its homogeneous normal map and slicing hold throughout the relative interiors, including the portions below $e=1/30$. These three cells were already null. [Triangle proof](../ball_entry/README.md#3-the-three-triangle-cells)

A cyclic shift of row 147 has the period-three face word of the denominator-25 patch. Its initial section lies in the affine two-plane $u=q_2+q_5-1=0$, $v=q_1-q_3-q_5=0$, $X=2q_1+q_5-2=0$. The checker verifies the tangent and the affine offsets. The strict comparisons and ties match throughout the relative interior, so the homogeneous normal map of Task 026 applies at each point, with a neighborhood depending on the reference. With $w=u_+$ and $k=v+w$, it is
$$
(u',v',X')=(X-w-k,\ |u|+2k_+,\ X+2(-k)_+-w).
$$
After one return $v\ge0$, and thereafter $v'\ge2v$. If $v_3=0$, tracing the formula back gives $v_2=v_1=u_2=u_1=u_0=0$. Then $v_0\le0,\ X_0=v_0,\ X_1=-v_0$ and $u_2=X_1=0$. Thus all initial normals vanish. Convergence forces a finite hit of the same affine two-plane independent of the limit point. Finite preimages give whole-cell basin nullity. [Task 026 patch and normal map](../README.md#3-normal-returns-at-the-denominator-25-family)

## Direct original-weight bound

Use $t_1<\cdots<t_5$ for the initial ordered folded lows $\min(q_i,1-q_i)$. The time-one state is $q'_1=1-t_1$ and $q'_j=1-t_{j-1}-t_j$ for $j=2,\ldots,5$. Each interior folded-low box has 32 disjoint original sheets up to null walls, each obtained by choosing $t_i$ or $1-t_i$ and ordering the resulting coordinates. Their Jacobians have absolute value one. [Initial fold](../../../docs/drafts/research-direction/astra-strategy-12/deductions.md#2-exact-original-weight-at-the-first-two-regular-stages)

The verifier divides $[0,1/2]$ into $M=128$ equal bins and uses boxes with five distinct bins in increasing order. It tracks all five coordinates on a box with shared initial offset variables in $[-1,1]^5$, scaled by $4M$. Each midpoint sign and adjacent comparison of ordered lows is verified against the sum of absolute coefficients of its affine form. A guard verified on the closed box is strict in its interior because the corresponding row is nonzero, by branch nonsingularity. The calculation omits boxes that need more than one branch word.

On strict states the positive interleaving order is the low list $(1-q_1,q_5,1-q_2,q_4,1-q_3)$. It is equivalent to positive within-level signs with strict pairing. Its action on gaps is the regular permutation of order five, so by definition five successive positive interleaving steps are equivalent to $\Pi$ at the start of the run. The program reads from time one and retains the first such run. Time-zero entry would imply time-one membership, so boxes recorded as first entering at time five or later are disjoint from the Task 130 population through time four.

On each original sheet the density is proportional to the product of six gaps between the sorted original differences and endpoints. Put $X_i=24M t_i$. Each $X_i$ has integer mean and variance 12 on its bin. On each side traverse the variables assigned to that side in increasing order, starting from zero. If $d$ is the expected product of consecutive differences so far and $e$ is the expectation of that product times the last coordinate, adding the independent $X$ replaces the pair by
$$
(d,e)\longmapsto(\mathbf E[X]d-e,\ \mathbf E[X^2]d-\mathbf E[X]e).
$$
The chains on the two sides use disjoint variables. Their terminal coordinates give the middle gap after subtraction from $24M$. The expected scaled product is $24M d_0d_1-e_0d_1-d_0e_1$. The native checker verifies this integral against exact tensor Simpson quadrature on three control boxes. The contribution to the original weight from the accumulated sum is that sum times $11!/((2M)^5(24M)^6)$.

The 4,879,817 boxes first entering at steps 5–19 carry exact original limiting mass
$$
w_{\rm box}=\frac{16984163369908523775525}{604462909807314587353088}.
$$
With Task 130's exact mass $w_4$ through step four this proves $m(D)\ge w_4+w_{\rm box}$, at least **16.2038159009%** rounded down. It also proves $m(D^c)\le83.7961840991\%$, which bounds the remaining unpaired union. Subtracting the positive Two lower bound yields a slightly stronger bound on that union. [Exact box transcript](mass.out) · [Integer verifier](folded_mass.cpp) · [Rational bound replay](reproduce.py)

If the remaining unpaired union is null, the conditional bound is still $m(D)\ge 1-5124467/80621568$, or at least 93.6438% rounded down. The centre's known joint contribution is $\int_{\{K^jx\to c\}} f(x,0,c)\,d\mu(x)$, where $c=(5,4,3,2,1)/6$. Its mass is open. The outstanding obligation for Two is transfer from its original basin into the bounded transverse physical landing where the $4a/9$ clock is proved. [Centre contribution](../../../docs/drafts/research-direction/astra-strategy-12/deductions.md#1-the-balanced-basin-has-a-known-clock-whatever-its-weight) · [Two basin](../entry_pi/README.md) · [Task 130 centre-or-Two bound](../en11b/README.md#where-the-step-four-remainder-can-go)

**Budget stop.** The return obstruction sits in independently composed adjacent branch matrices. The direct bound is limited by guard-straddling boxes and $M^5$ grid growth. At $M=128$ the exact scan has 252,399,461 ambiguous boxes out of 264,566,400 tested. The ambiguity weight requires integration over those boxes.

Run `python3 -B reproduce.py --recompute-mass` from this folder with Yices and Clang on the path. It replays the catalogue, expansion and mass checks in a private snapshot and takes the shared compute lock for the native scan. Omitting the flag replays the rational calculations from the retained mass transcript.
