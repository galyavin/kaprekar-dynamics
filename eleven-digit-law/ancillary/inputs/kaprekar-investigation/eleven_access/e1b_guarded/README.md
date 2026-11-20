# Guarded returns at unpaired leading families

**Outcome (b).** Of the 86 open families, 85 have null basins. The only remaining unpaired family is row 181. Its fixed first difference gives a new proved lower bound of **77.3995317809%** on the original entrance mass.

Task 139 concerns the limiting original ordered nonconstant eleven-digit population, with finite-base denominator $B^{11}-B$. The difference law has density $11!\prod_{i=0}^5 g_i$ in five independent gap coordinates. Here $D$ is finite entry into the open leading regular chamber $\Pi$. [Population and catalogue](../entry_pi/README.md) · [Earlier basin ledger](../e1_basins/README.md)

## Exact sorting language

For a row of least period $p$ let $q_j^*$ be its rational witness at phase $j\bmod p$ and $T_j$ its tangent. A phase branch $b=(j,k)$ consists of its unimodular Jacobian $A_b$ and a finite list of integral row guards $G_b$. The branch IDs $k$ are the zero-based order used in `phases_for`. The checker enumerates permutations within the groups of equal $|q_i^*-1/2|$, with both signs at zero. In a nonzero tied group the closed phase itinerary keeps its signs $s_i$, while the strict nearby branch orders the rows $s_i x_i$ for offsets $x=q-q_j^*$. At a zero level, require the last signed row to be positive after ordering the signed offsets. Strict differences between groups and all strict signs at the witness hold throughout a sufficiently small neighborhood of any given point of the relative interior. The branch maps the reference point to its next phase and transports $T_j$ to $T_{j+1}$ independently of the branch. Both statements also hold with the actual limit point as reference.

For a word $b_0\cdots b_{h-1}$ starting at phase $j$, put $P_0=I$, $P_{i+1}=A_{b_i}P_i$. Define its source cone and matrix by
$$
C(b_0\cdots b_{h-1})
=\{x: G_{b_i}P_i x>0\text{ for }0\le i<h\},\qquad P_h.
$$
This is an exact cone representation of the language at every finite length. A word occurs arbitrarily near a reference orbit if and only if its cone is nonempty. Sufficiency follows by scaling a feasible offset so that each of the finitely many iterates remains in its phase neighborhood. A return word starts and ends at the same cyclic phase. For strict homogeneous feasibility use the equivalent rational constraints $G_{b_i}P_i x\ge1$. Sorting-wall hits belong to a countable union of finite preimages of proper affine hyperplanes. We can omit them in a null-basin proof because every branch is unimodular. [Construction and feasibility](guarded.py)

Free concatenation of individually feasible one-period words can introduce infeasible longer words. The exact cone of a product retains both the initial cone and the pullback of the continuation's cone. This distinction also matters for a finite alphabet of period words with a chosen phase.

## Finite normal rows on the language

The finite certificate is a collection of triples $(j,u,C)$, with a primitive integer covector $u$ up to sign, annihilating $T_j$, and a full-dimensional open rational cone $C$. For a seed $u_0$ at phase zero include $(0,u_0,\mathbb R^5)$. For each $(j,u,C)$ and previous branch $(A,G)$, compute
$$
C'=\{x:Gx>0,\ Ax\in C\}.
$$
Whenever $C'$ is nonempty the certificate contains $(j-1,\pm uA,\widetilde C)$ with $C'\subseteq\widetilde C$. The replay checks this for every triple and branch, including cyclic wraparound. To check a cone inclusion, it adds the reversed strict inequality of each containing guard in turn and verifies infeasibility. A zero value of a nonzero row on a nonempty open cone would admit negative values arbitrarily nearby there, so the reversed strict query suffices. Search caches narrower cones under already covered ones for the same phase and row. The stored certificate only needs the displayed inductive closure check.

For any jointly admissible word ending at phase zero, backwards induction puts $\pm u_0P_h$ in the finite cover at the source phase. Along a tail converging by phases to an arbitrary $z$ in the relative interior, the values of $u_0$ at the returns relative to $z$ tend to zero. They belong to the finite set $\{\pm u(q-q^*_0)\}$ computed at the starting section. Each $u$ annihilates the tangent, so the finite set of proper zero hyperplanes is independent of $z$. At least one value is zero. The original basin is null after taking the countable union of finite preimages of those hyperplanes. [Search and separate certificate check](covectors.py)

A strictly feasible repeated word $w$ with matrix $W$ gives a necessary linear constraint on any rational phase-zero seed with a finite orbit. If $w^n$ is feasible for every $n$, the rows $uW^n$ are finite up to sign only if $u(W^{120}-I)=0$, by the degree-five cyclotomic argument of Task 138. We certify repetitions by rational open cones invariant under $W$, or by a strict positive-eigenvalue ray. A ray certificate uses a rational isolating interval for a real characteristic root and rational interval bounds on a signed column of $\operatorname{adj}(xI-W)$ at it. For eigenvalue one the check is a linear feasibility query. On an expanding ray the word can repeat any finite number of times in the local neighborhood after rescaling. On a contracting ray the infinite repeated word stays near and tends to the reference when started sufficiently small. A repeated contracting ray gives a direction of convergence in a Lebesgue null set. An expanding ray only supplies arbitrarily long finite local words. [Repeated-word checks](repeated.py)

## Limit constraints that can close a basin despite a repeated word

**Anchor and smallest-gap sums.** For every leading input, the first difference is nondecreasing and obeys
$$
q_1(Kq)=\max(q_1,1-q_5).
$$
Also each $S_k(g)$, the sum of the $k$ smallest gaps for $1\le k\le5$, is nondecreasing by global gap majorization. Both quantities must be at most their limiting values on a convergent tail. At any phase where the reference gaps have a strict cut between the $k$th and $(k+1)$st smallest, $S_k$ is the sum over that fixed index set on a sufficiently small neighborhood. Its bound and the anchor bound are then homogeneous linear inequalities on offsets from the limit, at that phase. If the bound's row annihilates the tangent, equality is a fixed proper hyperplane for all limit points of the cell, so use a strict inequality after discarding its finite preimages. Otherwise use the weak inequality. Impose these inequalities jointly with all sorting guards, at the start, all intervening phases and the end of a proposed finite word. If there is no feasible word of that length, no generic convergent tail can occur there. For dimension one the verifier splits the entire relative-interior line interval at the rational parameters where two gaps agree, retaining the breakpoint profiles as well. The interval itself uses all strict simplex and signed magnitude comparisons of the periodically transported reference. This checks the whole continuum. [Global majorization](../../odd_global_majorization.md#global-leading-convergence-and-the-periodic-geometry) · [Exact limit constraints](anchor.py)

For an isolated reference, one immediate case is a phase whose two smallest gaps are the tied endpoint gaps $g_0=1-q_1$ and $g_5=q_5$, both separated from every other gap. The anchor at that phase and the next forces both gaps to be at least their limits. The smallest-two-gap sum forces their sum to be at most its limit. Thus both gaps equal their limits at that finite phase. In particular $q_1+q_5=1$ there, a proper fixed hyperplane.

**Finite transverse noncontraction.** Let $V(x)=\max_i |L_i x|$ where the rows span the normal space of the phase-zero family. Test whether a strict word allows $V(x_s)<V(x_0)$ at all the next $k$ returns. On a generic section the signed row $\ell$ realizing $V(x_0)$ is unique after duplicate rows are removed. Fix it and impose $\ell x_0>\pm L_i P_{sp}x_0$ for every tested return and every $i$, together with the initial maximizing inequalities. An empty set of feasible words means every section has some later section within $k$ returns with $V\ge V(x_0)$. Repeating this along a convergent tail produces arbitrarily late sections with $V\ge V(x_0)$, so $V(x_0)=0$. This forces a finite hit of the same affine tangent plane throughout the cell. Ties of maximizing forms and sorting ties again have null finite preimages. [Exact finite-return queries](noncontraction.py)

## Reconstructed results

The Task 138 premise reproduces exactly. Of 161 unpaired cells, 54 archived and 17 additional finite covector orbits close under all independently composed adjacent returns. The upper space is zero on the other 90 rows. The four previously proved null basins (147, 176, 177, 179) account for the difference between those 90 rows and the 86 open rows. [Premise reconstruction](verification.log) · [Task 138 ledger](../e1_basins/README.md#finite-normal-covectors)

| Result | Exact rows (zero-based catalogue IDs) |
|---|---|
| Finite guard-cone cover | 29, 31, 33, 35, 38, 40, 41, 42, 47, 51, 52, 53, 55, 58, 61, 62, 63, 64, 65, 67, 69, 73, 83, 89, 90, 91, 92, 93, 95, 97, 98, 99, 100, 101, 103, 106, 107, 108, 109, 115, 119, 121, 122, 123, 125, 126, 127, 129, 131, 132, 133, 134, 135, 136, 138, 139, 141, 144, 145, 149, 150, 151, 152, 155, 157, 162, 167, 169, 170, 171, 172, 174, 180, 182 |
| Anchor and smallest-gap sums | 114, 142, 159, 163, 178, 183, 185, 188, 189 |
| Transverse noncontraction | 66, 143 |
| Basin open, zero finite-covector upper space on the actual sorting language | 181 |

For all 86 rows, even freely composing the individually feasible guarded one-period words still gives zero covector upper dimension. The guard-cone construction across successive periods closes 74 rows, with 8,130 exact cover triples in total and at most 46 distinct rows up to sign in a phase-zero orbit. The stored fixture records the seed, cone cover, exact counts and repetition exhibits by row. A missing exhibit is a bounded search result. For the nine anchor and gap-sum proofs the exact stopping depths in primitive steps, in table order, are $3,6,4,3,2,4,3,2,1$. The two noncontraction proofs require six and three whole returns. [Row-by-row certificate](certificate.json.gz) · [Replay transcript](verification.log)

### The obstruction at row 181

There are three strictly repeated contracting words at phase zero, with branch IDs $(0,1,0)$, $(1,1,0)$ and $(1,1,1)$. Each has an exact algebraic eigenray satisfying its joint sorting guards; the rays give genuine converging leading tails near every relative-interior reference. The 120th-power conditions from the three matrices intersect the four-dimensional normal row space in zero. Rotating the words gives the same zero-space result at each cyclic phase. The first two words preserve the first difference, and the third preserves the fifth difference at this section. Convergence on a tail repeating one of these words thus forces a finite hit of $q_1=7/8$ or $q_5=1/4$, both proper fixed hyperplanes across the family. A positive-measure basin would require more than these individual repeated itineraries. [Exact cyclic obstruction and isolated repetition exclusions](row181.py)

The anchor and gap-sum guard search also has surviving strict prefixes in all three exact gap-order profiles of row 181. The counts at primitive depths 3, 9 and 15 are exactly $8,70,310$ in each profile. These counts concern finite feasibility. They give no measure of infinite, possibly switching, survivor histories. The surviving contracting words and the finite-prefix growth locate the limit of the current test. [Prefix calculation](row181.py)

## Proved original entrance mass

The original leading orbit converges to a periodic orbit and almost surely enters the invariant leading ball. In that ball, Task 138 and the present null-basin certificates leave $\Pi$, the centre, Two and row 181. A limit in the open chamber $\Pi$ gives entry in finite time. The centre and Two have $q_1\le5/6$ and row 181 has constant $q_1=7/8$ throughout all phases and its whole affine family. By monotonicity of $q_1$, almost all of $D^c$ is therefore in $\{q_1(K^2q)\le7/8\}$. [Global convergence and entry ball](../entry_pi/README.md#chamber-entry-ball-and-periodic-inventory) · [Task 130 anchor](../en11b/README.md#where-the-step-four-remainder-can-go)

For original folded lows $\ell_1\le\cdots\le\ell_5$, the second-step first difference is $\max(1-\ell_1,\ell_4+\ell_5)$. At a cutoff $U$, the first constraint is equivalent to $g_0,g_5\ge1-U$. In each of six cells defined by the number of original differences above $1/2$, all folded lows are affine. The second constraint is the ten inequalities saying each pair of distinct lows sums to at most $U$. Exact rational integration of the original density over these cells reproduces $5124467/80621568$ at $U=5/6$ and gives $3882730875/17179869184$ at $U=7/8$. This uses the six-cell decomposition of Task 130. Hence
$$
m(D)\ge 1-\frac{3882730875}{17179869184}
       =\frac{13297138309}{17179869184}
       \ge77.3995317809\%.
$$
The right-hand decimal is rounded down. If row 181 has a null basin, the conditional lower bound is still $75497101/80621568\ge93.6438\%$. No new positive basin is proved. [Exact original-mass integration](mass.py) · [Replay transcript](verification.log)

Run `kaprekar-investigation/.venv-research/bin/python -B kaprekar-investigation/eleven_access/e1b_guarded/reproduce.py` from the project root. The replay uses Yices exact rational linear arithmetic and SymPy rational root isolation. It runs in a private snapshot while holding the shared compute lock, reconstructs the premise, replays every stored cone cover and witness, and recomputes both original anchor integrals. The cover replay verifies closure directly from its finite states. It does not depend on the search caps.
