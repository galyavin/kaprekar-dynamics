# Original-weight bound for the basin of row 181

**Outcome (b).** Row 181's basin has original weight at most **0.6516518719%** (rounded up). The proved entrance mass is **$m(D)\ge92.9921494990\%$** (rounded down). The measure-zero question for switching tails remains open. [Exact rational totals](verification.log)

Let $K$ be the continuous leading Kaprekar map on eleven-digit differences in decreasing order, and let $g(q)=(1-q_1,q_1-q_2,\ldots,q_4-q_5,q_5)$. Use the original ordered nonconstant input population with denominator $B^{11}-B$. Its leading normalized gaps have density $11!\prod_{i=0}^5g_i$ on the five-dimensional simplex. Let $D$ be finite entry into the open leading chamber $\Pi$. [Population and periodic inventory](../entry_pi/README.md) · [Task 139 basin ledger](../e1b_guarded/README.md)

## A finite-time necessary event

Every phase of the entire relative-interior family has $q_1=7/8$, three gaps equal to $1/8$, and a gap equal to $1/4$. The verifier computes those constant gaps at the rational witnesses and checks that their gap rows annihilate the respective phase tangents. Therefore the limit of the sum of the three smallest gaps is at most $3/8$, and the limit of the largest gap is at least $1/4$. Global gap majorization makes the former nondecreasing and the latter nonincreasing. At every finite time in the basin the same bounds hold. [Global majorization](../../odd_global_majorization.md#reflection-and-alternation-prove-majorization) · [Phase and tangent checker](mass_bound.py)

The first difference $q_1$ is nondecreasing and $q_1(Kq)=\max(q_1,1-q_5)$. At every finite time both endpoint gaps $g_0$ and $g_5$ are consequently at least $1/8$. At time three, every orbit tending to row 181 satisfies the explicit event
$$
 E_3=\left\{q:\ g_0,g_5\ge1/8,\quad
       \sum_{i=1}^3g_{(i)}\le3/8,\quad g_{(6)}\ge1/4\right\},
$$
where $g=g(K^3q)$ and $g_{(i)}$ is its $i$-th smallest gap. This includes all switching trajectories tending to the family. [Anchor identity](../e1b_guarded/README.md#limit-constraints-that-can-close-a-basin-despite-a-repeated-word)

## Exact integration

Intersect the time-three event with the image of $K$. Partition by the unordered index triple of the smallest gaps and the index of the largest. All constraints on each cell are rational linear inequalities. Ties between cells have null original weight because each of the finitely many branches of the fixed iterate has determinant of absolute value one. Pull back each cell three times, clipping by the exact image at each inverse step. The final inverse step integrates the sum of the 32 original branch densities. Inverse paths with identical integration domains retain their multiplicities. Task 130's rational simplicial quadrature integrates the resulting degree-six densities. [Inverse paths and quadrature](../en11b/README.md) · [Bound checker](mass_bound.py)

There are exactly 25 cells after image clipping. Of these, 23 have positive original weight. The replay integrates 989 surviving domains across the cells with 31,648 inverse paths counted with multiplicity. The transcript gives all cell weights and their exact sum $m(E_3)$. [Verification transcript](verification.log)

The centre and Two basins together have original weight at most $5124467/80621568$ by Task 130's second-step anchor bound. By Task 139, outside these basins and row 181, almost every orbit enters $\Pi$. This gives $m(D)\ge 1-5124467/80621568-m(E_3)$. The additive bound allows overlap of the necessary events. [Task 130 anchor bound](../en11b/README.md#where-the-step-four-remainder-can-go) · [Task 139 verified non-entry union](../e1b_guarded/README.md#proved-original-entrance-mass)

Run `kaprekar-investigation/.venv-research/bin/python -B kaprekar-investigation/eleven_access/e1c_row181/reproduce.py` from the project root. The runner first replays the 85 closures, the row 181 obstruction and the previous original-weight bound from Task 139. It then computes the new bound with Task 130's native integer density evaluator. Both calculations use private snapshots and the shared compute lock. The Task 139 replay matches the archived transcript byte for byte.
