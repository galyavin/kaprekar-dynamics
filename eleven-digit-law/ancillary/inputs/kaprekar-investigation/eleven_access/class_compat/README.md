# Base periods after projection of eleven-digit entrant routes

The measure is the original nonconstant word law on the leading entrance event $D$, with original denominator $B^{11}-B$. The retained marks are the original normalized differences, the first cyclic-word clock divided by $B$, and the normalized differences of that word. The last mark tends to $(5,4,3,2,1)/6$ by Theorem B′. [Truncation argument](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md#1-the-restricted-law-and-its-arithmetic-parameter)

**Result.** All fine-base projected joint laws agree within their class modulo 12 on $D$. The minimal common joint-law modulus is exactly 12. The proof covers the full entrance mass, which Task 141 bounds below by $92.9921494990\%$ on the original denominator.

## Premises

The marked-box check repeats the independent determinant pullback of Task 136. Re-execution gives identical fractions at all twelve eligible residues modulo 72. Each of the two labels has conditional box mass $1/28512$. The conditional profile of L has two weights $2/3$ and $1/3$.

| $B\bmod24$ | $r_2=603$ | $r_2=639$ | Pooled profile |
|---:|---|---|---|
| 5 | M2 | L | $(\mathrm{M2}+\mathrm{L})/2$ |
| 17 | L | M2 | $(\mathrm{M2}+\mathrm{L})/2$ |
| 11 | L | S | $(\mathrm{L}+\mathrm{S})/2$ |
| 23 | S | L | $(\mathrm{L}+\mathrm{S})/2$ |

The two labels both fill the same leading box with the same spatial density. Pooling therefore preserves the joint law with original position, clock and terminal mark. The raw box reader enumerates $5,11,17,23,29,35,41,47,53,59,65,71$ modulo 72. The comparison proves that the labeled profiles repeat modulo 24. The full seven-graph box projection was already proved to repeat modulo 12, including all cohorts and the alignment collars. [Task 136](../r2_routes/README.md#mass-accounting-and-period) · [Box lift and seven-graph table](../box_structure/README.md)

For a line sheet $X=p+kT+nL$, $B=b+kP$, write $L=\ell u$ with $u$ integral and primitive. Choose an integral normal projection $N$ with kernel $\mathbb R u$, and set $a=NT/P$. The shift $6m$ is compatible with the same physical normal equation exactly when $mT/P$ differs from an integer vector by a scalar multiple of $u$. Equivalently $m a\in\mathbb Z^5$, since $N\mathbb Z^5=\ker(c)\cap\mathbb Z^5$ for $N=I-uc$, $cu=1$. The least shift is $6\operatorname{lcm}_i\operatorname{den}(a_i)$. The direct congruence check $mT+c'u\equiv0\pmod P$, $c'\in\{0,\ldots,P-1\}$, agrees exactly with the canonical normal strides $6,6$ on two archived class-6 rank-one sections, and with $12,18,72$ on three constructed nonhomogeneous controls. [Premise checker](premises.py) · [Exact premise ledger](premises.json)

## Entire rank-one continuation uses the centered paired map

Let $K$ be the continuous leading difference map, which sorts the leading borrowed digit values. At every noncentral rank-one arrival or subsequent macroscopic exit, the leading point is periodic for $K$. Indeed, at an exit of a $p$-phase translating word the limiting endpoint is fixed by the same leading composite by continuity, including when the last guard is the first wall reached. A bounded short update applies $K$ to the leading position, so all shorts before the next macroscopic word run around this same finite leading orbit. A positive rank-one leading point cannot reach the centre solely by these updates because its orbit is periodic.

Task 142 proves coverage by a graph that includes all admitted inherited rank-one continuations. Its positive-duration rank-one words all have paired phases, and every noncentral starting line is homogeneous. The leading point after the intervening shorts belongs to a periodic orbit whose every phase is at least weakly paired: a positive-duration paired word covers the whole orbit. The phase at the start of the shorts is on that same orbit. On the homogeneous arrival line an adjacent interlevel slack is of the form $1+f(x)$, with $f$ a centered linear form. Its zero set meets the one-dimensional line in at most one point. Leading positivity of the last difference has the same property, with constant $1/6$. Away from those walls the arrival is strictly paired. Its next leading chart is homogeneous because the paired centered map is piecewise linear and homogeneous in the leading limit. The same argument applies inductively to each of the short phases. The excluded parameter values pull back to proper walls of the original source because the rank-one chart has nonzero variation. Countable union over finite records preserves the original null mass. Each repeat is already paired by the certificate. This proves physical pairing on almost every fixed rank-one route through each positive terminal collar for sufficiently large bases, even on the outer sections. Inside the collar use centered permanence, including the cycle. [Task 142 closure and physical guards](../z1b_rank_one/README.md#exact-graph-certification) · [Paired base covariance](../../odd_base_covariance.md)

The active rank-one exit reader is also autonomous. Write its winning guard as $g(X)=vX+bB+c$, with drift $d$ and $s=-vd>0$. The outgoing leading chart has $v x'=-b$ everywhere, by $x'=x+d(vx+b)/s$. The verified outgoing line is homogeneous or the exit is central, so its affine hull contains zero and forces $b=0$. A competitor tied identically at leading order must have the same ratio of base coefficient to decrement, since the incoming line is homogeneous. Thus the tied winner also uses a centered comparison. Other leading walls are proper and have zero original mass. This proves that exact rank-one switches, winning guards and their congruence splits can be read on the restored centered coset at every finite depth.

The significance for arithmetic is that the rank-one continuation needs only its **incoming centered coset**. For any admissible exact rank-one arrival $a=0$, so $T=\tau u$ with $\tau$ integral. The least shift by a multiple of its inherited stride $P$ that restores its exact centered coset is

$$
P_{\rm cent}=P\,\frac{\ell}{\gcd(\ell,\tau)}.
$$

This follows from $kT\in L\mathbb Z$ precisely when $\ell\mid k\tau$. The duration and intercept of the original prefix remain as marks. On a fixed route, for a shift by $P_{\rm cent}$, the two original source cosets are translates. The difference of the centered outputs of corresponding points belongs to the image lattice $L\mathbb Z$. The original direction lattice maps onto $L\mathbb Z$ by definition, so a fixed extra original lattice displacement cancels it. The combined translation thus matches the centered incoming state exactly. It changes original positions by $O(B^{-1})$ after normalization and the physical elapsed prefix by $O(1)$, by the exact affine route formulas. The boundaries of each strict chart and the alignment collars have asymptotically zero original weight. The two exact centered continuations then agree because the physical paired map is autonomous. First cyclic-difference time differs from first cyclic-word time by at most one, and the leading terminal mark is common. The bounded displacement preserves limiting weight on each original chart by the same weighted Riemann counting as Task 136.

Equivalently, on a progression restoring the whole incoming coset, choose its original representative so that the centered base-shift vector vanishes. Subsequent autonomous returns then require no new base refinement. The paired suffix argument applies at every collar depth, so it yields equality of whole clock laws under the same shift after the Task 135 truncation. There is no need to require a uniform number of inner stages. [Original route weights](../r2_routes/README.md#the-original-lattice-calculation) · [Inherited tangent order](../s1b_signed/README.md#a-period-for-a-whole-paired-continuation)

## Exact stride inventory

The exact ingress traversal repeats the six Task 137 inherited rank-two closure summaries. It records every stride read at ranks five through two, and independently computes the full centered period before canonicalizing each emitted rank-one arrival. This includes every inherited sheet needed for a given entrance chart. The rank-two domain compression shares route verification only for identical physical sheets and covered leading geometry. No probability is taken from that compression.

| Root $B\bmod6$ | $P=6$, $P_{\rm cent}=6$ | $P=6$, $P_{\rm cent}=12$ | $P=6$, $P_{\rm cent}=66$ | $P=12$, $P_{\rm cent}=12$ |
|---:|---:|---:|---:|---:|
| 0 | 109,770 | 29,689 | 3,248 | 35,426 |
| 1 | 119,196 | 58,281 | 0 | 74,812 |
| 2 | 109,803 | 29,689 | 3,248 | 35,420 |
| 3 | 119,346 | 58,281 | 0 | 75,073 |
| 4 | 109,802 | 29,689 | 3,248 | 35,434 |
| 5 | 119,269 | 58,281 | 0 | 75,010 |

The counts are exact records of closure, not weights. Ranks five, four and three read only stride 6. Rank two reads 6 and 12. Every rank-one arrival has normal stride 6. The only additional centered stride is 66, and it occurs exclusively in the even roots. In particular the joint weak law already exists modulo 12 in all six odd classes, and modulo 132 in all six even classes. For each $r\in\{0,2,4,6,8,10\}$ the potentially distinct classes are exactly the eleven candidate residues $r+12j\pmod{132}$, $j=0,\ldots,10$. Each $r\in\{1,3,5,7,9,11\}$ has one projected class $r\pmod{12}$ at this step. The separate marked box checkpoint uses the two arms $r=5,11$ only and projects as above. [Exact inventory](stride-inventory.json) · [Reproduction](stride_inventory.py)

The factor 11 occurs on a realized original route. In root 0 one uncompressed $5\to4\to3\to2\to2\to1$ route uses an unpaired ten-update third word whose winning guard has decrement 66 and explicit base coefficient $-1$. The third split has index 11. The inherited rank-one output at $B=6k$ is

$$
(6,18,-6,-30,-24)+6k(1,1,0,-1,-1)+66n(1,1,0,-1,-1).
$$

On the fixed strict original polytope of the replay the joint original index is 5,544. The proved direct-source mass of this one route exceeds $4.18\times10^{-13}$ on the original denominator. Its full lattice rotates through all eleven phases under successive twelve-shifts. The weight is for the arithmetic occurrence before reading the final clock. [Uncompressed original replay](period66-witness.out) · [Witness checker](period66_witness.py)

## The projected clock loses the prime 11

The second exact traversal retains every period-66 arrival and closes the canonical graph from their allowed intervals. Every such arrival has rank history `543221`. Each of the three even roots has 3,248 inherited arrivals and 98 canonical starting sections. The same 98 sections close the continuation. The exact winning returns are all terminal.

| Root $B\bmod6$ | Terminal decrement and count | Nonterminal winners |
|---:|---|---:|
| 0 | $66$ on $56$, $96$ on $30$, $192$ on $12$ | 0 |
| 2 | $96$ on $98$ | 0 |
| 4 | $24$ on $70$, $96$ on $28$ | 0 |

The certified short prefix and winning return on a strict canonical cell depend on its exact affine normal offset and tangent line. In a structural comparison the row annihilates the tangent, and the return drift is also constant under tangent translation because the return fixes the tangent. A proper leading wall has null original weight. All eleven centered phases of a period-66 arrival are thus covered by the same cell and terminal word on the common leading input chart. A winning return with word length $p$ and decrement $s$ adds $pD_a(y)/s$ to the scaled clock. Its physical remainder changes only the bounded terminal point and integer rounding, so it cannot change the leading clock or terminal mark. The remaining physical clock and any early-capture discrepancy are $O(1)$ on each fixed terminal record, by the centered and nonzero-drift bounds. The count is an exact overcover of guard choices. Each compared route retains its own weight from the original lattice. [Exact terminal reads](prime-reader.json) · [Verifier](prime_reader.py) · [Centered clock tolerance](../z1b_rank_one/README.md#clock-tolerance)

On a period-66 arrival, temporarily enlarge the **inherited** lattice $L\mathbb Z$ to $(L/11)\mathbb Z$. The base slope belongs to the larger lattice because it has order 11 in the original quotient. Thus the eleven input subcosets inside the enlargement exhaust the fine phases met by twelve-shifts. A centered remainder read $aX+c\pmod s$ gives the same support and weights on each of the eleven input subcosets whenever $\gcd(s,aL)=\gcd(s,aL/11)$. This holds automatically for $11\nmid s$. For each shared remainder, the exact affine exit carries an original subcoset into the enlarged exit coset. Their index is either 11 or one because an affine image of a quotient of order 11 has order dividing 11. Bounded paired switches obey the same index statement. Structural controls depend only on the common leading line and exact normal offset. If the gcd test holds at every read before the terminal or collar outcome, induction matches the finite projected laws with the same weight, leading clock and input. A rank-dropping exit to the centre requires no gcd test because its duration does not read its final remainder, and the residual scaled clock tends to zero. Collar truncation transfers the finite projected equalities to the final clock law. [Exit congruence and source weights](../r2_routes/README.md#the-original-lattice-calculation)

On a canonical overcover let $L_*=6u$ be the primitive physical tangent column at the beginning of a section. A sufficient test independent of the inherited tangent stride is that the largest power of 11 dividing the winning decrement $s$ divides the guard's increment on $L_*$ transported through the short prefix. In fact each inherited or enlarged tangent column through that section is an integer multiple of $L_*$, so multiplying an inherited index by 11 cannot change its gcd with $s$. This verifies a probability equality on the inherited cosets without treating the canonical overcover as a source distribution.

## Compatibility on all twelve classes

Fix an original strict entry prefix and one post-entry history through the first rank-one arrival. Its compatible progression has stride 6 or 12 and its original coset index and leading position and clock functions are fixed. On a class modulo 12, the rank-one branches of centered period 6 or 12 have identical centered cosets across all refined bases and therefore identical projected laws by the autonomous paired suffix argument. On the remaining branches, the period-66 phases have identical projected laws by the direct terminal read. The associated source density and coset weight are the same by Task 136. This compares the original input and the clock together on every such history.

Carry this equality through finite-record collar truncations. Both refinements have full coverage up to vanishing original weight by Tasks 137 and 142. The remaining-clock error is at most $(16/3)\eta+o(1)$ at the verified positive collar of radius $\eta$. Exhaust fixed leading entry prefixes in $D$ and send $\eta$ to zero. The finite measure restricted to $D$ therefore converges along every class modulo 12. The joint distinction between different classes can be localized in original position, as proved below. [Finite-record transfer](../r2_routes/README.md#the-original-lattice-calculation) · [Restricted-law truncation](../../../docs/drafts/research-direction/astra-strategy-13/deductions.md#1-the-restricted-law-and-its-arithmetic-parameter)

## A necessary divisor for the joint law

Any finite modulus governing the full entrant **joint** law must be divisible by 12. First, on the fixed positive-weight repeated-failure box the two arms $5$ and $11$ modulo 12 have distinct projected clock laws, with an exact positive difference of the means. Restricting a joint law to this box is legitimate because the marginal of original position is absolutely continuous and the spatial boundary has measure zero. If $4\nmid M$, then $\gcd(M,12)\mid6$, and there is a single class modulo $M$ meeting both arms in infinite progressions. Their distinct localized limits contradict a common joint law on that class. Second, on the fixed direct-input paired cone of Task 140 the projected means for residues $3$ and $5$ modulo 6 differ exactly. Its position boundary is likewise null, and the paired covariance gives the same cone law for $9$ as for $3$ modulo 12. If $3\nmid M$, then $\gcd(M,12)\mid4$, and a single $M$ class meets both $9$ and $5$ modulo 12. This again contradicts a common joint law. Original position allows both localizations. An integrated marginal can cancel distinctions on different original regions, so these two deductions concern the joint law. [Box signed clock](../r2_routes/README.md#mass-accounting-and-period) · [Exact paired cones](../s1_selection/README.md)

## Comparison with rule L

The translation order on the inherited unprojected whole-line cosets supplies the candidate common modulus 132. Retaining the auxiliary box label shows an additional mod-24 alternation, verified on the mod-72 checkpoint refinements. Both effects disappear under the required projection. The partner box labels exchange weights, and each phase of the rank-one index-11 quotient has the same projected terminal clock. Equivalently the nonconstant Fourier components within these coarse classes vanish for the projected measure. The two localized separations witness the remaining factors 4 and 3. This evaluates both the candidate torsion and its reduction by the clock reader in rule L. [Rule L and Lemmas S and N](../../../docs/drafts/research-direction/explore-ec-ha/PROPOSAL.md#established-results)

## Reproduction

Run from the repository root. The first command checks the archived box table and direct normal stride controls. The next two compile in private directories and hold the shared computation lock. The stride run covers all six roots. The decrement run revisits the three even roots and only closes the canonical successors of period-66 arrivals.

```sh
python3 -B kaprekar-investigation/eleven_access/class_compat/premises.py | cmp - kaprekar-investigation/eleven_access/class_compat/premises.json
python3 -B kaprekar-investigation/eleven_access/class_compat/stride_inventory.py | cmp - kaprekar-investigation/eleven_access/class_compat/stride-inventory.json
python3 -B kaprekar-investigation/eleven_access/class_compat/prime_reader.py | cmp - kaprekar-investigation/eleven_access/class_compat/prime-reader.json
python3 -B kaprekar-investigation/eleven_access/class_compat/period66_witness.py | cmp - kaprekar-investigation/eleven_access/class_compat/period66-witness.out
```
