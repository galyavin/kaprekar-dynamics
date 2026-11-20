# Eleven-digit stage graph from the regular chamber

**Task 097, outcome (b).** The first two macroscopic stages are exact on the entire uniform-\(\Pi\) source, up to null leading walls. The first branch that needs the full borrowed map is the third-stage section reached through an interlevel tie. Its limiting weight under uniform-\(\Pi\) entrance is **\(2007825833/29378565679380\), approximately \(6.83432\cdot 10^{-5}\)**, in all twelve classes. There is an exact seven-stage continuation on an open paired cell in every class. Finite global closure and a no-accumulation bound at nonzero radius remain open.

## Section and renewal calculation

Let \(q\) be the five ordered differences, \(v=(5,4,3,2,1)\), \(X=6q-Bv\), and \(x=X/B\). Keep the original normalized differences \(q^\circ\) as an input mark. The full positive-inner borrowed digits in units of six are \(6B-6\) and the five pairs
\[
 (X_i+v_iB-6\mathbf 1_{i=5},\;-X_{6-i}+v_iB-6\mathbf 1_{i<5}).
\]
A phase sorts the full list and takes symmetric differences, subtracting \(Bv\) from the output. It records the positive-inner guard, the comparisons between adjacent sorted entries, and the weak paired-level guards. The shorter borrowed list implements zero inner differences. When all levels stay ordered, the update agrees with the autonomous centered map. [Borrowed expansion](../../odd_gap_transport.md#exact-states-and-their-input-weights) · [Paired covariance](../../odd_base_covariance.md#an-autonomous-integer-map-beneath-the-base-parameter)

A section stores an affine lattice \(p+L\mathbb Z^r\) in centered coordinates. The implementation puts \(L\) in integer column echelon form and reduces \(p\) in the pivot rows. It tests a candidate phase word on the real face and exact affine sheet. The return must translate by \(d\) on the sheet, and its linear part must fix \(d\). For every comparison \(a(X)\) in every phase its value under the \(j\)-th return is then \(a(X)+j\,a_{\rm lin}d\). On the strict cell of a leading first winner, put \(s=-a_{\rm lin}d>0\). Comparisons of \(a/s\) with equal leading forms resolve by the exact affine intercepts. Given a remainder \(j\) of \(a\) modulo \(s\), the stage takes \(1+(a(X)-j)/s\) returns. The outgoing sheet is the image of the residue sublattice under this exact affine map. All arithmetic classes are compatible with the actual incoming image and with the stated base. [Implementation](graph.cpp)

The lattice calculation extends the renewal argument to the rank-preserving sections met below. Put \(h=\gcd(s,a_{\rm lin}L_1,\ldots,a_{\rm lin}L_r)\). There are \(s/h\) remainders \(j\) in one coset modulo \(h\), each of weight \(h/s\). For a fixed strict leading cell, the incoming mass has this split by Riemann sums on the cosets. Each outgoing lattice is the image of one of those cosets, with direction matrix
\[
 \left(I+\frac{d\,a_{\rm lin}}s\right)L_{\rm split}.
\]
It is exact on the split sheet even when \(d\) is transverse to its linear span. For subsequent residue tests, every compatible class in the image has an equal-index preimage in the source coset. Applying weighted Riemann sums to those preimages with any continuous test of the original position and leading clock retains the joint mark. This supplies the renewal statement required by the graph, including the case where the drift in the original Lemma R has to be allowed transverse to the incoming lattice. It also shows why testing independent full physical lattices at every section would give incorrect weights. [Lemma R and the clock reduction](../../../docs/drafts/research-direction/o2l-task088.md#2-pivot-and-bridge) · [Passage with a varying refinement](../../../docs/drafts/research-direction/astra-strategy-review-8.md#1-what-can-be-carried-from-eleven-digits)

## First two stages and first unpaired section

The regular return is five steps with centered displacement \((-30,-24,-18,-12,-6)\). In ordinary gap coordinates its displacement is \((5,-1,-1,-1,-1,-1)\). Over all five phases, the within-level rows fall by 36 in centered units and the interlevel rows rise by 30. The first exit is therefore positively paired. Of fifteen distinct within-level normals, exactly ten have positive-volume strict winner cells. The five singleton-gap normals have zero strict winner volume. The exact volume of \(\Pi\) in normalized difference coordinates is \(884293/727806643200\). The original ordered-word mass starting in \(\Pi\) without an entry prefix is about \(0.000806610935\). [First-face table and exact integrals](pi-geometry.json) · [Geometry check](pi_geometry.py)

The first-winner labels in the following table use the first canonical representative of a tied normal. Each displayed weight is for one of its five rotated phases. Residues \(r\) are in centered units, for the minimum-intercept comparison in the first return.

| First-winner class | Uniform-\(\Pi\) weight per phase | \(B\bmod12\) | First remainders \(r\) and weight per remainder | Second return |
| --- | ---: | --- | --- | --- |
| Level 2 comparison | \(3119581/44214650\) | All | \(0,6,12,18,24,30\), \(1/6\) each | 5 steps at \(r=18\), 10 otherwise |
| Level 3 comparison | \(5723349/44214650\) | \(0,2,4,6,8,10\) | \(6,18,30\), \(1/3\) each | 5 steps at \(r=18\), 10 otherwise |
| Level 3 comparison | \(5723349/44214650\) | \(1,3,5,7,9,11\) | \(0,12,24\), \(1/3\) each | 10 steps |

Every rank-four sheet above starts its second return immediately. The candidate affine return and all its primitive comparisons are verified throughout the entire first-winner cell by exact vertex minima. Splitting by the minimum of the declining guard quotients partitions that cell exactly in volume. Every second-stage drift is tangent to the incoming face, so each exit has rank three. The staged ledger includes the drift, winning leading time row, centered decrement, remainder set, and the probability of each remainder. Exact intercepts select the actual comparison when multiple rows give the same leading time. The comparisons with a nonzero \(B\) coefficient that win have compatible residue shift six. [Second-stage lane table in all twelve classes](sections.json) · [Exact phase and vertex calculation](sections.cpp) · [Volume and residue calculation](sections.py)

The first unpaired section is reached when a level-3 first winner is followed by a second-stage interlevel exit. Its weight conditional on any such first face is \(2007825833/19014487788834\), independent of the incoming remainder. The second-stage winners there include a weak physical level-order guard with positive decrement. All distinct leading competitors are strictly later in the interior of that winner cell. If a sorting comparison fails at equality before the weak physical guard fails, exchanging equal digits leaves the numerical update unchanged. The weak guard then fails at the same phase within one more return unless an interlevel failure occurs earlier. This identifies a first real post-\(\Pi\) pairing loss, as well as the first contact. No such contact through stage two occurs for a level-2 first winner. The leading portion of the third-stage entrance through this interlevel exit has the exact uniform-\(\Pi\) weight given above. Integrating the same strict cells against \(11!\prod_i g_i\) gives a limiting **original direct-entry mass of approximately \(4.66615569\cdot10^{-8}\)** in every class. It is a lower bound for the full original-weight population that ever loses pairing after \(\Pi\), which also contains trajectories with a nonzero entry prefix. [Exact contact masses](sections.json)

The first coverage gap in the all-\(\Pi\) atlas is the continuation from the rank-three sections. In particular, the interlevel exit requires a third-stage word that can alternate paired and unpaired phases. The accepted box route gives a concrete such word of length 20 with displacement \((-120,-24,-12,12,-24)\) in centered units. Its failing phase is zero and may recur at phase ten depending on the overshoot. The full-sort implementation admits this word. The box certificate and observed words give continuations on their covered cells. A finite exhaustive catalogue for all rank-three exits and their descendants, including the purely paired descendants, has yet to be proved. [Box third-stage return](../box_mod12/README.md#routes-from-the-first-block-exits) · [Whole-box failure clock](../box_structure/README.md#repeated-pairing-failures)

## Skeleton maps and speed lotteries

For any certified word of length \(p\), let \(D=a_{\rm lead}(x)\) on its current chart, including the base coefficient of the winning row. The update contributes \(\Delta S=(p/s)D\) and projects the leading position to \(x+(D/s)d\). The graph composes those affine maps. The successive remainder splits give the exact conditional distribution of the vectors \(p/s\). The total clock is the sum of their contributions, with zero leading cost for each bounded short switch. At an exit with zero leading state, the actual centered terminal tail is bounded on each certified finite affine sheet. The actual first cyclic-word time differs from the difference time by at most one, and the terminal mark is \(v/6\) after normalization. [Clock and envelope](../../eleven_centered_clock.md#3-global-capture-of-the-autonomous-physical-system) · [Explicit frozen law](PREDICTIONS.md)

For a fully enumerated example use the open centered cube of radius \(1/100000\) in maximum norm about
\[
 x_*=\frac{(2938,3225,2886,1913,387)}{50000}.
\]
All six tables certify strict radius \(1/50000\). This cube is inside the permanent source in the large-base limit and all calculated descendants are paired. The first four stages drop ranks from five to one. The fifth preserves rank one with conditional weight one in **each** class. A second rank-preserving stage occurs on some lanes. All paths reach rank zero by the seventh stage. This gives an explicit positive-weight counterexample to an argument based on strict rank decrease at every exit. Its conditional uniform-\(\Pi\) mass is at least \(350987/1036280859375000000000000000\), the probability of the certified cube. It leaves open a uniform bound on the number of rank-preserving stages. [Exact edge ledger and terminal affine forms](cell-lotteries.json) · [Cell verifier](lottery.cpp)

All terminal clock forms on this cube can be written \(S(x)+wD(x)\) with
\[
 S(x)=\frac{x_1}{6},\qquad D(x)=2x_2-2x_3+x_4-\frac45x_1 .
\]
The literal rows and weights in the ledger give the speeds by reading the coefficient of \(x_4\). The six distinct tables are indexed by \(B\bmod6\) and used for both compatible residues modulo 12. For example, the row for residue 5 or 11 is

| \(w\) | Weight |
| ---: | ---: |
| \(5/54\) | \(25/1296\) |
| \(5/64\) | \(31/1296\) |
| \(35/288\) | \(841/2592\) |
| \(55/459\) | \(11/324\) |
| \(35/648\) | \(4499/15552\) |
| \(15/272\) | \(2329/20736\) |
| \(65/1152\) | \(12241/62208\) |

The outgoing lattices merge by their current section and affine residue data. The ledger still retains the affine clock coefficients and weights of all leaves in each class. Covers of the other rank-three charts are a separate problem from this closed cell. [Exact all-class table](cell-lotteries.json) · [Reproduction driver](exact_cells.py)

## Prospective test and pairing loss

Both source laws were frozen with the exact-lottery sampler and independent schoolbook whole-word replay. At each of \(12000+r\) and \(48000+r\), \(0\le r<12\), the test used 8,000 genuine original words and 3,000 model draws. A direct uniform-\(\Pi\) control used 1,000 each. All model draws completed, with at most seven macroscopic stages, two rank-preserving stages and 251 short steps. This is a sampled coverage result. [Freeze and hashes](PREDICTIONS.md) · [Comparison ledger](holdout-results.json)

| Source and scale | Classes passing primary gates | Largest mean error | Largest CDF distance |
| --- | ---: | ---: | ---: |
| Original, \(12000+r\) | 12/12 | 0.004009 | 0.05602 |
| Original, \(48000+r\) | 12/12 | 0.000580 | 0.01840 |
| Uniform \(\Pi\), \(12000+r\) | 12/12 | 0.002605 | 0.059 |
| Uniform \(\Pi\), \(48000+r\) | 10/12 | 0.001391 | 0.059 |

All five position-and-clock moments and all tested paired terminal marks passed. The uniform-\(\Pi\) misses were the frozen 0.055 CDF bound at 48,001 and 48,007, with distances 0.059 and 0.057. A separately frozen larger follow-up with 5,000 draws per arm and model passed both, with distances 0.04292 and 0.04330. The two initial misses remain in the primary record. The six high-base original-weight shift comparisons passed, with largest actual CDF distance 0.00584 and mean difference \(7.664\cdot10^{-5}\). The model uses a shared lane stream across draws. An altered unpaired route can consume a different number of lane variates and alter later sampled paired lanes. The resulting finite model shift mismatch in class 3 was 0.00101 in CDF distance. [Frozen follow-up](FOLLOWUP.md) · [Follow-up record](followup-results.json)

On direct \(\Pi\) starts, a shift by six takes \(q\) to \(q+v\). It leaves \(X\) and the raw centered paired trajectory unchanged, and increases every physical adjacent-level slack by one. The lower-base paired orbit remains paired after this shift. On fixed strict leading entry prefixes of bounded length, one can pull this coupling back through the unimodular affine branch of the borrowed update. The resulting original displacement is bounded on the fixed prefix cell and preserves its limiting density. Thus any difference in the matched clocks after a six-shift requires a pairing failure in at least one arm, apart from prefix or guard collars. This is a localized covariance of the compared population. A separate coverage result is still required to pass to the bulk limit, and conditioning each arm on its own survivors can change the population being compared. [Covariance and its probability scope](../../odd_base_covariance.md#what-covariance-doesand-does-notimply-for-a-probability-law) · [Full leading map](../../odd_global_majorization.md#the-full-leading-map)

The actual original-weight sample recorded 25 losses in 96,000 draws at the smaller scale and 36 in 96,000 at the larger scale. In all 61 records the first failing guard was between levels 2 and 3, with no zero-inner failure. The first failing slack was \(-1\) or \(-2\). The large-base rates were \(3/8000\) in every class, with 15 to 82 failing updates through the replay. These pooled counts share random inputs across bases, so the 36 events at the larger scale represent three distinct underlying input draws, and the rates do not resolve a general mod-12 contrast. Forty-four high-base records had no observed \(\Pi\) entry before cycle detection and are charged separately. The archived box proves the signed class 5 minus 11 contrast **\(17/7128\)** conditional on its original weight \(8.51239\cdot10^{-14}\), which gives a localized contribution \(2.03017\cdot10^{-16}\) in the original weight. All 48 fresh box controls exhibited post-\(\Pi\) interlevel failures and landed on the frozen clock-test plateaus. [Actual failure rows and controls](holdout-results.json) · [Proved whole-box contrast](../box_structure/README.md#the-two-limit-laws)

## Reproduction

Verify the frozen sources and reproduce the two exact geometric stages from this folder. The second-stage run treats the source inputs as read-only and compiles in a private directory. The selected-cell and held-out campaigns also take the project compute lock.

~~~~sh
shasum -a 256 -c FREEZE.sha256
shasum -a 256 -c FOLLOWUP.sha256
python3 -B pi_geometry.py | cmp - pi-geometry.json
python3 -B sections.py | cmp - sections.json
python3 -B exact_cells.py | cmp - cell-lotteries.json
python3 -B holdout.py
python3 -B summarize.py /path/printed/by/holdout.py | cmp - holdout-results.json
python3 -B holdout.py --followup
python3 -B summarize.py /path/printed/for/followup --followup | cmp - followup-results.json
~~~~

The clipped-geometry helper is Task 088's `geometry.py`, with SHA-256 `2e53ea48fb81fc751936bc2760a10a18b08b158cce93a7dc241c6c835eb5d633`. The Task 095-owned box note was read via the stable coordinator snapshot, SHA-256 `b38f245530d1ef7e4d7450783507d83b1761447c9b58faab9ed7d967f40a76c2`.

