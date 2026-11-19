# Independent inverse-basin reconstruction

Run `python3 -B independent/basin_reverse.py` from the supplement root. The calculation concerns the four-cycle at $B=5m$.

For each target, enumerate the number $h=1,\ldots,5$ of positive source differences. The leading $11-2h$ output digits are $B-1$. Solve for the remaining $2h$ sorted digits by exact rational elimination, imposing the target differences and the opposite digit sums $B-1+r_i$. The integer row $r$ ranges over zero-sum rows with entries from $-2$ through $2$, which includes every row allowed by the mirror-error bound in Appendix F. Discard any inconsistent row and any solution with out-of-range, nonintegral or unordered digits.

Pair the recovered digits by the sums implied by positional borrowing. Recover the ordered positive source differences, append zeros as needed and check each proposed predecessor by positional borrowing. Check every accepted edge again by forming a single integer subtraction value and extracting its base-$B$ digits. Repeat to exhaustion starting with the four cycle states.

For $m=20,21,50,1000,20001$, this yields 235 states with layers (4,27,31,38,57,36,25,17,0). The minimum full gaps (0,1,2,3) have multiplicities (141,48,30,16), giving weighted sum 391. There are 45 zero-innermost states at each base. At $m=20$, the reconstructed state-to-depth map matches the accompanying exhaustive enumeration in `basin_m20.txt`. These are finite checks of the affine certificate.
