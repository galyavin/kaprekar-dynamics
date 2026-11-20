# Finite parity forecast

For each order 3 through 20, use 512 fresh uniform draws with replacement from each integer cube ([0,L)^r) at (L=20{,}000) and (L=200{,}000). Both scales are held out from the archived parity panels. Orders 18, 19 and 20 have no archived panel in the assigned inputs. This forecast was informed by the earlier scales and orders. Use Python 3's `random.Random` in two independent reproducible streams per order, seeded with `17620261120 + 100*r + j`, where (j=1,2) indexes the two scales. Use `randrange(L)` per coordinate and do not replace draws that are censored.

The outcome is recurrence with zero sum of the nonnegative integer losses around one root cycle. Use exact state equality and at most 1,000,000 explicit map updates per draw, including steps used to reconstruct projected return words. Guarded skips may represent longer clocks. Mark cap exhaustion as unknown. Cubes for this test are uniform source measures. This panel tests no constraint on physical entrance, original word weight or linear lift budgets.

The following strict inequalities are forecast at **both** scales. The low threshold of 8% at orders 14 through 20 is a detection forecast, with no proposed odd/even distinction there. The earlier high even orders are included explicitly.

| Orders | Finite cube prediction |
|---|---|
| 3, 4, 5, 7, 8 | Above 10% |
| 9, 11, 13 | Above 2% |
| 6 | Below 6.5% |
| 10, 12 | Below 2.5% |
| 14, 15, 16, 17, 18, 19, 20 | Below 8% |

Use simultaneous confidence of at least 95% across the 36 panels by assigning (1/1440) to each of the 72 one-sided exact binomial tails. For an above-threshold prediction, test the binomial upper tail at the threshold using the confirmed successes. For a below-threshold prediction, test the lower tail using confirmed successes plus all unknowns. The opposite tail with the opposite worst-case completion can contradict a prediction. A panel is confirmed or contradicted only when the corresponding worst-case tail is at most (1/1440). All other panels are unresolved. The full forecast is held if each panel confirms, contradicted if any contradicts, and otherwise unresolved. Exact arithmetic determines tail comparisons.

The pointwise limiting conjecture for odd (r) asks for each fixed odd order a positive lower density of zero-loss sources with order-dependent fixed linear clock, loss and spread budgets. The even conjecture takes the cube-size limsup first and asks whether it tends to zero along even orders. A single fixed even order may have positive limiting lower density. This finite panel neither supplies nor requires an order-uniform lower bound on the odd orders.
