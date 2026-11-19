# Reproduction bundle

Companion to *Arithmetic selection of arrival times in eleven-digit Kaprekar dynamics* by Rustem Galyavin.

Contact: rustem@galyavin.com

Code is licensed under MIT. Text, figures, and data are licensed under CC BY 4.0. See [LICENSE.md](LICENSE.md).

## Reproduce the certificates

From a fresh copy of this directory, run

```sh
shasum -a 256 -c SHA256SUMS
bash reproduce.sh
```

The recorded run used Bash 3.2.57, Python 3.13.15, and Apple Clang 21.0.0 on arm64 macOS. The command uses Python's standard library, Bash, `patch`, a C compiler (`cc`), and a C++17 compiler with 128-bit integer support. Set `PYTHON` and `CXX` to select alternative executables. Keep Python assertions enabled. Temporary files are removed on exit. For a shared compute lock, set `KAPREKAR_COMPUTE_LOCK` to its absolute path. The command holds that lock throughout its run.

| Stage | Expected result |
| --- | --- |
| Parity cohort and even exclusion | The two guarded itineraries and the integer-subtraction witness replay. The even tree has 361,560 leaves. Its decoded leaf inventory has SHA-256 `1c62f809a5019babd80d429ec50f5a743fb9531ea064fee519c05c1bd33df63b`. The separate partition audit verifies density one, all 703 cross-lattice tests, the cut clocks, and 118 constant-state tails. |
| Six residue classes on the small gap rectangle | The C++ verifier reproduces 3,288,862 route leaves. Separate arithmetic checks each class's exact route weights and clock forms. A first-residue leaf dump matches the archived sample. |
| Second-box prefix, matchings, and residual routes | The two prefix records and both affine matching records agree byte for byte, with 4,576 matched labels per class. The residual record contains 6,336 entries, 50,196 affine leaves, and the guard cutoff 86,496,744. |
| Structural twelve-shift | The section calculation gives 1,584 canonical sections, a 30-update prefix increment, a 40-update checkpoint alignment, and a word clock difference of −10. The separate rational calculation gives the same translation. |
| Full second-box coverage and exact limit weights | The continuation closes 99 families in 6,555 leaves. Its 99,076 nonempty pullbacks all have base step 12. The full coverage cutoff for cohorts 1 and 2 is 86,496,843. Both Table 2 columns sum to 28,512 units and their signed difference has positive mass 68 units. |
| Locality and marked clocks | The two off-box pairs reproduce their recorded clocks. The marked weights reproduce from 96 selected residual leaves, and eight direct clock controls agree with the archived plateau values. |

The fresh copy passed in 762.71 seconds under the lock, including 569.14 seconds through the residual routes. The second-box guarded coverage took 71.45 seconds and the exact pullback weights took 92.74 seconds. An earlier replay through the residual routes took 464.31 seconds. [Fresh-copy log](results/reproduce.log) · [Earlier replay](results/baseline.log)

The inherited centered decrease, pairing condition, and cycle envelope are checked in the supplement to *Linear-time access for eleven-digit Kaprekar dynamics*. Its Appendix E gives the integer inequalities checked with Z3 5.4.1 and Yices 2.7.0. [Paper and supplement](https://doi.org/10.5281/zenodo.24358476)

## Separate checks

The `independent-review/` directory contains finite-check code and compact records. The earlier coverage, matching, and residual campaign summaries are in `independent-review/box/` beside the residual weight checks. The unstratified original-word sampling and literal digit subtraction check the seven graph weights and the local twelve-shift. The same command checks the parity cohort at a fresh base and the ramp on the small gap rectangle. Run it in a fresh copy with

```sh
python3.13 -B independent-review/run_direct.py
python3.13 -B independent-review/translation_controls.py
```

The recorded direct campaign took 274.37 seconds. The four main bases contribute 5,000 draws each, with 32 additional draws above the proof cutoffs. The twelve-shift agrees on all 424 matched inputs in both base classes and all three cohorts. [Direct campaign log](results/direct-recompute.log)

The independent residual replay reuses the primary symbolic route generator and checks chart extrema, fixed-base fiber weights, and 558 terminal states using full borrowed-digit sorting. Its exact outputs include the sharp fractional clock window used for the residual plateaus in Appendix D. The recorded separate replay took 346.57 seconds while the direct campaign was also running.

```sh
python3.13 -B independent-review/box/sharp_audit.py > sharp-rerun.json
cmp sharp-rerun.json independent-review/box/sharp-audit.json
```

To regenerate the residual leaf ledger and compare the separate signed-weight calculation, start in a fresh writable copy. The leaf driver uses the shared route generator and emits its arithmetic leaves. The weight script recomputes each fiber determinant and phase without importing that generator.

```sh
mkdir -p independent-review/box/results
cd box
PYTHONPATH=. python3.13 -B ../independent-review/box/leafdump_driver.py > ../independent-review/box/results/leaves.jsonl
cd ../independent-review/box
python3.13 -B weights.py
cmp results/weights.json weights.json
cmp results/chart_check.json chart_check.json
```

The earlier independent even, cohort, and cut-state experiments use `ww.c`, `cen.c`, and `spot.c`. In a fresh copy, run the following from `independent-review/`. The last command uses the dump from the next section.

```sh
for name in ww cen spot; do clang -std=c11 -O3 "$name.c" -o "$name"; done
bash even_campaign.sh
python3 -B cohort_check.py ./cen ./ww all
python3 -B forms_check.py ./ww 1201440 2000160 2998560
python3 -B terminal_check.py ../parity ./spot ./ww 300 20261119
python3 -B terminal_check.py ../parity ./spot ./ww 150 57 ../leaves.jsonl.gz
```

Use `locked.py LOCKFILE LOGFILE -- COMMAND` to give a separate experiment the same shared compute lock. The exact six-residue arithmetic is in `mixture/mixture.py` and the recorded clocks are in `mixture/tables.log`. The residual and full-coverage scripts and compact expected records are in `box/`.

## Regenerate a larger record

The archived bundle includes the compact records needed for the checks. It leaves out the full even leaf inventories, the full six-residue per-leaf listings, the 9.7 MiB residual leaf ledger, and the large raw per-draw and per-class campaigns for the second box, including the repeated-seed runs and point and sample archives. The selected marked leaves retain every ledger row used for the marked weights.

In a private copy, save the 361,560 even terminal classes with

```sh
python3 -B dump_tree.py ./parity leaves.jsonl.gz > exclusion-rerun.txt
cmp exclusion-rerun.txt parity/exclusion.txt
gzip -dc leaves.jsonl.gz | shasum -a 256
```

In a private copy, save the six-residue leaf inventory with

```sh
python3 -B mixture/dump_leaves.py 0 0 leaves-0.txt.gz > mixture-0.txt
```

The two integers select the first and last residue. Add `--all` before them to request every leaf for the chosen residues. For the full second-box pullback ledger, run `python3 -B limit_weights.py coverage.json --leaves` in `box/` of a private copy.
