"""Direct-subtraction checks for the joint law; finite evidence, not its proof.

The twelve reference bases below were fixed after the joint formulas and proof
were written. Generate missing records with --holdout-executable.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess

from five_digit import (finite_graph, classical_cycles, multiplicity,
                        canonical_cycle)
from five_limit_checks import raw_step
from five_joint import (core_cycle_label, cycle_weights, joint_cdf,
                        rare_pairs, rare_arrival_counts, terminal_label)
from five_clock import core_region
from five_limit_law import cdf

ROOT = Path(__file__).resolve().parent
HOLDOUTS = (*range(431, 437), *range(1973, 1979))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--holdout-executable')
    args = parser.parse_args()
    checked = core = 0
    for base in range(30, 111):
        cycles, fate = finite_graph(base, raw_step)
        labels = {c: label for label, c in classical_cycles(base).items()}
        assert set(cycles) == set(labels)
        rare_hist = Counter()
        observed_rare = set()
        for q, (cycle, h) in fate.items():
            label = labels[cycle]
            assert terminal_label(base, q) == label
            checked += 1
            if core_region(base, q):
                core += 1
                assert core_cycle_label(base, q) == label
            if label == 'T2':
                observed_rare.add(q)
                rare_hist[h+1] += multiplicity(base, q)
        assert observed_rare == rare_pairs(base)
        if observed_rare:
            rare_hist[1] -= 2
            rare_hist[0] += 2
        assert dict(rare_hist) == rare_arrival_counts(base)
    print(f'Passed direct subtraction: {checked} pairs, {core} core labels, '
          'complete rare basins and actual histograms.', flush=True)

    thresholds = (Q(1, 25), Q(1, 6), Q(3, 10))
    limits = {t: cdf(t) for t in thresholds}
    for r in range(3):
        for t in thresholds:
            assert sum(joint_cdf(t, label, r) for label in
                       ('T1','T2','T3','T4','T5')) == limits[t]
    references = []
    for base in HOLDOUTS:
        path = ROOT/f'five_joint_holdout_{base}.json'
        if args.holdout_executable and not path.exists():
            data = json.loads(subprocess.check_output(
                [args.holdout_executable, str(base)], text=True))
            assert data['base'] == base
            path.write_text(json.dumps(data, indent=2)+'\n')
        data = json.loads(path.read_text())
        labels = {c: label for label, c in classical_cycles(base).items()}
        observed = {}
        for cycle in data['cycles']:
            c = canonical_cycle(tuple(tuple(q) for q in cycle['pairs']))
            hist = {int(t): w for t, w in cycle['input_arrival_counts'].items()}
            observed[labels[c]] = hist
        total = base**5-base
        assert data['total_input_weight'] == total
        assert sum(sum(h.values()) for h in observed.values()) == total
        assert observed.get('T2', {}) == rare_arrival_counts(base)
        row = {'base': base, 'pair_count': data['pair_count'], 'cycles': {}}
        for label, hist in observed.items():
            mass = sum(hist.values())
            pi = cycle_weights(base % 3).get(label, Q(0))
            row['cycles'][label] = {
                'basin_probability': float(Q(mass, total)),
                'limiting_probability': str(pi),
                'joint_cdf_minus_limit': {
                    str(t): float(Q(sum(w for h,w in hist.items()
                                       if Q(h,base) <= t), total)-pi*limits[t])
                    for t in thresholds},
                'conditional_mean_over_base': float(Q(
                    sum(t*w for t,w in hist.items()), mass*base))}
        references.append(row)
        print(f'Checked fresh reference B={base}: {sorted(observed)}; '
              'exact rare histogram agrees.', flush=True)
    result = {'direct_subtraction_pairs': checked, 'core_labels': core,
              'small_bases': [30, 110], 'references': references,
              'reference_pair_count': sum(r['pair_count'] for r in references),
              'scope': 'Core and rare formulas checked exactly against direct sorting; '
                       'finite joint CDF discrepancies are observations, not rate bounds.'}
    (ROOT/'five_joint_results.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
