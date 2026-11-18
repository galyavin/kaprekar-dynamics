"""Replay the archived review, including the long complete-graph sweeps and sampling."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

source = Path(__file__).resolve().parent
out = source.parent / 'results' / 'review-replay'
out.mkdir(parents=True, exist_ok=True)


def run(*args, cwd=source):
    return subprocess.check_output(list(map(str, args)), cwd=cwd, text=True)


def record(name, output, compare=False):
    (out / name).write_text(output)
    if compare and output.split() != (source / name).read_text().split():
        raise ValueError(f'{name} differs from the archived result')
    print(name, 'passed', flush=True)


def require(value, message):
    if not value:
        raise ValueError(message)


def replay(work):
    for name in ['pairgraph', 'strings_exhaustive', 'mc_fresh', 'overshoot_mc', 'overshoot_fast']:
        subprocess.run([os.getenv('CC', 'cc'), '-O2', source / (name + '.c'),
                        '-o', work / name, '-lm'], check=True)
    # The archived acquisition order and whitespace can differ between runs.
    for name, lo, hi in [('sweep_30_1200.log', 30, 1200), ('sweep_1201_3000.log', 1201, 3000)]:
        output = run(sys.executable, 'sweep_pairgraph.py', work / 'pairgraph', lo, hi)
        record(name, output, compare=True)
    for name, bases in [('pairgraph_fresh.log', [997, 1500, 2401, 2999, 3000, 5003, 10007]),
                        ('pairgraph_10009.log', [10009])]:
        flags = [] if bases == [10009] else ['checks']
        record(name, '\n'.join(run(work / 'pairgraph', b, *flags) for b in bases), compare=True)
    record('strings_exhaustive.log',
           '\n'.join(run(work / 'strings_exhaustive', b) for b in range(4, 37)), compare=True)
    record('mc_fresh.log', '\n'.join(run(work / 'mc_fresh', b, 1000000, 7919*b + 17)
                                     for b in [997, 1500, 2401, 2999, 3000]), compare=True)
    record('fresh_table.md', run(sys.executable, source / 'summarize_fresh.py', cwd=out), compare=True)
    record('compare_branches.log', run(sys.executable, 'compare_branches.py', source.parent / 'source'),
           compare=True)
    z3 = run(sys.executable, 'z3_independent.py', '--probe')
    require('=== 43 checks, 0 failures ===' in z3 and '[FAIL]' not in z3, 'Independent Z3 check failed')
    record('z3_with_probes.log', z3)
    exact = run(sys.executable, 'limit_exact.py', 300)
    require('equal: True' in exact and 'paper interval contains ours: True' in exact
            and 'False' not in exact, 'Exact integration check failed')
    record('limit_exact.log', exact)
    record('tau_mc.log', run(sys.executable, 'tau_mc.py', 20000000), compare=True)
    rows = []
    for band in [(0.5, 1.5, 0.6, 0.9), (1, 2.5, 0.55, 0.95)]:
        rows.append(run(sys.executable, 'overshoot_limit.py', *band))
        for b, size in [(10007, 20000), (40009, 20000), (160001, 20000),
                        (160001, 100000), (640007, 100000), (2560009, 50000)]:
            # The archived logs omit their seeds. Replay the design with a stated seed.
            sample = [b, *band, size, 42]
            fast = run(work / 'overshoot_fast', *sample)
            if size == 20000:
                require(fast == run(work / 'overshoot_mc', *sample), f'Fast/slow disagreement at {b}')
            m = re.search(r'max \|e_B\| in \(44\) = ([\d.]+).*violations of \(41\): 0; high-branch mismatches: 0', fast)
            require(m and float(m[1]) <= 8, f'Overshoot check failed at {b}')
            rows.append(fast)
    record('overshoot_seed42.log', '\n'.join(rows))


if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as folder:
        replay(Path(folder))
