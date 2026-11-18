"""Reproduce the finite computations in Arrival times for the four-digit Kaprekar map."""
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
from time import monotonic
import hashlib
import json
import os
import subprocess
import sys

from independent import run as independent_checks

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT/'results'
sys.path.insert(0,str(ROOT/'source'))

from arithmetic_scaling_checks import four_checks
from arrival_checks import direct_input_histogram
from arrival_limits import limiting_law
from cycle_model import weighted_cycle_layers
from investigate import published_maximum_distance
from modular_checks import cycles
from symbolic_arrivals import symbolic_distribution


def compiled(source,work):
    compiler = os.environ.get('CXX','c++')
    program = Path(work)/source.stem
    subprocess.run([compiler,'-std=c++17','-O2',str(source),'-o',str(program)],check=True)
    return program


def section_eight(work):
    comparisons = 0
    for m in range(5,32,2):
        targets = list(map(set,cycles(m)))
        if len(targets)>1:
            targets.append(set.union(*targets))
        for n in range(1,4):
            for target in targets:
                lifted = {(u<<n,v<<n) for u,v in target}
                observed = weighted_cycle_layers(m<<n,lifted)
                actual = symbolic_distribution(m,n,target)
                assert observed['input_arrival_counts'] == actual['input_arrival_counts']
                comparisons += 1
    assert comparisons == 300

    inputs = []
    for m,n,which in [(5,0,'all'),(5,1,'all'),(7,0,'all'),(7,1,'all'),
                      (13,0,'first'),(13,0,'all'),(7,1,'empty')]:
        parts = list(map(set,cycles(m)))
        target = {'all':set.union(*parts),'first':parts[0],'empty':set()}[which]
        base = m<<n
        lifted = {(u<<n,v<<n) for u,v in target}
        hist = direct_input_histogram(base,lifted)
        assert hist == weighted_cycle_layers(base,lifted)['input_arrival_counts']
        inputs.append(dict(m=m,n=n,target=which,eligible=base**4-base,histogram=hist))
    assert sum(row['eligible'] for row in inputs) == 146904

    program = compiled(ROOT/'source/arrival_holdout.cpp',work)
    large = json.loads(subprocess.check_output([program,'5','10','fixed'],text=True))
    claimed = symbolic_distribution(5,10,{(3,1)})
    assert large['checked_pairs'] == 13109759
    assert {int(k):v for k,v in large['input_arrival_counts'].items()} == claimed['input_arrival_counts']
    for key in ('basin_input_count','complement_input_count','first_moment_sum','second_moment_sum'):
        assert large[key] == claimed[key]
    maxima = []
    for n in range(1,41):
        value = symbolic_distribution(5,n,{(3,1)})['maximum_finite_input_time']
        assert value == published_maximum_distance(5<<n)
        maxima.append(value)
    limits = [limiting_law(5,r,{(3,1)}) for r in range(4)]
    uniform = four_checks()
    assert all(Fraction(row['tv'])*row['odd_part']<26 for row in uniform['finite_distributions'])
    assert len(uniform['finite_distributions']) == 20
    return dict(comparisons=comparisons,inputs=inputs,base5120=large,
                maxima=maxima,limits=limits,uniform=uniform)


def record(name,data):
    raw = (json.dumps(data,default=str,indent=2,sort_keys=True)+'\n').encode()
    path = RESULTS/name
    if path.exists() and path.read_bytes() != raw:
        raise ValueError(f'Recomputed values differ from {path.name}')
    if not path.exists():
        path.write_bytes(raw)


def manifest():
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(ROOT.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts
            and p != RESULTS/'artifact_sha256.json' and not p.name.startswith('.')}


def main():
    if not __debug__:
        raise RuntimeError('Run with Python assertions enabled.')
    RESULTS.mkdir(exist_ok=True)
    start = monotonic()
    with TemporaryDirectory() as work:
        evidence = section_eight(work)
        print(f"Section 8 passed in {monotonic()-start:.2f} s",flush=True)
        results = independent_checks(compiled(ROOT/'graph_check.cpp',work))
    for paper,independent in zip(evidence['limits'],results['table']):
        assert [paper['atoms'].get(j,0)*1125 for j in range(1,6)] == independent['numerators']
        assert str(paper['basin_limit']) == independent['basin']
        assert [str(paper[key]) for key in ('conditional_mean_coefficient',
                'conditional_second_moment_coefficient')] == independent['moments']
    record('section8.json',evidence)
    record('independent.json',results)
    record('artifact_sha256.json',manifest())
    print(f"All checks passed in {monotonic()-start:.2f} s")
    print(sys.version.splitlines()[0])
    print(subprocess.check_output([os.environ.get('CXX','c++'),'--version'],text=True).splitlines()[0])


if __name__ == '__main__':
    main()
