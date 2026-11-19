"""Replay the original archived SMT-LIB input with the two command-line solvers."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from itertools import zip_longest
from pathlib import Path
import gzip
import subprocess

from exhaust import write_json


def run(command, text):
    process = subprocess.run(command, input=text, capture_output=True, text=True, check=True)
    if process.stderr.strip():
        raise RuntimeError(process.stderr)
    return process.stdout


def archive(evidence, stem):
    smt_path = evidence / f"{stem}.smt2.gz"
    status_path = evidence / f"{stem}.statuses.gz"
    source = gzip.decompress(smt_path.read_bytes()).decode("utf-8")
    expected = gzip.decompress(status_path.read_bytes()).decode("utf-8").splitlines()
    queries = [number for number, line in enumerate(source.splitlines(), 1)
               if line.strip() == "(check-sat)"]
    if len(queries) != len(expected):
        raise ValueError((stem, len(queries), len(expected)))
    hashes = {p.name: sha256(p.read_bytes()).hexdigest() for p in (smt_path, status_path)}
    return source, expected, queries, hashes


def compare(expected, actual, query_lines):
    mismatches = [{"query": i, "input_line": line, "recorded": exp, "observed": obs}
                  for i, (exp, obs, line) in enumerate(zip_longest(expected, actual, query_lines), 1)
                  if exp != obs]
    unknowns = [i for i, status in enumerate(actual, 1) if status == "unknown"]
    return {"statuses": dict(Counter(actual)), "mismatches": mismatches, "unknowns": unknowns}


def replay(evidence, output, z3, yices):
    commands = {"Z3": [z3, "-smt2", "-in"], "Yices": [yices, "--incremental"]}
    versions = {"Z3": run([z3, "-version"], "").splitlines()[0],
                "Yices": run([yices, "--version"], "").splitlines()[0]}
    archives = {stem: archive(evidence, stem) for stem in ("positive", "zero")}
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {(stem, name): pool.submit(run, command, archives[stem][0])
                for stem in archives for name, command in commands.items()}
        output_data = {"versions": versions}
        for stem, (_, expected, queries, hashes) in archives.items():
            scores = {}
            for name in commands:
                statuses = jobs[stem, name].result().split()
                scores[name] = compare(expected, statuses, queries)
                print(f"{stem} {name}: {len(statuses)} results, "
                      f"{len(scores[name]['mismatches'])} mismatches, "
                      f"{len(scores[name]['unknowns'])} unknowns", flush=True)
            output_data[stem] = {"queries": len(queries), "recorded": dict(Counter(expected)),
                                 "sha256": hashes, "results": scores}
    write_json(output, output_data)
    return all(not value["mismatches"] and not value["unknowns"]
               for stem in archives for value in output_data[stem]["results"].values())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("evidence", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("z3")
    parser.add_argument("yices")
    args = parser.parse_args()
    raise SystemExit(not replay(args.evidence, args.output, args.z3, args.yices))
