"""Validate the five supplemental certificates and tally the completed census."""

from collections import Counter
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "arith"))
from census import validate


def read_rows(path):
    with path.open() as stream:
        return [json.loads(line) for line in stream]


root = Path(__file__).resolve().parent
initial = {row["n"]: row for row in read_rows(root / "arith/certificates.jsonl")}
supplemental = read_rows(root / "arith/planar-certificates.jsonl")
assert sorted(row["n"] for row in supplemental) == [811, 3343, 6967, 8167, 9931]
for row in supplemental:
    original = initial[row["n"]]
    assert original["cls"] == "unknown" and row["cls"] == "scarred"
    assert all(row[key] == original[key] for key in ("t", "r", "band", "bin"))
    validate(row)
    initial[row["n"]] = row

census = [row for row in initial.values() if row["band"] == "census"]
assert sorted(row["n"] for row in census) == list(range(7, 10_000, 2))
classes = ("strict", "scarred", "full-orbit", "zero")
expected = {
    "all": (2507, 1629, 701, 160),
    "primes": (226, 139, 701, 160),
    "composites": (2281, 1490, 0, 0),
    "index-two primes": (2, 45, 0, 68),
}
groups = {
    "all": census,
    "primes": [row for row in census if row["bin"].startswith("P")],
    "composites": [row for row in census if row["bin"].startswith("C")],
    "index-two primes": [row for row in census if row["bin"] == "P2"],
}
print("group: strict, scarred, full-orbit, zero")
for label, rows in groups.items():
    counts = Counter(row["cls"] for row in rows)
    actual = tuple(counts[name] for name in classes)
    assert sum(actual) == len(rows) and actual == expected[label], (label, counts)
    print(f"{label}: " + ", ".join(map(str, actual)))
