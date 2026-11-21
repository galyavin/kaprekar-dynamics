"""Recompute the finite measurements quoted in the paper."""
from collections import Counter
from csv import DictReader
from fractions import Fraction
from gzip import open as gzip_open
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "inputs/kaprekar-investigation/cross_length"


def check_original_panels():
    expected = ((111000, 100, 1428, 793), (1111000, 128, 1902, 826),
                (11110000, 128, 1866, 810))
    for scale, (anchor, per_base, off, rounded) in enumerate(expected, start=1):
        path = ROOT / "x1_terminal_map/raw" / f"n17-s{scale}.tsv.gz"
        with gzip_open(path, "rt") as stream:
            rows = list(DictReader(stream, delimiter="\t"))
        assert all(int(row["period"]) > 0 for row in rows)
        assert Counter(int(row["base"]) for row in rows) == {
            b: per_base for b in range(anchor, anchor + 18)
        }
        count = sum(1000 * int(row["radius"]) > 9 * int(row["base"]) for row in rows)
        assert count == off and round(Fraction(1000 * count, len(rows))) == rounded
        print(f"Early panel {scale}: {count}/{len(rows)} off-centre")


def check_holdout():
    panels = json.loads((ROOT / "ker17/holdout-counts.json").read_text())
    digests = ROOT.parents[2] / "KER17-PANEL-DIGESTS.tsv"
    with digests.open() as stream:
        digest_rows = [line.split() for line in stream]
    assert all(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest)
               for _, _, digest in digest_rows)
    digest_counts = {name: int(count) for name, count, _ in digest_rows}
    assert len(digest_counts) == len(digest_rows) == 108
    assert digest_counts == {
        f"{prefix}{base}.tsv": count
        for start in (90000, 360000, 1080000)
        for base in range(start, start + 18)
        for prefix, count in (("", 2048), ("literal-", 2))
    }
    for anchor, off, early in ((90000, 29264, 750), (1080000, 29690, 246)):
        panel = panels[str(anchor)]
        assert panel["diagnostics"]["cap"] == 0
        cells = panel["cells"]
        observed = Counter()
        for key, cell in cells.items():
            base = int(key.split(":")[0])
            assert 0 <= cell.get("early", 0) <= cell.get("off", 0) <= cell["n"]
            observed[base] += cell["n"]
        assert observed == {b: digest_counts[f"{b}.tsv"] for b in range(anchor, anchor + 18)}
        per_mark = tuple(sum(c.get(mark, 0) for c in cells.values()) for mark in ("n", "off", "early"))
        assert per_mark == (36864, off, early)
        print(f"Held-out panel {anchor}: {off}/36864 off-centre, {early}/{off} early")


def check_sections():
    expected = ((90000, 4096, 3225, (1393, 1663, 1702, 1639), 85),
                (1080000, 4096, 3264, (1262, 1510, 1516, 1446), 39),
                (10800000, 2048, 1652, (627, 759, 760, 738), 8))
    for base, count, off, episodes, early in expected:
        with gzip_open(ROOT / "s17_sect" / f"panel-{base}.tsv.gz", "rt") as stream:
            rows = [line.split() for line in stream]
        assert len(rows) == count
        assert all(int(row[0]) == base and int(row[6]) > 0 for row in rows)
        off_rows = [row for row in rows if 1000 * int(row[7]) > 9 * base]
        assert len(off_rows) == off
        first_absorb = tuple(sum(int(row[13 + 3 * j].split(",")[1]) == 1 for row in off_rows)
                             for j in range(4))
        assert first_absorb == episodes and first_absorb[0] < first_absorb[1] < first_absorb[2]
        assert sum(int(row[13].split(",")[1]) == 1 and int(row[5]) - int(row[8]) <= 1024
                   for row in off_rows) == early
        print(f"Section panel {base}: {off}/{count} off-centre; first-episode counts {first_absorb}")


if __name__ == "__main__":
    check_original_panels()
    check_holdout()
    check_sections()
