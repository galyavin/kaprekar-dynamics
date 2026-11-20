"""Active comparison forms in the five-step folded rotating frame."""
import json
from pathlib import Path

P = (0, 5, 1, 4, 2, 3)
FORMS = (
    (-1, 0, 0, 0, 0, 1),
    (-1, -1, 0, 0, 1, 1),
    (-1, -1, -1, 1, 1, 1),
)
DIRECTIONS = (
    (-1, 1, 0, 0, -1, 1),
    (0, -1, 1, -1, 1, 0),
    (0, 0, -1, 1, 0, 0),
)


def dot(row, v):
    return sum(a * b for a, b in zip(row, v))


def active_at_boundary(word):
    assert len(word) == 5 and all("-" not in phase for phase in word)
    index = tuple(range(6))
    active = []
    for phase, signs in enumerate(word):
        for j, sign in enumerate(signs):
            if sign == "0":
                a = tuple(FORMS[j][index.index(i)] for i in range(6))
                v = tuple(DIRECTIONS[j][index.index(i)] for i in range(6))
                active.append((phase, j, a, v))
        index = tuple(index[i] for i in P)
    assert index == tuple(range(6))
    return active


def main():
    rows = json.loads(Path("catalogue.json").read_text())
    for row in sorted(rows, key=lambda r: (-r["dimension"], r["word"])):
        word = row["word"]
        if not row["paired"] or len(word) != 5 or any("-" in s for s in word):
            continue
        if row["dimension"] == 5:
            continue
        active = active_at_boundary(word)
        cross = [[dot(a, v) for _, _, _, v in active] for _, _, a, _ in active]
        print(row["dimension"], " ".join(word), [(p, j + 1) for p, j, _, _ in active], cross)


if __name__ == "__main__":
    main()
