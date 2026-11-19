"""Check the twelve-shift label formulas on both fresh base classes."""
from random import Random

LOW = (94040, 73770, 57940, 20, 200)
HIGH = (94130, 73850, 58015, 80, 300)
SHIFT = (0, -12, 18, -6, -18)


def add(q, w, scale=1):
    return tuple(x + scale * y for x, y in zip(q, w))


def chart(q, base):
    a, b, c, d, e = q
    return a, b, c, 2 * (b - e) - base, base - c - d


def inside(q, base):
    return all(l * base <= 100000 * x <= h * base for l, x, h in zip(LOW, chart(q, base), HIGH))


def labels(q, base):
    cohort = (q[1] - q[4]) % 3
    first, first_rem = divmod(2 * (q[1] - q[4]) - base + 1 - (-4, -2, -6)[cohort], 6)
    assert first_rem == 0
    aligned = add(add(q, (-5, -4, -3, -2, -1), first),
                  ((-5, 0, 1, 2, -1), (0, 0, 0, 0, 0), (-5, 2, 3, 4, -1))[cohort])
    slack = base - aligned[2] - aligned[3] - 1
    overshoot = 2 - slack % 2
    second = (slack + overshoot) // 2
    prepared = add(aligned, (-10, -2, 0, 2, -2), second)
    third, rem = divmod(prepared[0] + prepared[2] - 3 * (base - 1) // 2 + 4, 22)
    pinned = add(prepared, (-20, -4, -2, 2, -4), third)
    p, r = pinned[2], pinned[1]
    assert min(first, second, third) >= 0
    return (cohort, overshoot, rem - 4, r % 6, first, second, third), p, r


def main():
    random = Random(2416711)
    for base in (4200185, 4800233, 4300175, 5000159):
        ranges = [range((l * base + 99999) // 100000, h * base // 100000 + 1)
                  for l, h in zip(LOW, HIGH)]
        counts = [0, 0, 0]
        for _ in range(150):
            a, b, c, h, g = (random.choice(r) for r in ranges)
            h += (base - h) % 2
            q = a, b, c, base - c - g, b - (base + h) // 2
            q2 = add(q, SHIFT)
            if not (inside(q, base) and inside(q2, base + 12)):
                continue
            tag, p, r = labels(q, base)
            next_tag, next_p, next_r = labels(q2, base + 12)
            assert tag == next_tag
            assert (next_p - p, next_r - r) == (18, -12)
            # The switch's offset depends on the preserved prefix label.
            # It cancels in the J increment even without its numeric value.
            assert (24 + 2 * (next_p - p) - 4 * (next_r - r)) % 72 == 36
            counts[tag[0]] += 1
        assert min(counts) > 0
        print(base, base % 12, counts)


if __name__ == "__main__":
    main()
