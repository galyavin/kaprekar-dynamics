"""Exact local maps about the Task 098 periodic witnesses."""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys


COMPARISONS = tuple(
    tuple(int(6 - i <= j) - int(j < i) for j in range(6))
    for i in (1, 2, 3)
)
IDENTITY = tuple(range(6))
BASIS = tuple(
    tuple(int(j == i) - int(j == 5) for j in range(6))
    for i in range(5)
)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def gaps(q):
    return tuple(x - y for x, y in zip((1, *q), (*q, 0)))


def word(g):
    return "".join("+" if d > 0 else "-" if d < 0 else "0"
                   for d in (dot(c, g) for c in COMPARISONS))


def sorted_step(g):
    q = [sum(g[j:]) for j in range(1, 6)]
    digits = sorted((1, *q, *(1 - x for x in q)), reverse=True)
    return gaps(tuple(digits[i] - digits[-1 - i] for i in range(5)))


def paired_step(g):
    d = (0, *(abs(dot(c, g)) for c in COMPARISONS))
    out = []
    for i in range(3):
        total = g[i] + g[5 - i]
        difference = d[i] - d[i + 1]
        out.extend(((total + difference) / 2, (total - difference) / 2))
    return tuple(out)


def primitive(word):
    return next(tuple(word[:p]) for p in range(1, len(word) + 1)
                if len(word) % p == 0 and word[:p] * (len(word) // p) == word)


def targets(rows):
    selected = []
    for row in rows:
        if not row["paired"]:
            continue
        w = primitive(row["word"])
        if row["dimension"] and not any("-" in phase for phase in w):
            continue
        selected.append((row["dimension"], row["period"], w,
                         gaps(tuple(map(Fraction, row["witness"])))))
    return sorted(selected)


def phase_options(w):
    return tuple(product(*(
        (left if left != "0" else right,) if left != "0" or right != "0"
        else ("+", "-")
        for left, right in zip("0" + w, w)
    )))


def reference_permutation(signs):
    return tuple(j for i, sign in enumerate(signs)
                 for j in ((i, 5-i) if sign == "+" else (5-i, i)))


def corrections(w, signs, frame):
    """Return oriented tie forms and pullback directions at one phase."""
    result = []
    for i, (left, right) in enumerate(zip("0" + w, w)):
        sign = 1 if signs[i] == "+" else -1
        for j, side in ((i, 1), (i+1, -1)):
            if j == 0 or w[j-1] != "0":
                continue
            normal = [0] * 6
            direction = [0] * 6
            for k, value in enumerate(COMPARISONS[j-1]):
                normal[frame[k]] = sign * value
            direction[frame[i]] = sign * side
            direction[frame[5-i]] = -sign * side
            result.append((j, tuple(normal), tuple(direction)))
    return result


def frame_schedule(w, signs):
    frame = IDENTITY
    schedule = []
    while True:
        for phase, side in zip(w, signs):
            schedule.append(corrections(phase, side, frame))
            frame = tuple(frame[j] for j in reference_permutation(side))
        if frame == IDENTITY:
            return schedule
        assert len(schedule) < 61 * len(w)


def positive_row(schedule):
    directions = [v for phase in schedule for _, _, v in phase]
    for t, phase in enumerate(schedule):
        for j, form, direction in phase:
            if (dot(form, direction) == 2
                    and all(dot(form, v) >= 0 for v in directions)):
                return t, j, form, direction
    return None


def orbit(witness, period):
    result = [witness]
    for _ in range(period):
        ref = result[-1]
        successor = sorted_step(ref)
        assert successor == paired_step(ref)
        result.append(successor)
    assert result[-1] == witness and len(set(result[:-1])) == period
    return result[:-1]


def verify(catalogue):
    rows = json.loads(catalogue.read_text())
    assert len(rows) == 191 and sum(r["paired"] for r in rows) == 30
    chosen = targets(rows)
    assert len(chosen) == 11
    certified = set()
    for dim, period, w, witness in chosen:
        phases = orbit(witness, period)
        assert tuple(map(word, phases)) == w
        centre = Fraction(1, 6)
        scaled = tuple(centre + (g - centre)/100 for g in witness)
        scaled_phases = orbit(scaled, period)
        assert all(min(g) > Fraction(1, 7) for g in scaled_phases)
        assert tuple(map(word, scaled_phases)) == w
        display = " ".join(w)
        print(f"{dim=} {period=} {display}: {', '.join(map(str, witness))}")
        for signs in product(*(phase_options(phase) for phase in w)):
            schedule = frame_schedule(w, signs)
            certificate = positive_row(schedule)
            if certificate:
                print(f"  rotating frame {signs}, {len(schedule)=}: {certificate}")
                certified.add(w)
                break
        else:
            print("  no monotone comparison in these reference frames")
    assert certified == {
        ("+++", "0--", "00+"),
        ("+++", "+++", "++0", "+0-", "++0", "+++", "0-0"),
        ("+++", "+++", "0--"),
    }


def local_offset_step(h, positive):
    """Leading offsets from ordered digits on a fixed sign branch."""
    q = tuple(sum(h[j:]) for j in range(1, 6))
    digits = [0]
    for i, sign in enumerate((*positive, positive[1], positive[0])):
        pair = q[i], -q[4 - i]
        digits.extend(pair if sign else pair[::-1])
    differences = tuple(digits[i] - digits[10 - i] for i in range(5))
    return tuple(x - y for x, y in zip((0, *differences), (*differences, 0)))


def verify_edge_returns():
    """Test each signed branch on all five basis offsets."""
    for first, second, extra in product((False, True), repeat=3):
        for h in BASIS:
            a, b, c, eta, _, f = h
            x, y, z = b-a, c-a, f-a
            s = 0 if first else -z

            # Canonical Two phase at zero anchor. The new tie is d3.
            u = -6*a - 2*x - 2*y if extra else 0
            t = 0 if second else -(y+s+u)
            first_out = local_offset_step(h, (first, False, extra))
            assert first_out[5] - first_out[0] == y + s + u
            k = local_offset_step(first_out, (second, True, True))
            assert (k[0], k[1]-k[0], k[2]-k[0], k[5]-k[0], k[3]) == (
                a-s-t, y+s+u+2*t, z+2*s+2*t, x-z+t, eta-u-t
            )

            # Other Two edge: in the return phase it is d2 that ties.
            w = eta-a
            t = 0 if second else -(y+s)
            r = 0 if extra else z-y-w
            first_out = local_offset_step(h, (first, False, False))
            assert dot(COMPARISONS[1], first_out) == y+w-z
            k = local_offset_step(first_out, (second, extra, True))
            assert (k[0], k[1]-k[0], k[2]-k[0], k[5]-k[0], k[3]-k[0]) == (
                a-s-t, y+s+2*t, z+2*s+2*t-r, x-z+t-r, w+s+r
            )
    print("Both boundary Two returns: 40 exact basis checks each")


def determinant(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return sum((-1)**j * value * determinant(
        tuple(row[:j] + row[j+1:] for row in matrix[1:])
    ) for j, value in enumerate(matrix[0]))


def verify_period_three():
    for first, second, third in product((False, True), repeat=3):
        for h in BASIS:
            a, b, c, d, e, f = h
            x, y, z, kappa = c-a, e-d, e+f-a-b, b+e-2*a
            r, s = (z, 0) if first else (0, -z)
            k, t = (x-r, 0) if second else (0, r-x)
            u = 0 if third else x-r-kappa
            j = local_offset_step(h, (True, first, False))
            assert j[5]-j[0] == x-r
            l = local_offset_step(j, (second, False, False))
            assert dot(COMPARISONS[1], l) == kappa-x+r
            result = local_offset_step(l, (True, third, True))
            aa, bb, cc, dd, ee, ff = result
            xp, yp, zp, kp = cc-aa, ee-dd, ee+ff-aa-bb, bb+ee-2*aa
            assert (xp, yp, zp, kp) == (
                x-r+2*t-u, -y+z+k, 2*r-y+t, kappa-y+z+k+u+2*t
            )
            df = (xp - yp + zp + kp) - (x - y + z + kappa)
            assert df == s + 5*t
            assert df - (xp-x) == r + s + 3*t + u
            assert 2*df + (kp-zp) - (kappa-z) == k + u + 11*t
    for second_level, centre, final_outer in product((False, True), repeat=3):
        for h in BASIS:
            a, b, c, d, e, f = h
            x, y, z, w = e-a, c-a, b-a, f-d
            s = 0 if second_level else w-y
            t = 0 if centre else w-y+z-x
            q = x-s+t
            u, v = (0, q) if final_outer else (-q, 0)
            j = local_offset_step(h, (True, True, True))
            l = local_offset_step(j, (True, second_level, centre))
            assert l[5]-l[0] == q
            aa, bb, cc, dd, ee, ff = local_offset_step(
                l, (final_outer, False, False)
            )
            xp, yp, zp, wp = ee-aa, cc-aa, bb-aa, ff-dd
            assert (xp, yp, zp, wp) == (
                y+s+u, x+z+2*u, q+2*u, w-s+v
            )
            assert (xp+yp+zp-wp) - (x+y+z-w) == 2*s + 4*u
            assert (xp+yp+zp+wp) - (x+y+z+w) == 4*u + 2*v
    a = ((1,-1,0), (-1,1,1), (0,1,0))
    b = ((1,-1,0), (0,1,-1), (1,-1,-1))
    assert abs(determinant(a)) == abs(determinant(b)) == 1
    rate = Fraction(2,3)**3 + Fraction(3,4)**3 + Fraction(1,2)**3
    assert rate == Fraction(1457,1728) < 1
    print(f"Period-three returns: 80 exact branch basis checks. Volume factor {rate}.")


def centre_branches():
    for outer, middle, centre in product((False, True), repeat=3):
        for h in BASIS:
            q = tuple(sum(h[j:]) for j in range(1, 6))
            a, b, c, d, e = -q[0], -q[1], q[2], -q[3], q[4]
            m, big = (a, e) if outer else (e, a)
            max_mid, max_reflected = (-b, -d) if middle else (d, b)
            abs_c = c if centre else -c
            out = local_offset_step(h, (outer, middle, centre))
            q_next = tuple(sum(out[j:]) for j in range(1, 6))
            assert q_next == (-m, -a-e, -big+max_mid, d-b, abs_c-max_reflected)
    print("Centre leading map: 40 exact basis checks.")


if __name__ == "__main__":
    verify(Path(sys.argv[1]))
    verify_edge_returns()
    verify_period_three()
    centre_branches()
