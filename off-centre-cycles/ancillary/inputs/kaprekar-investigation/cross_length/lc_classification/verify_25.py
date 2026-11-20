"""Check the forced-25 transfer and its exact constants by literal subtraction."""
import random

from verify_root import literal, roots

N = 25
M = 13
K = frozenset((5, 10))
J = roots(N)
R = len(J)-1


def fold(j):
    return min(2*j, N-2*j)


def wall(g, i):
    return sum(g[-i:]) - sum(g[:i])


def diff(g):
    return tuple(sum(g[i:]) for i in range(1, M))


def forced(g):
    z = wall(g, 1)
    e = (z > 0) - (z < 0)
    h1, h2, h6 = (max(0, -1-wall(g,i)) for i in (1,2,6))
    nxt = [0]*M
    for j, x in enumerate(g):
        nxt[fold(j)] = x
    nxt[0] += e-h1
    nxt[1] += -2*e+h1
    nxt[2] += e+h1-h2
    nxt[3] += -h1+h2
    nxt[4] += h2
    nxt[5] -= h2
    nxt[10] -= h6
    nxt[11] += h6
    assert sum(nxt) == sum(g)
    assert nxt[5]+nxt[10] == g[5]+g[10]-h2-h6
    assert walls_profile == (0,0,1,1,1,0)
    if not h2:
        assert h1*(-g[0]+g[12]+g[1]-g[11]) <= 0
    return tuple(nxt), tuple(h > 0 for h in (h1,h2,h6))


walls_profile = tuple(wall(tuple(int(j in K) for j in range(M)), i) for i in range(1,7))


def main():
    assert R == 10 and set(J) == set(range(M))-K
    assert fold(5) == 10 and fold(10) == 5 and fold(11) == 3
    max_kick = R+2
    assert (constant := 14*max_kick*R*(R-1)+4*max_kick**2*R) == 20880
    assert (root_radius := constant//2+2*R*max_kick) == 10680
    assert 1+10*root_radius == 106801
    mask_set = set()
    grid = random.Random(251610)
    for trial in range(2048):
        shift = trial % 7
        g = tuple(10_000 + (1000 if j in K else 0) + shift + grid.randrange(-35, 36)
                  for j in range(M))
        b = sum(g)
        assert (M+1)*min(g) >= b+1 and all(wall(g,i) >= -1 for i in (3,4,5))
        image, mask = forced(g)
        mask_set.add(mask)
        assert diff(image) == literal(diff(g), b)
    assert len(mask_set) == 8
    for bump in (1, 99, 401):
        g = tuple(10_000 + (bump if j in K else 0) for j in range(M))
        assert forced(g)[0] == g and literal(diff(g), sum(g)) == diff(g)
    expected = ((11,11,11),(10,12,11),(9,11,11),(8,10,13),(7,9,11),
                (6,8,9),(5,8,7),(4,5,3),(3,4,4),(3,4,4))
    fixed = (10000, *(10100 if j in K else 10000 for j in range(1,M)))
    ray = (-10, *(0 if j in K else 1 for j in range(1,M)))
    predecessor = tuple(x+y for x,y in zip(fixed, ray))
    offsets = ray
    source = [next(i for i in range(M) if fold(i) == j) for j in range(M)]
    image = predecessor
    for comparisons in expected:
        assert tuple(wall(image,i) for i in (1,2,6)) == comparisons
        assert tuple(wall(offsets,i) for i in (1,2,6)) == (11,11,11)
        b = sum(image)
        nxt, mask = forced(image)
        assert not any(mask) and literal(diff(image),b) == diff(nxt)
        image = nxt
        offsets = tuple(offsets[i] for i in source)
    assert image == fixed and offsets == ray
    print("ten-step feeder, affine comparison slopes 11 and exact terminal family")
    print("2048 integer forced-tube transitions, all eight loss masks, three composite fixed states")


if __name__ == "__main__":
    main()
