"""Verify the order-12 root certificate and its all-base lifts."""
import argparse
from fractions import Fraction
import importlib.util
import json
from math import ceil, lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
TASK164 = HERE.parent / "lc_classification"
LIFTS = {45: {3: 1, 5: 1, 9: 2, 15: 2}, 91: {3: 1, 7: 1, 9: 1, 13: 1}}


def folded(length, index):
    return min(2*index, length-2*index)


def orbits(length):
    remaining = set(range(length//2+1))
    while remaining:
        orbit = []
        index = min(remaining)
        while index in remaining:
            remaining.remove(index)
            orbit.append(index)
            index = folded(length, index)
        assert index == orbit[0]
        yield tuple(orbit)


def margins(state):
    return tuple(sum(state[-i:])-sum(state[:i]) for i in range(2,len(state)//2+1))


def suffixes(state):
    total = 0
    suffix = []
    for item in reversed(state[1:]):
        total += item
        suffix.append(total)
    return tuple(reversed(suffix))


def lift(length, budget):
    root = None
    orbit_sets = tuple(orbits(length))
    m = length//2+1
    weights = [0]*m
    for orbit in orbit_sets:
        if orbit[0] == 1:
            root = (0, *orbit)
        for index in orbit:
            weights[index] = LIFTS[length].get(orbit[0], 0)
    assert root is not None and len(root) == 13 and 3 not in root
    assert set(LIFTS[length]) == {orbit[0] for orbit in orbit_sets} - {0, 1}
    mean = Fraction(sum(weights), m)
    epsilon = Fraction(1,2*m*(m+1))/(mean-min(weights))
    profile = tuple(Fraction(1,m)+epsilon*(w-mean) for w in weights)
    assert sum(profile) == 1
    assert all(profile[i] == profile[folded(length,i)] for i in range(m))
    spare = min(profile)-Fraction(1,m+1)
    margin = min(margins(profile))
    amplitude = max(abs(x-Fraction(m-i,m)) for i,x in enumerate(suffixes(profile),1))
    assert min(spare, margin, amplitude) > 0
    error = m+2*budget+1
    room = min(spare, margin/m, amplitude/m)
    d = ceil(4*error/room)
    assert Fraction(error,d) <= room/4
    print(f"strict length={length} margin={margin} delta={amplitude/2} D={d} "
          f"eta={length}!/({2*(m+1)}^{m}*{d}^{m-1})")
    return root, profile, amplitude/2, d


def screen():
    order_twelve = []
    for n in range(5,102,2):
        orbit_by_representative = {orbit[0]: orbit for orbit in orbits(n)}
        if len(orbit_by_representative[1]) != 12:
            continue
        candidates = tuple(orbit for rep,orbit in orbit_by_representative.items() if rep not in (0,1))
        zero_rows = tuple(i for i in range(2,(n//2+1)//2+1)
                          if all(margins(tuple(j in orbit for j in range(n//2+1)))[i-2] == 0
                                 for orbit in candidates))
        if n in LIFTS:
            assert not zero_rows
        else:
            assert zero_rows
        order_twelve.append((n,zero_rows))
    assert order_twelve == [(35,(2,)),(39,(4,)),(45,()),(91,())]
    print("order twelve through 101:", order_twelve)


def physical(root_module, record, lift_data):
    """Independently subtract whole words at stage boundaries and every terminal phase."""
    capture, stop = record["cycle"]
    times = set(range(12)) | set(range(capture,stop+1))
    for s,t in record["insertions"]:
        times.update((s-1,s,t-1,t))
    retirement = record["insertions"][-1][1]
    times.update((retirement-1,retirement,retirement+1))
    samples = {n:{} for n in lift_data}
    x = tuple(record["source"])
    models = {}
    for n,(root,profile,delta,d) in lift_data.items():
        m = len(profile)
        base = lcm(*(v.denominator for v in profile)) * (2*d)
        initial = [int(base*p) for p in profile]
        for j,value in zip(root[1:],x):
            initial[j] += value
        initial[max(set(range(m))-set(root))] -= sum(x)
        root_sum = sum(initial[j] for j in root)
        inverse = [0]*m
        for j in range(m):
            inverse[folded(n,j)] = j
        old_indices = [tuple(range(m))]
        for _ in range(12):
            old_indices.append(tuple(old_indices[-1][inverse[j]] for j in range(m)))
        assert old_indices[-1] == old_indices[0]
        models[n] = root, initial, root_sum, old_indices, base, delta
    loss = 0
    for t in range(stop+1):
        if t in times:
            for n,(root,initial,root_sum,old_indices,base,delta) in models.items():
                m = n//2+1
                g = [initial[j] for j in old_indices[t%12]]
                g[0], residue = divmod(root_sum+loss-sum(x), len(root))
                assert residue == 0
                for j,value in zip(root[1:],x):
                    g[j] = g[0]+value
                assert sum(g) == base and (m+1)*min(g) > base+1 and min(margins(g)) > 0
                q = suffixes(g)
                nxt = root_module.gap_step(g,n)
                assert root_module.literal(q,base) == suffixes(nxt)
                assert root_module.step(x)[0] == tuple(nxt[j]-nxt[0] for j in root[1:])
                radius = max(abs(Fraction(value,base)-Fraction(m-i,m)) for i,value in enumerate(q,1))
                assert radius > delta
                samples[n][t] = tuple(g)
        x, h = root_module.step(x)
        loss += h
        if h and t < stop:
            for root,initial,root_sum,old_indices,base,delta in models.values():
                initial[old_indices[(t+1)%12][3]] -= h
    for n, observed in samples.items():
        assert observed[capture] == observed[stop]
        assert len({observed[t] for t in range(capture,stop)}) == stop-capture
        print(f"literal length={n} base={models[n][4]} comparisons={len(observed)} "
              f"least physical period={stop-capture}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task164", type=Path, default=TASK164)
    parser.add_argument("--literal", action="store_true")
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("verify_root", args.task164 / "verify_root.py")
    root_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(root_module)

    old = next(record for line in (args.task164 / "root-certificates.jsonl").read_text().splitlines()
               if (record := json.loads(line))["order"] == 13)
    root_module.verify(old)
    record = json.loads((HERE / "certificate.jsonl").read_text())
    assert record["order"] == 12
    _, rays, cycle, budget = root_module.verify(record)
    assert all(v >= 0 for ray in rays for v in ray)
    assert len(cycle) == record["cycle"][1]-record["cycle"][0] == 24
    screen()
    lift_data = {n:lift(n,budget) for n in LIFTS}
    if args.literal:
        physical(root_module,record,lift_data)


if __name__ == "__main__":
    main()
