"""Reduce the frozen paired run using only the Python standard library."""
from collections import Counter
import hashlib
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist
import sys

N = 60_000
BASES = [scale + r for scale in (12000, 48000, 192000) for r in (3, 5, 15, 17)]
Z = NormalDist().inv_cdf(0.975)
SIM = NormalDist().inv_cdf(1 - 0.05 / (2 * 22))
ONE_THIRD = 0.00015


def wilson(n, total=N):
    p = n / total
    denominator = 1 + Z**2 / total
    center = (p + Z**2 / (2 * total)) / denominator
    half = Z * sqrt(p * (1 - p) / total + Z**2 / (4 * total**2)) / denominator
    return dict(count=n, denominator=total, fraction=p,
                Wilson95=[max(0, center-half), min(1, center+half)])


def estimate(values, selected):
    n = sum(selected)
    assert n > 0
    mean = sum(v for v, i in zip(values, selected) if i) / n
    return mean, [(v-mean)*N/n if i else 0 for v, i in zip(values, selected)]


def combine(*terms):
    mean = sum(k * m for k, (m, _) in terms)
    influence = [sum(k * v for (k, _), v in zip(terms, row))
                 for row in zip(*(i for _, (_, i) in terms))]
    return mean, influence


def interval(estimate, z=Z):
    point, influence = estimate
    se = sqrt(sum(a*a for a in influence)/(N*(N-1)))
    return dict(estimate=point, SE=se, halftime=z*se, CI=[point-z*se, point+z*se])


def cohort_summary(left, right):
    masks = zip(left["complete"], right["complete"], left["capture"], right["capture"],
                left["entry"], right["entry"], left["lost"], right["lost"],
                left["model_loss"], right["model_loss"])
    chart = Counter()
    match, paired = [], []
    for ml, mr, cl, cr, el, er, pl, pr, gl, gr in masks:
        exclusion = (int(not (ml and mr)) + 2*int(not (cl and cr)) +
                     4*int(not (el and er)) + 8*int(pl or pr) + 16*int(gl or gr))
        chart[exclusion] += 1
        match.append(exclusion & 7 == 0)
        paired.append(exclusion == 0)
    return match, paired, {
        "retained_entry": wilson(sum(match)),
        "retained_without_later_loss": wilson(sum(paired)),
        "exclusion_bit_counts": {key: wilson(n) for key, n in sorted(chart.items())},
    }


def compare(left, right):
    match, paired, output = cohort_summary(left, right)
    details = {}
    for name, mask, zscore in [("entry", match, SIM), ("entry_no_loss", paired, Z)]:
        physical = combine((1, estimate(left["actual"], mask)),
                           (-1, estimate(right["actual"], mask)))
        model = combine((1, estimate(left["graph"], mask)),
                        (-1, estimate(right["graph"], mask)))
        gap = combine((1, physical), (-1, model))
        details[name] = {"physical": physical, "model": model, "gap": gap}
        output[name] = {k: interval(v, zscore if k=="gap" else Z) for k, v in details[name].items()}
    all_rows = [True]*N
    uncond_mask = [a and b for a, b in zip(left["capture"], right["capture"])]
    output["all_physical_captured"] = interval(combine(
        (1, estimate(left["actual"], uncond_mask)),
        (-1, estimate(right["actual"], uncond_mask))))
    for arm in (left, right):
        arm["own_mask"] = [m and c and e for m, c, e in
                           zip(arm["complete"], arm["capture"], arm["entry"])]
    for key in ("actual", "graph"):
        scores = [[x if i else 0 for x, i in zip(arm[key], arm["own_mask"])]
                  for arm in (left, right)]
        output[f"{key}_original_denominator_scores"] = interval(combine(
            (1, estimate(scores[0], all_rows)), (-1, estimate(scores[1], all_rows))))
    return output, details, match


def read_data(root, count):
    by_base = [[] for _ in BASES]
    hashes = {}
    index = 0
    leading = Counter()
    for path in sorted(Path(root).glob("rows-*.txt")):
        digest = hashlib.sha256()
        with path.open("rb") as file:
            for line in file:
                digest.update(line)
                row = line.split()
                assert int(row[0]) == index and len(row) == 98
                leading[int(row[1])] += 1
                index += 1
                for k, dataset in enumerate(by_base):
                    start = 2+8*k
                    status, model_loss, graph, actual, strict, strict_loss, five, five_loss = row[start:start+8]
                    dataset.append((int(status), int(model_loss), float(graph), int(actual),
                                    int(strict), int(strict_loss), int(five), int(five_loss)))
        hashes[path.name] = digest.hexdigest()
    assert index == count
    return by_base, hashes, leading


def to_columns(dataset, base, reader):
    status, graph_loss, graph, clock, strict, strict_loss, five, five_loss = zip(*dataset)
    enter, loss = (strict,strict_loss) if reader=="strict" else (five,five_loss)
    physical = dict(
        no_entry=wilson(sum(t<0 for t in enter)),
        late_entry=wilson(sum(50*t>base for t in enter)),
        interlevel_loss=wilson(sum(bool(i&1) for i in loss)),
        zero_inner_loss=wilson(sum(bool(i&2) for i in loss)),
        both_losses=wilson(sum(i==3 for i in loss)),
        uncaptured=wilson(sum(t<0 for t in clock)),
        loss_union=wilson(sum(i>0 for i in loss)),
        leading_entered_but_physical_never_entered=wilson(sum(
            s!=5 and t<0 for s,t in zip(status,enter))),
    )
    model = dict(
        statuses={i: wilson(n) for i,n in sorted(Counter(status).items())},
        any_loss=wilson(sum(graph_loss)),
        complete_with_loss=wilson(sum(s==0 and l for s,l in zip(status,graph_loss))),
    )
    columns = dict(actual=[t/base for t in clock], graph=graph,
                   complete=[s==0 for s in status], capture=[t>=0 for t in clock],
                   entry=[t>=0 for t in enter], lost=[i>0 for i in loss], model_loss=graph_loss)
    return columns, dict(physical=physical, model=model)


def main(directory):
    by_base, hashes, leading = read_data(directory, N)
    result = dict(independent_vectors=N, bases=BASES,
                  nominal_simultaneous_z=SIM, sha256_rows=hashes,
                  leading_entry=dict(no_entry_by_cap=wilson(leading[-1])),
                  readers={})
    shrinking, persistent, widths, shrink_direction, persistence_direction = [], [], [], [], []
    for reader in ("strict", "five"):
        bundles = [to_columns(data, base, reader) for data, base in zip(by_base,BASES)]
        arms = [item[0] for item in bundles]
        book = dict(bases={base: item[1] for base,item in zip(BASES,bundles)}, scales={})
        result["readers"][reader] = book
        adjacent = []
        high_pairs = None
        for j in range(3):
            low_pair, low_effects, low_mask = compare(*arms[4*j:4*j+2])
            high_pair, high_effects, high_mask = compare(*arms[4*j+2:4*j+4])
            pooled_effects = {cohort: {k: combine((0.5, low_effects[cohort][k]),
                                                 (0.5, high_effects[cohort][k]))
                                      for k in low_effects[cohort]} for cohort in low_effects}
            pooled = {c: {k: interval(value, SIM if c=="entry" and k=="gap" else Z)
                          for k,value in effects.items()}
                      for c,effects in pooled_effects.items()}
            both = [a and b for a,b in zip(low_mask,high_mask)]
            shifts = {}
            for name, offset in (("class3",0),("class5",1)):
                left, right = arms[4*j+offset], arms[4*j+2+offset]
                shifts[name] = {k: interval(combine((1,estimate(right[k],both)),
                                                   (-1,estimate(left[k],both))))
                                for k in ("actual","graph")}
            scale = BASES[4*j]-3
            book["scales"][scale] = dict(
                pair_3_5=low_pair, pair_15_17=high_pair, equally_weighted=pooled,
                both_pair_cohorts=wilson(sum(both)),
                either_pair_cohort=wilson(sum(a or b for a,b in zip(low_mask,high_mask))),
                high_minus_low_gap=interval(combine((1,high_effects["entry"]["gap"]),
                                                    (-1,low_effects["entry"]["gap"]))),
                twelve_shifts_high_minus_low=shifts,
            )
            adjacent.append(pooled_effects["entry"]["gap"])
            widths.extend(record["entry"]["gap"]["halftime"]
                          for record in (low_pair,high_pair,pooled))
            high_pairs = [low_effects["entry"]["gap"][0],high_effects["entry"]["gap"][0]]
        changes = [interval(combine((1,adjacent[i]),(-1,adjacent[i-1])),SIM) for i in (1,2)]
        book["pooled_gap_changes"] = changes
        widths.extend(change["halftime"] for change in changes)
        midpoint = [i[0] for i in adjacent]
        sign = 1 if midpoint[0]>=0 else -1
        shrink_direction.append(sign)
        shrinking.append(
            abs(midpoint[2])<=abs(midpoint[0])/3 and
            any(
                0<=sign*midpoint[i]<sign*midpoint[i-1] and
                max(sign*v for v in changes[i-1]["CI"])<0
                for i in (1,2)))
        sign = 1 if midpoint[1]>=0 else -1
        persistence_direction.append(sign)
        persistent.append(
            all(min(sign*v for v in interval(estimate,SIM)["CI"])>0.0003
                for estimate in adjacent[1:]) and
            all(sign*v>0 for v in high_pairs) and
            all(abs(v)<ONE_THIRD for v in changes[1]["CI"]))
    meets = max(widths)<=ONE_THIRD
    result["largest_primary_half_width"] = max(widths)
    result["interpretation"] = (
        "shrinking" if meets and all(shrinking) and len(set(shrink_direction))==1 else
        "persistent" if meets and all(persistent) and len(set(persistence_direction))==1 else
        "inconclusive")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
