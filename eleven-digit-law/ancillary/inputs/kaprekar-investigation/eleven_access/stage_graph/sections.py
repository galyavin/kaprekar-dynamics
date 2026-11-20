"""Certify the first rank-four stage on every strict first-exit cell."""
from collections import defaultdict
from fractions import Fraction as Q
import json
from math import gcd, lcm
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

from pi_geometry import Poly, dirichlet_mass, within_interlevel, minus, relative_volume, dot


def leading_gap_row(xrow):
    a, b, den = xrow[:5],xrow[5],xrow[6]
    return tuple(Q(b+sum(c*(6*int(j>i)-(5-i)) for i,c in enumerate(a)),den)
                 for j in range(6))


def representative(poly, rng):
    w = [rng.randrange(1,997) for v in poly.vertices]
    s = sum(w)
    y = [sum(v[i]*t for v,t in zip(poly.vertices,w))/s for i in range(6)]
    den = lcm(*(v.denominator for v in y))
    return den, *(int(den*sum(y[i+1:])) for i in range(5))


def candidate_rows(rep, source, base):
    data = subprocess.check_output(["./sections"], input=" ".join(map(str,(base,*rep)))+"\n",
                                   text=True, cwd=source)
    lanes={}
    for row in data.splitlines():
        fields=row.split()
        kind=fields[0]
        if kind=="face":
            count=int(fields[2])
        elif kind=="remainder":
            residue=int(fields[1])
            lane={"guards":[]}
            lanes[residue]=lane
        elif kind=="stage":
            lane["stage"]=list(map(int,fields[1:]))
        elif kind=="guard":
            g=tuple(map(int,fields[1:]))
            lane["guards"].append(g)
        elif kind=="gap":
            return None
        else:
            raise ValueError(row)
    assert len(lanes)==count
    return lanes


def full_cell_law(poly, lanes, cache):
    # All strict guards of the bounded prefix and the candidate whole word
    # have the claimed sign on the entire relative first-winner interior.
    result={}
    for k,lane in sorted(lanes.items()):
        guards=sorted(set(lane["guards"]))
        falling=defaultdict(list)
        for g in guards:
            row=leading_gap_row(g)
            lower=min(dot(row,v) for v in poly.vertices)
            upper=max(dot(row,v) for v in poly.vertices)
            assert lower>=0, (lane["stage"],g,lower)
            if g[8]:
                assert upper>0, (lane["stage"],g,upper)
                win=tuple(a/g[8] for a in row)
                falling[win].append(g)
        profile=tuple(sorted(falling))
        if profile not in cache:
            table={}
            for a in profile:
                c=poly.restrict(minus(b,a) for b in profile if b!=a)
                vol=relative_volume(c)
                if vol: table[a]=vol
            assert sum(table.values())==relative_volume(poly)
            cache[profile]=table
        table=cache[profile]
        entries=[]
        for win,volume in table.items():
            contact=any(g[7] for g in falling[win])
            if contact:
                assert all(g[7] == 1 and g[11] in (1,2) for g in falling[win])
                assert any(g[11] == 2 for g in falling[win])
                if ("original",win,profile) not in cache:
                    c=poly.restrict(minus(b,win) for b in profile if b!=win)
                    cache["original",win,profile]=dirichlet_mass(c)
            else:
                assert all(g[7] == 0 for g in falling[win])
            # Select the least exact affine intercept on this physical sheet.
            # Exact ties give the same repeat count throughout the sheet.
            winner = min(falling[win],key=lambda g: Q(g[9],g[8]))
            classes=[]
            for g in (winner,):
                stride = g[10]
                classes.append({
                    "x_numerator":g[:6], "x_denominator":g[6],
                    "base_coefficient":g[7],"decrement":g[8],
                    "exit_remainders":list(range(g[9]%stride,g[8],stride)),
                    "each_lane":Q(stride,g[8]),
                    "compatible_base_shift":stride//gcd(g[7],stride)})
            entries.append({"gap_time_row":win,
                            "conditional_volume":volume/relative_volume(poly),
                            "leads_to_interlevel_contact":contact,
                            "original_contact_weight":cache["original",win,profile] if contact else Q(0),
                            "forms":classes})
        result[k]={"short_steps":lane["stage"][0], "period":lane["stage"][1],
                   "failures":lane["stage"][2], "drift":lane["stage"][3:],
                   "stage2_winners":entries}
    return result


def produce(snapshot):
    initial,inter=within_interlevel()
    pi=Poly.root().restrict([*sorted(inter),*sorted(initial)])
    pi_volume=relative_volume(pi)
    results=[]
    rng=random.Random(9704741)
    contact={b:Q(0) for b in range(12)}
    original_contact={b:Q(0) for b in range(12)}
    seen=set()
    for a,(phase,pair,_) in sorted(initial.items()):
        poly=pi.restrict(minus(b,a) for b in initial if b!=a)
        if poly is None: continue
        cache={}
        by_base={}
        frac=relative_volume(poly)/pi_volume
        for b in range(12):
            for attempt in range(12):
                rep=representative(poly,rng)
                sections=candidate_rows(rep,snapshot,48000+b)
                if sections is not None: break
            else:
                raise RuntimeError(f"no generic direction for {phase} {pair} {b}")
            law=full_cell_law(poly,sections,cache)
            for lane in law.values():
                seen.add((lane["short_steps"],lane["period"],lane["failures"]))
                contact[b]+=frac*sum(x["conditional_volume"] for x in lane["stage2_winners"]
                                       if x["leads_to_interlevel_contact"])/len(law)
                original_contact[b]+=sum(x["original_contact_weight"]
                                         for x in lane["stage2_winners"])/len(law)
            by_base[b]=law
        results.append({"first_phase":phase,"first_pair":pair,"first_normal":a,
                        "uniform_pi_weight":frac,"stage2_by_B_mod_12":by_base})
    assert len(results)==10 and sum(r["uniform_pi_weight"] for r in results)==1
    return {"stage2_shapes":sorted(seen),"contact_weight_by_B_mod_12":contact,
            "original_direct_contact_weight_by_B_mod_12":original_contact,"faces":results}


def main():
    here=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="task097-sections-") as temp:
        snapshot=Path(temp)
        for name in ("sections.cpp","lottery.cpp","graph.cpp"):
            shutil.copyfile(here/name,snapshot/name)
        subprocess.run(["clang++","-O3","-std=c++17","-Wno-return-type","sections.cpp",
                        "-o","sections"],cwd=snapshot,check=True)
        print(json.dumps(produce(snapshot),default=str,indent=None,separators=(",",":")))


if __name__=="__main__":
    main()

