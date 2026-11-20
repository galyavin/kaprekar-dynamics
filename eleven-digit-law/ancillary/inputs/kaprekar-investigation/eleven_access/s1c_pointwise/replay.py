"""Build in a private snapshot and verify lane means and local sign cells."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
from itertools import product
import fcntl
import json
from math import lcm
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ACCESS = HERE.parent
ROOT = ACCESS.parents[1]
sys.path[:0] = [str(ACCESS / name) for name in ("z1_stages", "s1_selection", "r2_routes")]
from cover import record_stage_probes  # noqa: E402
from exact import density, in_pi, integrate, add, multiply, scale, constant, coordinate  # noqa: E402
from larger_cell import ANCHOR, GRID, integration  # noqa: E402

ARMS = (3, 5, 9, 11)
SMALL = tuple(Q(x, 50000) for x in (2938, 3225, 2886, 1913, 387))
LARGE = tuple(Q(x, GRID) for x in ANCHOR)
BOX = tuple((Q(a,100000), Q(b,100000))
            for a,b in ((94040,94130),(73770,73850),(57940,58015),(20,80),(200,300)))
BOX_SOURCE = [coordinate(i) for i in range(3)] + [
    add(constant(1),scale(add(coordinate(2),coordinate(4)),-1)),
    add(coordinate(1),scale(add(constant(1),coordinate(3)),Q(-1,2))),
]
MASK = (1 << 64) - 1
UPPER_CLIP = Q(1)


def box_point(point):
    a,b,c,h,g = point
    return a,b,c,1-c-g,b-(1+h)/2


def q_from_x(x):
    return tuple((v + 5 - i) / 6 for i,v in enumerate(x))


def x_from_q(q):
    return tuple(v * 6 - (5 - i) for i,v in enumerate(q))


def build(snapshot):
    for source in (ACCESS/"stage_graph/graph.cpp",ACCESS/"stage_graph/lottery.cpp",
                   HERE/"pointwise.cpp"):
        shutil.copy2(source,snapshot)
    record_stage_probes(snapshot/"graph.cpp", paired_fast=True)
    if shutil.which("pkg-config"):
        flags = subprocess.check_output(["pkg-config","--cflags","--libs","gmpxx"],text=True).split()
    else:
        prefix = Path(os.getenv("S1C_GMP_PREFIX","/opt/homebrew"))
        flags = ["-I"+str(prefix/"include"),"-L"+str(prefix/"lib"),"-lgmpxx","-lgmp"]
    subprocess.run(["c++","-std=c++17","-O3","-Wno-return-type",
                    "-Wno-deprecated-literal-operator","pointwise.cpp","-o","pointwise",*flags],
                   cwd=snapshot,check=True)
    return snapshot/"pointwise"


def run(executable, base, q, dump_guards=False):
    den = lcm(*(x.denominator for x in q))
    data = " ".join(map(str,(base,den,*(int(x*den) for x in q))))+"\n"
    cmd = [str(executable),*(["--guards"] if dump_guards else [])]
    lines = subprocess.check_output(cmd,input=data,text=True,timeout=180).splitlines()
    fields = lines[0].split()
    if fields[0] != "clock":
        return {"status": lines[0]}
    b = fields.index("states")
    c = fields.index("branches")
    return dict(value=Q(fields[1]), row=tuple(map(Q,fields[2:8])),
                guard_count=int(fields[9]), radius=Q(fields[11]),
                upper=Q(fields[12]), slope=Q(fields[13]), rays=int(fields[15]),
                canonical_reads=int(fields[17]), shorts=int(fields[19]), unpaired=bool(int(fields[21])),
                base_winner=bool(int(fields[22])),
                states=list(map(int,fields[b+1:c])), branches=list(map(int,fields[c+1:])),
                guards=[tuple(map(int,line.split()[1:])) for line in lines[1:]])


def affine_value(row, x):
    return sum(a*v for a,v in zip(row[:5],x))+row[5]


def cell_is_covered(read, vertices):
    for vertex in vertices:
        if any(affine_value(g,vertex)<=0 for g in read["guards"]):
            return False
    return True


def controls(evaluate_at):
    q = {"small": q_from_x(SMALL), "large": q_from_x(LARGE),
         "box": box_point(tuple((lo+hi)/2 for lo,hi in BOX))}
    laws = {name:evaluate_at(point) for name,point in q.items()}
    for name, panel in laws.items():
        for base, read in panel.items():
            if "status" in read:
                raise RuntimeError((name,base,read))
            assert read["upper"] < 1
    small = json.loads((ACCESS/"s1_selection/exact.json").read_text())["cone"]
    raw_large = json.loads((ACCESS/"s1_selection/larger-cell.json").read_text())
    for base in ARMS:
        speed = Q(small["mean_speeds"][base%6])
        assert laws["small"][base]["row"] == (Q(1,6)-Q(4,5)*speed, 2*speed,-2*speed,speed,0,0)
        expected = tuple(sum(Q(leaf["row"][i])*Q(leaf["weight"])
                             for leaf in raw_large["laws"][str(base % 6)]["leaves"])
                         for i in range(6))
        assert laws["large"][base]["row"] == expected
    small_integral = (laws["small"][3]["row"][3]-laws["small"][5]["row"][3])*Q(small["distance_integral"])
    assert small_integral == Q(small["original_3_minus_5"])
    large_arms = {}
    for base in (3,5):
        archive = raw_large["laws"][str(base)]
        stages = {int(i):list(map(Q,row)) for i,row in archive["stages"].items() if i!="5"}
        stages[5] = list(laws["large"][base]["row"])
        for i in range(6):
            stages[5][i] -= sum(row[i] for stage,row in stages.items() if stage!=5)
        large_arms[base] = dict(stages=stages,
                               leaves=[dict(row=list(map(Q,l["row"])), weight=Q(l["weight"]))
                                       for l in archive["leaves"]])
    large_integral = integration(large_arms)["contribution"]
    assert large_integral == Q(raw_large["cone"]["contribution"])

    box_q = BOX_SOURCE
    # Substitute the original box variables into the signed centered row.
    row = tuple(a-b for a,b in zip(laws["box"][5]["row"],laws["box"][11]["row"]))
    signed = constant(row[5])
    for i, qi in enumerate(box_q):
        signed = add(signed,scale(add(scale(qi,6),constant(i-5)),row[i]))
    box_integral = integrate(scale(multiply(density(box_q),signed),Q(1,2)),BOX)
    expected = json.loads((ACCESS/"r2_routes/control-weights.json").read_text())["clock_marginal_class5_minus11"]
    assert box_integral == Q(expected["original_contribution"]), (box_integral,expected)
    vertices = [x_from_q(box_point(v)) for v in product(*BOX)]
    return dict(laws=laws,integrals={"small":small_integral,"large":large_integral,"box":box_integral},
                box_common_guard_cell={base:cell_is_covered(read,vertices)
                                       for base,read in laws["box"].items()})


def splitmix(s):
    s = (s+0x9e3779b97f4a7c15) & MASK
    z = s
    z = ((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK
    z = ((z^(z>>27))*0x94d049bb133111eb)&MASK
    return s,z^(z>>31)


def original_probe(index):
    seed = indexed_seed(14078212001,index)
    digits = []
    for _ in range(11):
        seed, z = splitmix(seed)
        digits.append(Q(z>>26,1<<38))
    digits.sort()
    q = tuple(digits[-i-1]-digits[i] for i in range(5))
    for depth in range(501):
        if in_pi(q):
            return depth,q
        w = sorted((Q(1),*q,*(1-x for x in q)))
        q = tuple(w[-i-1]-w[i] for i in range(5))
    raise RuntimeError(f"no strict entrance for exploration row {index}")


def indexed_seed(seed, index):
    return splitmix((seed + index * 0x9e3779b97f4a7c15) & MASK)[1]


def strict_ball(q, reads, contrast):
    radius = min(r["radius"] for r in reads.values())
    for r in reads.values():
        assert r["upper"] < UPPER_CLIP
        if r["slope"]:
            radius = min(radius,(UPPER_CLIP-r["upper"])/r["slope"])
    x = x_from_q(q)
    row = tuple(a-b for a,b in zip(*(reads[b]["row"] for b in contrast)))
    value = affine_value(row,x)
    norm = sum(map(abs,row[:5]))
    if norm and value:
        radius = min(radius,abs(value)/norm)
    radius /= 2
    while not all(in_pi(q_from_x(tuple(z+t*radius for z,t in zip(x,s))))
                  for s in product((-1,1),repeat=5)):
        radius /= 2
    bounds = [(z-radius,z+radius) for z in x]
    source = [scale(add(constant(5-i),coordinate(i)),Q(1,6)) for i in range(5)]
    mass = integrate(scale(density(source),Q(1,6**5)), bounds)
    signed = constant(row[5])
    for i,a in enumerate(row[:5]):
        signed = add(signed,scale(coordinate(i),a))
    return dict(center=x,radius=radius,value=value,bounds=(value-radius*norm,value+radius*norm),
                direct_mass=mass,direct_contribution=integrate(
                    scale(multiply(density(source),signed),Q(1,6**5)),bounds))


def screen(evaluate_at, kind):
    if kind == "adversarial":
        rows = [(f"probe-{i}",q) for i in (149,88,161,131,8,6,0,2,5,9,1,3)
                for _,q in [original_probe(i)]]
    elif kind == "marked":
        marks = json.loads((ACCESS/"box_structure/marked-clock.json").read_text())
        rows = [(f"mark-{r['base']}-{r['r2']}",tuple(Q(x,r["base"]) for x in r["source"])) for r in marks]
    else:
        raise ValueError(kind)
    output = {}
    for name,q in rows:
        arms = evaluate_at(q)
        contrasts = {}
        if all("status" not in r for r in arms.values()):
            for pair in ((3,5),(9,11),(3,11),(9,5),(3,9),(5,11)):
                a,b = pair
                contrasts[f"{a}-{b}"] = dict(value=arms[a]["value"]-arms[b]["value"])
                if arms[a]["value"] < arms[b]["value"]:
                    contrasts[f"{a}-{b}"].update(strict_ball(q,arms,pair))
        output[name] = dict(q=q,arms=arms,contrasts=contrasts)
        print(name,{k:float(v["value"]) for k,v in contrasts.items()},file=sys.stderr,flush=True)
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument("scope",choices=("controls","adversarial","marked"))
    p.add_argument("--lock",action="store_true")
    p.add_argument("--jobs",type=int,default=4)
    args = p.parse_args()
    with (ROOT/".task-monitor/compute.lock").open() as lock:
        if args.lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix="task146-pointwise-") as tmp:
            executable = build(Path(tmp))
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                def evaluate_at(point):
                    values = pool.map(lambda b:run(executable,b,point,True), ARMS)
                    return dict(zip(ARMS,values))
                result = controls(evaluate_at) if args.scope=="controls" else screen(evaluate_at,args.scope)
    print(json.dumps(result,default=str,indent=2))


if __name__ == "__main__":
    main()
