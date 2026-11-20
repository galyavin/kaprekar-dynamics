"""Carry each marked physical box witness through Task 097's exact lane choices."""
from fractions import Fraction as Q
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from control_weights import box_profiles

HERE = Path(__file__).resolve().parent

def matched_profile(row, leading, table):
    base = row["base"]
    a,b,c,d,e = row["source"]
    a,b,c,h,g = (Q(u,base) for u in (a,b,c,2*(b-e)-base,base-c-d))
    n = (a+c-Q(3,2)-Q(11,2)*h-5*g)/22
    p,r = c-h/2-2*n,b-Q(3,2)*h-g-4*n
    profiles = table[base%72][str(row["r2"])]
    matches = []
    for record in profiles:
        s0,sp,sr = map(Q,record["coefficients"])
        if leading==a+s0+sp*p+sr*r:
            matches.append(record["coefficients"][0])
    assert len(matches)==1
    assert row["F"]==int(matches[0] in ("139/144","863/918"))
    return matches[0]


def main():
    with tempfile.TemporaryDirectory(prefix="r2-box-route-") as snap:
        for source in (HERE.parent/"stage_graph/graph.cpp", HERE/"diagnostic.cpp", HERE/"box_route.cpp"):
            shutil.copy2(source,snap)
        subprocess.run(["clang++","-std=c++17","-O3","-Wno-return-type",
                        "box_route.cpp","-o","box_route"],cwd=snap,check=True,timeout=40)
        rows = json.loads((HERE.parent/"box_structure/marked-clock.json").read_text())
        data = "".join(" ".join(map(str,(r["base"],r["r2"],r["word_clock"],r["N"],*r["source"])))+"\n"
                       for r in rows)
        result = subprocess.check_output(["./box_route"],cwd=snap,input=data,text=True,timeout=30)
    assert len(result.splitlines())==len(rows)
    table = box_profiles()
    print("base pin actual_time route_time leading_clock stages short_steps failures "
          "macro_failures certified_N X1 X2 X3 X4 X5 unpaired_repeat_length:count profile_intercept")
    for row,line in zip(rows,result.splitlines()):
        f = line.split()
        assert (int(f[0]),int(f[1]))==(row["base"],row["r2"])
        print(line,matched_profile(row,Q(f[4]),table))


if __name__ == "__main__":
    main()
