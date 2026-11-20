"""Enumerate the same exact guarded cell in each physical class."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def run(snapshot, beta):
    cmd = ["./lottery", str(48000+beta), "300000", "252938", "203225",
           "152886", "101913", "50387"]
    result = subprocess.check_output(cmd, cwd=snapshot, text=True)
    out = {"clock": [], "edge": [], "gap": []}
    for row in result.splitlines():
        f = row.split()
        if f[0]=="base":
            out["states"]=int(f[3])
            out["inequalities"]=int(f[5])
            out["strict_radius"]=f[7]
        elif f[0]=="clock":
            out["clock"].append({"numerator":list(map(int,f[1:7])),
                                 "denominator":int(f[7]),"weight":f[8]})
        elif f[0]=="gap":
            out["gap"].append([int(f[1]),f[2]])
        elif f[0]=="edge":
            out["edge"].append({"stage":int(f[1]),"period":int(f[2]),
                                "ranks":[int(f[3]),int(f[4])],"failures":int(f[5]),
                                "decrement":int(f[6]),"residue_index":int(f[7]),
                                "drift":list(map(int,f[8:13])),
                                "distance_numerator":list(map(int,f[13:19])),
                                "distance_denominator":int(f[19]),"weight":f[20]})
        else:
            raise ValueError(row)
    return out


def main():
    source=Path(__file__).resolve().parent
    repo=source.parents[2]
    with tempfile.TemporaryDirectory(prefix="task097-cell-") as directory:
        snapshot=Path(directory)
        for name in ("graph.cpp","lottery.cpp"):
            shutil.copyfile(source/name,snapshot/name)
        subprocess.run(["clang++","-O3","-std=c++17","-Wall","-Wextra",
                        "-Wno-return-type","lottery.cpp","-o","lottery"],
                       cwd=snapshot,check=True)
        with (repo/".task-monitor"/"compute.lock").open("r+") as lock:
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
            with ThreadPoolExecutor(max_workers=4) as pool:
                tables=list(pool.map(lambda b:run(snapshot,b),range(6)))
        print(json.dumps({"center_grid":300000, "center_q":[252938,203225,152886,101913,50387],
                          "mod12_to_table":{r:r%6 for r in range(12)},"tables":tables},indent=2))


if __name__=="__main__":
    main()

