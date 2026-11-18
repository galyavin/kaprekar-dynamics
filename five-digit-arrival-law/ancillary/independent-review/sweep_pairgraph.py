"""Run pairgraph (with finite-base reduction checks) on every base in a range
and summarise.  usage: sweep_pairgraph.py BIN lo hi"""
import math, re, subprocess, sys
binp, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
bad = []
worst = {"lo9": (1e9, None), "hi9": (-1e9, None), "six": 0, "glob_slack": 1e9,
         "c23lo": 99, "c23hi": -99, "c36": -10**9}
labels_bad = []
sec = []
for B in range(lo, hi + 1):
    out = subprocess.run([binp, str(B), "checks"], capture_output=True, text=True).stdout
    ok = True
    for key in ("ok3 1", "ok9 1", "check23", "checksix", "checkglobal", "check36"):
        pass
    for line in out.splitlines():
        if line.startswith("check") and " ok 0" in line: ok = False
        if line.startswith("maxT") and ("ok3 0" in line or "ok9 0" in line): ok = False
        if line.startswith("rare") and "match 0" in line: ok = False
        if line.startswith("witness37") and " ok 0" in line: ok = False
        if line.startswith("cycle") and "label OTHER" in line: ok = False
        if "never 1" in line: ok = False
    m = re.search(r"maxT (\d+)", out); mT = int(m.group(1))
    sec.append((B, mT - B / 3 - math.sqrt(B) / 3))
    m = re.search(r"3h-2M in \[(-?\d+),(-?\d+)\]", out)
    worst["c23lo"] = min(worst["c23lo"], int(m.group(1))); worst["c23hi"] = max(worst["c23hi"], int(m.group(2)))
    worst["six"] = max(worst["six"], int(re.search(r"J_B->H_B u T2 = (\d+)", out).group(1)))
    g = int(re.search(r"->H_B u T2 = (\d+)  ceil", out).group(1))
    worst["glob_slack"] = min(worst["glob_slack"], math.ceil(math.sqrt(B)) + 13 - g)
    worst["c36"] = max(worst["c36"], int(re.search(r"\(2u0-B\) = (-?\d+)", out).group(1)))
    # labels must match B mod 3 (+T2 iff even)
    labs = sorted(re.findall(r"label (\S+)", out))
    exp = {0: ["T1"], 1: ["T3", "T4"], 2: ["T5"]}[B % 3] + (["T2"] if B % 2 == 0 else [])
    if labs != sorted(exp): ok = False; labels_bad.append((B, labs))
    if not ok: bad.append(B)
print(f"bases {lo}..{hi}: failing bases {bad}; label mismatches {labels_bad}")
print(f"  3h-2M on H_B ranges over [{worst['c23lo']},{worst['c23hi']}] (claim [-9,18])")
print(f"  max J_B->H_B u T2 steps {worst['six']} (claim 6)")
print(f"  min slack of global entry vs ceil(sqrtB)+13: {worst['glob_slack']}")
print(f"  max 2max(a_j,b_j)-(2u0-B) in (36): {worst['c36']} (claim 40)")
mn = min(sec, key=lambda t: t[1]); mx = max(sec, key=lambda t: t[1])
print(f"  maxT - B/3 - sqrt(B)/3: min {mn[1]:.3f} at B={mn[0]}, max {mx[1]:.3f} at B={mx[0]}  (claim in [-5, 60])")
for B, d in sec:
    if B % 100 == 0 or B == lo: print(f"    B={B}: maxT-B/3-sqrtB/3 = {d:.3f}")
