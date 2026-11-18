"""Build the fresh-base comparison table from pairgraph_fresh.log,
pairgraph_10009.log, mc_fresh.log and the exact limit values."""
import math, re
from fractions import Fraction as Fr
F = {0.04: Fr(713903, 39062500), 0.1: Fr(51653, 400000), 0.2: Fr(31645397267, 54000000000),
     0.3: Fr(430102558494119377593841787, 444140127977914368000000000)}
txt = open("pairgraph_fresh.log").read() + "\n" + open("pairgraph_10009.log").read()
blocks = [b for b in txt.split("\nB ") if b.strip()]
ex = {}
for b in blocks:
    b = b if b.startswith("B ") else "B " + b
    B = int(re.match(r"B (\d+)", b).group(1))
    d = {}
    for t in (0.04, 0.1, 0.2, 0.3):
        m = re.search(r"cdf t=%.6f P\(T<=tB\)=([0-9.]+)" % t, b); d[t] = float(m.group(1))
    d["mean"] = float(re.search(r"E\[T/B\]=([0-9.]+)", b).group(1))
    d["max"] = int(re.search(r"maxT (\d+)", b).group(1))
    d["wit"] = int(re.search(r"witness37 n \d+ pair \(\d+,\d+\) T (\d+)", b).group(1))
    d["cyc"] = re.findall(r"label (\S+) len \d+ .*? mass ([0-9.]+)", b)
    d["cond"] = re.findall(r"cdf t=0.200000 .*", b)[0]
    ex[B] = d
mc = {}
for b in open("mc_fresh.log").read().split("\n\n"):
    m = re.match(r"B (\d+)", b.strip())
    if not m: continue
    B = int(m.group(1)); d = {}
    for t in (0.04, 0.1, 0.2, 0.3):
        mm = re.search(r"cdf t=%.4f p=([0-9.]+) se=([0-9.]+)" % t, b); d[t] = (float(mm.group(1)), float(mm.group(2)))
    d["max"] = int(re.search(r"max T observed (\d+)", b).group(1))
    d["cyc"] = re.findall(r"cycle period (\d+) .* frac ([0-9.]+) se ([0-9.]+)", b)
    d["mean"] = re.search(r"mean T/B ([0-9.]+) se ([0-9.]+)", b).groups()
    mc[B] = d
print("| B | B mod 3 | t | F(t) limit | exact P(T<=tB) | exact - F | B*(exact-F) | MC (1e6) +- se | (MC-exact)/se |")
print("|---|---|---|---|---|---|---|---|---|")
for B in sorted(ex):
    for t in (0.04, 0.1, 0.2, 0.3):
        e = ex[B][t]; f = float(F[t])
        mcs = ""
        if B in mc:
            p, se = mc[B][t]; mcs = f"{p:.6f} +- {se:.6f} | {(p-e)/se:+.2f}"
        else:
            mcs = "- | -"
        print(f"| {B} | {B%3} | {t} | {f:.7f} | {e:.7f} | {e-f:+.2e} | {B*(e-f):+.2f} | {mcs} |")
print()
print("| B | exact max T | witness (37) T | B/3+sqrt(B)/3 | max-(B/3+sqrtB/3) | paper bounds [-5,+60] | MC max (1e6) | E[T/B] exact |")
print("|---|---|---|---|---|---|---|---|")
for B in sorted(ex):
    c = B/3 + math.sqrt(B)/3
    print(f"| {B} | {ex[B]['max']} | {ex[B]['wit']} | {c:.2f} | {ex[B]['max']-c:+.2f} | [{c-5:.2f}, {c+60:.2f}] ok={c-5 <= ex[B]['wit'] and ex[B]['max'] <= c+60} | {mc[B]['max'] if B in mc else '-'} | {ex[B]['mean']:.6f} |")
print()
for B in sorted(ex):
    print(B, "cycles exact:", ex[B]["cyc"], "| MC:", mc.get(B, {}).get("cyc"))
    print("   ", ex[B]["cond"])
