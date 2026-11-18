"""Summarise the fresh-base campaign logs (logs/*.txt) against the predicted limits."""
import glob, math, re, os
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
grid = [0.01 * k for k in range(26)]

# ---- limit law ----
lim = {}
txt = open("logs/limitstats.txt").read()
print(txt.strip())
for beta in (0, 1):
    m = re.search(r"beta=%d P\(fast\)=(\S+) mean=(\S+) var=(\S+) E e\^\{4t\}=(\S+) E e\^\{-4t\}=(\S+)" % beta, txt)
    cdf = [float(t) for t in re.search(r"beta=%d CDF (.*)" % beta, txt).group(1).split()]
    lim[beta] = dict(pf=float(m.group(1)), mean=float(m.group(2)), var=float(m.group(3)), e4=float(m.group(4)), em4=float(m.group(5)), cdf=cdf)

print("\n== Two clocks / parity weights (boxes of 512 residue classes) ==")
print("%-9s %4s %5s | %-28s | %8s %8s %6s" % ("B", "B%8", "pts", "fast fraction (side: value x count)", "maxdev", "min gapB", "far"))
for f in sorted(glob.glob("logs/clocks_*.txt"), key=lambda s: int(re.findall(r"\d+", s)[-1])):
    B = int(re.findall(r"\d+", f)[-1])
    rows = [l for l in open(f) if l.startswith("xi")]
    stats = {}
    maxdev = 0
    mingap = 1e18
    far = 0
    bad = []
    for l in rows:
        side = int(re.search(r"side ([+-]\d)", l).group(1))
        frac = float(re.search(r"frac (\S+)", l).group(1))
        md = float(re.search(r"maxdev (\S+)", l).group(1))
        gap = float(re.search(r"gapB (\S+)", l).group(1))
        far += int(re.search(r"far (\d+)", l).group(1))
        maxdev = max(maxdev, md)
        mingap = min(mingap, gap)
        stats.setdefault((side, frac), 0)
        stats[(side, frac)] += 1
        exp = 0.375 if B % 2 else (0.5 if side > 0 else 0.75)
        if abs(frac - exp) > 1e-12:
            bad.append(l.strip())
    desc = " ".join("%+d:%.4g x%d" % (s, fr, n) for (s, fr), n in sorted(stats.items()))
    print("%-9d %4d %5d | %-28s | %8.1f %8.0f %6d" % (B, B % 8, len(rows), desc, maxdev, mingap, far))
    for l in bad:
        print("    off-prediction box:", l)

print("\n== Population law by uniform random words (word-level digit subtraction) ==")
print("%-7s %7s | %9s %9s | %9s %9s | %9s %9s | %8s | %7s" % ("B", "N", "mean", "limit", "var", "limit", "Ee^{4T/B}", "limit", "sup|dF|", "maxT/B"))
for f in sorted(glob.glob("logs/pop_*.txt"), key=lambda s: int(re.findall(r"\d+", s)[-1])):
    B = int(re.findall(r"\d+", f)[-1])
    T = np.loadtxt(f) / B
    beta = B % 2
    L = lim[beta]
    se = math.sqrt(T.var() / len(T))
    cdf = [float((T <= g).mean()) for g in grid]
    sup = max(abs(a - b) for a, b in zip(cdf, L["cdf"]))
    print("%-7d %7d | %9.5f %9.5f | %9.6f %9.6f | %9.5f %9.5f | %8.4f | %7.4f   (SE of mean %.5f)" % (
        B, len(T), T.mean(), L["mean"], T.var(), L["var"], np.exp(4 * T).mean(), L["e4"], sup, T.max(), se))

print("\n== Exact weighted law at fresh bases (all B^7-B words; difference-level DP) ==")
print("%-6s | %9s %9s | %9s %9s | %9s %9s | %8s" % ("B", "mean", "limit", "var", "limit", "Ee^{4T/B}", "limit", "sup|dF|"))
for f in sorted(glob.glob("logs/exhaust_*.txt")):
    for blk in re.finditer(r"B=(\d+) weights total=\S+ expected=\S+ rel=(\S+)\s+mean\(T/B\)=(\S+) var=(\S+)\s+E e\^\{4T/B\}=(\S+) E e\^\{-4T/B\}=(\S+)\nB=\d+ CDF (.*)", open(f).read()):
        B = int(blk.group(1))
        L = lim[B % 2]
        cdf = [float(t) for t in blk.group(7).split()]
        sup = max(abs(a - b) for a, b in zip(cdf, L["cdf"]))
        print("%-6d | %9.5f %9.5f | %9.6f %9.6f | %9.5f %9.5f | %8.4f   (weight total rel. err %s)" % (
            B, float(blk.group(3)), L["mean"], float(blk.group(4)), L["var"], float(blk.group(5)), L["e4"], sup, blk.group(2)))
print("\n== Exhaustive maxima / destination / lemma checks ==")
for f in sorted(glob.glob("logs/exhaust_*.txt")):
    for l in open(f):
        if "max h=" in l or "fails=" in l:
            pass
    lines = open(f).read().splitlines()
    for l in lines:
        if "max h=" in l:
            print(l)
