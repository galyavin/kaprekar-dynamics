"""Summarize sim.c logs: per base, quantiles of sigma/B, T/B, (T-sigma)/B, max B*|Z-v/6|,
and counts of any failed structural checks.  Usage: python3 summarize.py LOGDIR"""
import glob
import os
import re
import sys

pat = re.compile(r"B=(\d+) s=(\d+) sigma=(-?\d+) T=(\d+) mu=(\d+) lam=(\d+) dist=([\d.]+) Rsig=(-?\d+) "
                 r"maxRafter=(-?\d+) maxjump=(-?\d+) maxhit=(\d+) bad\[pair=(\d+) jump=(\d+) hit=(\d+) env=(\d+) cyc=(\d+) consist=(\d)\]")


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


def main():
    rows = []
    tot = dict(n=0, capped=0, never=0, bad=0, inconsist=0)
    gmax = dict(sig=0, T=0, Tms=0, dist=0, Rafter=0, jump=-99, hit=0, lam=0)
    for fn in sorted(glob.glob(os.path.join(sys.argv[1], "sim_B*.log")), key=lambda f: int(re.findall(r"\d+", os.path.basename(f))[0])):
        sig, T, Tms, dist, lams = [], [], [], [], {}
        B = None
        capped = never = bad = 0
        Rafter = jump = hit = -99
        for line in open(fn):
            if "CAPPED" in line:
                capped += 1
                continue
            m = pat.search(line)
            if not m:
                continue
            g = m.groups()
            B = int(g[0])
            s, t, mu, lam = int(g[2]), int(g[3]), int(g[4]), int(g[5])
            if s < 0:
                never += 1
            else:
                sig.append(s / B)
                Tms.append((t - s) / B)
            T.append(t / B)
            dist.append(float(g[6]))
            lams[lam] = lams.get(lam, 0) + 1
            Rafter = max(Rafter, int(g[8])); jump = max(jump, int(g[9])); hit = max(hit, int(g[10]))
            if any(int(x) for x in g[11:16]):
                bad += 1
            if g[16] != "1":
                tot["inconsist"] += 1
        n = len(T)
        tot["n"] += n; tot["capped"] += capped; tot["never"] += never; tot["bad"] += bad
        gmax["sig"] = max(gmax["sig"], max(sig)); gmax["T"] = max(gmax["T"], max(T)); gmax["Tms"] = max(gmax["Tms"], max(Tms))
        gmax["dist"] = max(gmax["dist"], max(dist)); gmax["Rafter"] = max(gmax["Rafter"], Rafter)
        gmax["jump"] = max(gmax["jump"], jump); gmax["hit"] = max(gmax["hit"], hit); gmax["lam"] = max(gmax["lam"], max(lams))
        rows.append(f"| {B} | {B % 10} | {B % 5} | {B % 2} | {n} | {q(sig, .5):.4f} | {q(sig, .99):.4f} | {max(sig):.4f} | "
                    f"{q(T, .5):.4f} | {q(T, .99):.4f} | {max(T):.4f} | {max(Tms):.4f} | {max(dist):.2f} | {never} | {capped} | {bad} | "
                    f"{sorted(lams)} |")
    hdr = ("| B | B mod 10 | mod 5 | parity | n | sigma/B median | sigma/B q99 | sigma/B max | T/B median | T/B q99 | T/B max | "
           "max (T-sigma)/B | max B*|Z-v/6| | sigma=inf | capped | failed checks | periods seen |\n" + "|---" * 17 + "|")
    print(hdr)
    print("\n".join(rows))
    print()
    print(f"TOTAL samples={tot['n']} sigma=inf={tot['never']} capped={tot['capped']} samples with any failed structural check={tot['bad']} "
          f"word/difference inconsistencies={tot['inconsist']}")
    print(f"GLOBAL max sigma/B={gmax['sig']:.4f} (cutoff 4509), max T/B={gmax['T']:.4f}, max (T-sigma)/B={gmax['Tms']:.4f}, "
          f"max B*|Z-v/6|={gmax['dist']:.3f} (bound 420), max R after first R<=120: {gmax['Rafter']} (bound 504), "
          f"max one-step R increase after S_B: {gmax['jump']} (bound 12; claimed 6), max lower-hit delay: {gmax['hit']} (bound 32), "
          f"max period: {gmax['lam']}")


if __name__ == "__main__":
    main()
