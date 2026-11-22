"""Integer theorem-lift controls at three nearby-base residues and repeat counts."""
import json
import sys
from fractions import Fraction
from pathlib import Path


def main(certificate):
    root_slots = (1,2,4,8,16,3,6,12,11,13,9,17)
    pi = [Fraction(x) for x in certificate["lift_profile"]]
    cut = Fraction(certificate["lift_cut"])
    D = certificate["lift_D"]
    for i,B in enumerate((2*D,2*D+1,3*D+17)):
        L = B//D
        label_root = [0] * len(root_slots)
        entry = certificate["terminal_entry_from_release"]
        export = certificate["terminal_export_to_entry"]
        for j, stage in enumerate(certificate["stages"]):
            count=stage["coefficient_floor"]*L + ((i*j+j+i)%L)
            for t in range(len(label_root)):
                label_root[t] += stage["shift"][t] + count*stage["ray"][t]
            entry += stage["prefix"] + (count+1)*stage["period"]
            export += stage["prefix_export"] + (count+1)*stage["word_export"]
        g = [int(B*x)+(i+j)%L for j,x in enumerate(pi)]
        g[0] = (B*pi[0]-Fraction(sum(label_root),len(root_slots)+1)) // 1 + i%L
        for j,x in zip(root_slots,label_root[::-1]):
            g[j] = g[0] + x
        g[15] += B-sum(g)
        print(B,L,entry,export,cut.numerator,cut.denominator,
              L*certificate["spread_coefficient"],L*certificate["export_coefficient"],
              L*certificate["time_through_a_period_coefficient"],*g)


if __name__ == "__main__":
    main(json.loads(Path(sys.argv[1]).read_text()))
