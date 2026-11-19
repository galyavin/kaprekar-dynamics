"""Runs in a private copy of the box certificate.  Re-implements the main loop of verify.py but records
every arithmetic leaf (lattice, lane, leading profile).  Writes JSON lines to stdout."""
import json, sys
import verify as V

def main():
    for label, guards, section, entry, clock in V.initial_sections():
        _, _, _, r2 = label
        for k in (0, 2, 4):
            for arm in (0, 1):
                base = 6*(k+arm) + 5
                p0 = ((r2 - section.clock_constant)//2 - base) % 36
                for residue in range(6):
                    guards.current = (label, k, arm, residue)
                    start = V.Lattice((base, p0 + 2*residue, residue), ((36, 0, 0), (0, 36, 0), (0, 12, 6)))
                    for end, jumps, tail, lattice in V.branches(entry, clock, section, start, guards):
                        leading = end.leading()
                        lane = V.profile(leading)[0]
                        print(json.dumps({'label': label, 'k': k, 'arm': arm, 'rho': residue, 'base': base,
                                          'point': lattice.point, 'columns': lattice.columns, 'lane': lane,
                                          'profile': leading.record()[:3], 'verifier_phases': lattice.base_phases(base)},
                                         separators=(',', ':')))
main()
