"""Post-freeze comparison with the pass's checkpoint sections."""
import json
from pathlib import Path

import orbit_matching as archived
from reconstruct import build_sections

records = json.loads(Path("endpoint_verification.json").read_bytes())["residues"]
independent = build_sections(records)
assert independent.keys() == archived.SECTIONS.keys()
origins = 0
for label, section in independent.items():
    entry = archived.SECTIONS[label]
    assert section.word == entry["word"], (label, section.word, entry["word"])
    assert section.pin == entry["u"], (label, section.pin, entry["u"])
    for K, (columns, origin) in entry["lattice"].items():
        assert section.columns == columns, (label, K, section.columns, columns)
        assert section.origin(K) == origin, (label, K, section.origin(K), origin)
        origins += 1
    cols, _, origin = archived.section(label, 0)
    delta = archived.period_displacement(section.word, cols, origin)
    assert section.displacement == delta, (label, section.displacement, delta)

print(f"Exact agreement for {len(independent)} words, phase displacements and pins")
print(f"Exact agreement for {origins} ordered columns and origins across K=0..11")
