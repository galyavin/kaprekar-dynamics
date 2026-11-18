"""Verify the canonical domains and every printed fast-residue column."""
import json
from pathlib import Path

from seven_return_checks import vertices
from seven_return_law import POSITIVE_TETRAHEDRA, NEGATIVE_TETRAHEDRA


def canonical_slacks(x, y, z, sign):
    a = 2*x-1
    b = 2*y-a
    c = 2*z-b
    return (1-a, a-b, b-c, c, 2*a-1, a+b-1, a-b-c, b+2*c-2*a+1,
            2-a-2*b-c, x+z-2*y, x-y-z, sign*(x-3*z))


def check_geometry(sign, tetrahedra):
    values = [canonical_slacks(*point, sign) for point in
              ((0,0,0), (1,0,0), (0,1,0), (0,0,1))]
    rows = [(x-d, y-d, z-d, d) for d,x,y,z in zip(*values)]
    assert vertices(rows) == {q for tetra in tetrahedra for q in tetra}


def main():
    check_geometry(1, POSITIVE_TETRAHEDRA)
    check_geometry(-1, NEGATIVE_TETRAHEDRA)
    root = Path(__file__).resolve().parent
    positive = json.loads((root/'seven_facet_results.json').read_text())['rate_tables']
    negative = json.loads((root/'seven_return_results.json').read_text())['negative_phase_tables']
    expected = [('12','03','03','03','03','12','12','12'),
                ('02','02','','02','02','02','','02'),
                ('0123','03','0123','12'), ('','023','','012')]
    tables = [t['rates_by_S_mod8_and_m_mod8'] for t in positive] + [t['rates'] for t in negative]
    for table, reference in zip(tables, expected):
        fast = tuple(''.join(str(i) for i,r in enumerate(row) if r == '3/2') for row in table)
        assert fast == reference
    print('Canonical polyhedra and all 96 residue entries agree.')


if __name__ == '__main__':
    main()
