"""Exclude a three-step positive prefix ending on the zero face."""
import hashlib
import json
from pathlib import Path

import z3
from seven_clock_audit import cvc_status, prefix, require


def main():
    base, states, constraints = prefix(3, positive_final=False)
    solver = z3.SolverFor('QF_LIA')
    solver.set(timeout=180000)
    solver.add(*constraints, base >= 18, states[-1][-1] == 0)
    expression = solver.to_smt2()
    zs = str(solver.check())
    cs, _ = cvc_status(expression)
    require(zs == cs == 'unsat', (zs, cs))
    out = Path(__file__).with_name('seven_first_access_results.json')
    out.write_text(json.dumps(dict(z3=zs, cvc5=cs, smtlib=expression,
                                   sha256=hashlib.sha256(expression.encode()).hexdigest()), indent=2))
    print("Three positive sources followed by a zero, B >= 18:", zs, cs)


if __name__ == '__main__':
    main()
