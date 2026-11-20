"""Select theorem checks from the archived verifiers."""
import importlib.util
import json
from pathlib import Path
import sys

from fractions import Fraction

ROOT = Path(__file__).resolve().parent
EVIDENCE = ROOT / "kaprekar-investigation/cross_length"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, EVIDENCE / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def check_period_eight():
    seventeen = load("verify", "x2_seventeen/verify.py")
    assert tuple(map(seventeen.check_route, seventeen.ROUTES)) == (
        (4, 0, 177), (16, 4, 179), (52, 12, 155), (48, 18, 124), (52, 13, 72),
    )
    seventeen.tail()
    seventeen.physical_controls()
    fifteen = load("fifteen", "x2b_fifteen/verify.py")
    fifteen.premise()
    fifteen.fifteen()


def check_shifted_order_twelve():
    module = load("verify_root", "lc_classification/verify_root.py")
    source = ROOT / "docs/drafts/research-direction/astra-strategy-16/period132-root.jsonl"
    _, _, cycle, budget = module.verify(json.loads(source.read_text()))
    assert len(cycle) == 11 and budget == 809022


def check_thirteen_digit_cycle_identity():
    # Compare the two least cycles in the same unscaled centered coordinates.
    module = load("class0_replay", "survival_reconstruction/class0.py")
    route = json.loads((EVIDENCE / "tl3_survival/route.json").read_text())
    _, cycle1 = module.analytic_tail(tuple(route["tail"]), 1)
    outer_q = tuple(-sum(module.OUTER[:i]) for i in range(1, 7))
    _, cycle2 = module.analytic_tail(outer_q, 7)
    assert len(cycle1) == len(cycle2) == 177
    assert {tuple(7 * v for v in x) for x in cycle1} == set(cycle2)
    print("Both thirteen-digit certificates reach the same translated cycle.")


def check_full_orbit_example():
    from common import difference_rows, literal

    # All formulas are affine on the parity lane once the listed guard signs
    # hold. The endpoints here check the intercept and first step in a.
    for a in (7, 9, 1001):
        B = 71*a
        x = (-a, -2*a, -3*a, -4*a+12, a-18, -18,
             -a-6, -2*a+6, -3*a+18, 2*a-12, a)
        # The gap offsets are the differences of successive q offsets.
        g = tuple(Fraction(B + left-right, 12) for left, right in
                  zip((0, *x), (*x, 0)))
        assert all(v.denominator == 1 for v in g)
        g = tuple(int(v) for v in g)
        assert min(g) == (B-5*a+30)//12
        assert 13*min(g) > B+1
        assert difference_rows(g) == (0,-1,(-a+3)//2,(-a+3)//2,-2,-3)
        assert literal(g, 1) == g
    print("Literal fixed-point controls at length 23.")


if __name__ == "__main__":
    sys.path.insert(0, str(EVIDENCE / "survival_reconstruction"))
    check_period_eight()
    check_shifted_order_twelve()
    check_thirteen_digit_cycle_identity()
    check_full_orbit_example()
