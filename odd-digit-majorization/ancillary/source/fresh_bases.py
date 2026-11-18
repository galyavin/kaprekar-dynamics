"""Run an independent direct subtraction test on the five displayed kernels."""
from fractions import Fraction
from itertools import product
import subprocess
import tempfile
from pathlib import Path

from nine_core_checks import route_domains
from nine_entry import evaluate

# Certify the sampling boxes independently of every measured outcome.
# All four original differences exceed one half on this inverse branch.
domain = next(item["domain"] for item in route_domains()
              if item["route"] == ("A0", "A:B1", "AB"))
for base in range(300035, 300040):
    for signs in product((-210, 210), repeat=4):
        digits = [base * n // 50 + s for n, s in zip((43, 39, 31, 28), signs)]
        radius = [2 * Fraction(x, base) - 1 for x in digits]
        if min(evaluate(row, (*radius, 1)) for row in domain) <= 0:
            raise ArithmeticError("Sampling box crosses a route boundary")

with tempfile.TemporaryDirectory() as temp:
    binary = Path(temp) / "separating-route"
    source = Path(__file__).with_name("separating_route.cpp")
    subprocess.run(["c++", "-std=c++17", "-O3", str(source), "-o", str(binary)], check=True)
    subprocess.run([str(binary)], check=True)
