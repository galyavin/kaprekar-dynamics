"""Run only the five-digit part of the unchanged arithmetic scaling check."""
import json
from fractions import Fraction
from pathlib import Path
from arithmetic_scaling_checks import five_checks
result = five_checks()
assert all(Fraction(row['max_actual_center_error']) <= 8 for row in result['clock_checks'])
Path('p2_overshoot_results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
