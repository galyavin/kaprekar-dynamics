"""Replay the proof from a disposable copy of the encodings and certificates."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

STAGES = {
    "entry": ["nine_global_entry_audit", "nine_global_checks"],
    "routes": ["nine_entry_checks", "nine_phase_checks", "nine_a_facet_checks",
               "nine_b_facet_checks", "nine_aa_checks", "nine_ab_checks",
               "nine_bb_checks", "nine_b_routes_checks", "nine_bbb_checks",
               "nine_aaa_checks", "nine_aab_checks", "nine_releases_checks",
               "nine_core_checks"],
    "weights": ["nine_kernel_checks"],
    "covariance": ["odd_base_covariance_checks"],
    "zeros": ["odd_zero_access_audit", "nine_zero_return_audit",
              "nine_clock_reduction_checks"],
    "descent": ["nine_clock_audit", "nine_clock_checks"],
    "fresh": ["fresh_bases"],
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = list(STAGES)
    parser.add_argument("stage", choices=["all", *choices, "independent"], default="all", nargs="?")
    args = parser.parse_args()
    selection = choices if args.stage == "all" else ([] if args.stage == "independent" else [args.stage])
    # The arithmetic replay of the descent reads the regenerated finite graph
    # controls, so a stand-alone descent first runs the reduction.
    if args.stage == "descent":
        selection = ["zeros", "descent"]
    source = Path(__file__).resolve().parent / "source"
    with tempfile.TemporaryDirectory(prefix="kaprekar-nine-") as temp:
        work = Path(temp) / "source"
        shutil.copytree(source, work)
        for stage in selection:
            for script in STAGES[stage]:
                print(f"Running {script}", flush=True)
                subprocess.run([sys.executable, "-B", "-u", script + ".py"],
                               cwd=work, check=True,
                               env={**os.environ, "PYTHONNOUSERSITE": "1"})
    if args.stage in ("all", "independent"):
        subprocess.run([sys.executable, "-B", str(source.parent / "independent-review" / "replay.py")], check=True)
    print("All selected checks passed.", flush=True)

if __name__ == "__main__":
    main()
