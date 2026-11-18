"""Robustness replay: re-run the entry and clock z3 encodings with the logic-specific
solvers (QF_LRA / QF_LIA) and two non-default random seeds."""
import z3, subprocess, sys, re
orig = z3.Solver
def make(logic):
    def S(*a, **k):
        s = z3.SolverFor(logic)
        s.set("random_seed", SEED)
        return s
    return S
for SEED in (7, 12345):
    for mod, logic in (("z_entry", "QF_LRA"), ("z_clock", "QF_LIA")):
        z3.Solver = make(logic)
        z3.set_param("smt.random_seed", SEED)
        print(f"--- {mod} with SolverFor({logic}), seed {SEED}", flush=True)
        src = open(mod + ".py").read()
        g = {"__name__": "__main__"}
        sys.argv = [mod]
        exec(compile(src, mod + ".py", "exec"), g)
z3.Solver = orig
