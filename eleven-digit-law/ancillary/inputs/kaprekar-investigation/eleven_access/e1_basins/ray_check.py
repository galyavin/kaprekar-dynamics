"""Rebuild the singular-ray expansion query in exact linear arithmetic."""
import argparse
from subprocess import run


def affine(a, b):
    return f"(+ {a} {b})"


def negate(a):
    return f"(- {a})"


def disjunction(expressions):
    return "(or " + " ".join(expressions) + ")"


def query(steps, bound):
    names = lambda t: tuple(f"{x}{t}" for x in "vcdz")
    commands = ["(set-logic QF_LRA)"]
    for t in range(steps+1):
        commands.extend(f"(declare-fun {x} () Real)" for x in names(t))
    def require(expression):
        commands.append(f"(assert {expression})")

    require("(>= c0 d0)")
    for x in names(0):
        for sign in (x,negate(x)):
            require(f"(<= {sign} 1)")
    require(disjunction(f"(= {sign} 1)" for x in names(0) for sign in (x,negate(x))))

    for t in range(steps):
        v,c,d,z = names(t)
        v1,c1,d1,z1 = names(t+1)
        other = negate(affine(c,d))
        require(f"(= {v1} {negate(affine(d,z))})")
        require(f"(= {z1} {negate(affine(v,c))})")
        require(disjunction((f"(and (>= {v} {other}) (= {c1} {v}) (= {d1} {other}))",
                             f"(and (<= {v} {other}) (= {c1} {other}) (= {d1} {v}))")))
    for x in names(steps):
        for sign in (x,negate(x)):
            require(f"(<= {sign} {bound})")
    return "\n".join((*commands,"(check-sat)","(exit)",""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", default="yices-smt2")
    args = parser.parse_args()
    for steps,bound,expected in ((1,2,"sat"),(10,2,"unsat")):
        result = run([args.solver], input=query(steps,bound), capture_output=True,
                     text=True, check=True).stdout.strip()
        assert result == expected, (steps,result)
        print(f"Unit input, output sup norm <= {bound} after {steps} step(s): {result}")


if __name__ == "__main__":
    main()
