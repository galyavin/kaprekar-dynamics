import sys,subprocess
n=int(sys.argv[1]); bound=sys.argv[2]
L=["(set-logic QF_LRA)"]
for i in range(4): L.append(f"(declare-const x0_{i} Real)")
for k in range(n):
    v,c,d,z=[f"x{k}_{i}" for i in range(4)]
    for i in range(4): L.append(f"(declare-const x{k+1}_{i} Real)")
    L.append(f"(assert (= x{k+1}_0 (- (- {d}) {z})))")
    L.append(f"(assert (= x{k+1}_1 (ite (>= {v} (- (- {c}) {d})) {v} (- (- {c}) {d}))))")
    L.append(f"(assert (= x{k+1}_2 (ite (>= {v} (- (- {c}) {d})) (- (- {c}) {d}) {v})))")
    L.append(f"(assert (= x{k+1}_3 (- (- {v}) {c})))")
L.append("(assert (>= x0_1 x0_2))")
for i in range(4): L.append(f"(assert (and (<= x0_{i} 1) (>= x0_{i} (- 1))))")
L.append("(assert (or "+" ".join(f"(= x0_{i} 1) (= x0_{i} (- 1))" for i in range(4))+"))")
for i in range(4): L.append(f"(assert (and (<= x{n}_{i} {bound}) (>= x{n}_{i} (- {bound}))))")
L.append("(check-sat)")
fn="query_tmp.smt2"; open(fn,"w").write("\n".join(L)+"\n")
print(subprocess.run(["yices-smt2",fn],capture_output=True,text=True).stdout.strip())
