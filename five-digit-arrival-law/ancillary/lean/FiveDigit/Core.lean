import FiveDigit.Map

namespace FiveDigit

theorem pair_valid (B : Int) (p : Pair) (h : valid B p) : valid B (pairStep B p) := by
  grind [pairStep, positiveStep, valid]

def corePairStep (B : Int) (p : Pair) : Pair :=
  if 1 ≤ (margins B p).1 then (p.1-1,p.1-p.2+1) else
  if -1 ≤ (margins B p).1 then (B-p.2,2*p.1-B) else (B-p.2,p.1-p.2-1)

theorem core_pair_rule (B : Int) (hB : 30 ≤ B) (p : Pair) (hp : core B p) :
    pairStep B p = corePairStep B p := by
  grind [pairStep, positiveStep, corePairStep, core, valid, margins]

theorem core_rule (B : Int) (hB : 30 ≤ B) (p : Pair) (hp : core B p) :
    core B (pairStep B p) ∧ margins B (pairStep B p) = marginStep (margins B p) := by
  rw [core_pair_rule B hB p hp]
  grind [corePairStep, marginStep, core, valid, margins]

theorem strip_bounds (B : Int) (p : Pair) (hp : valid B p) (hv : 0 < p.2) :
    valid B (pairStep B p) ∧ B-1 ≤ 2*(pairStep B p).1 ∧
      2*((pairStep B p).1-(pairStep B p).2) ≤ B+2 := by
  refine ⟨pair_valid B p hp, ?_⟩
  grind [pairStep, positiveStep, valid]

theorem two_step_descent (p : Pair) (ha : -2 ≤ p.1) (hb : -2 ≤ p.2) (h : 10 ≤ height p) :
    height (marginStep (marginStep p)) = height p - 3 + (if min p.1 p.2 = -2 then 1 else 0) := by
  grind [marginStep, height]

theorem rare_predecessors (m : Int) (hB : 30 ≤ 2*m) (p : Pair) (hp : valid (2*m) p)
    (hi : rare m (pairStep (2*m) p)) : rare m p := by
  rcases p with ⟨u,v⟩
  grind [rare, pairStep, positiveStep, valid]

theorem rare_images (m : Int) (hB : 30 ≤ 2*m) :
    pairStep (2*m) (m,m-1) = (m+1,0) ∧ pairStep (2*m) (m+1,0) = (m,m-1) ∧
    pairStep (2*m) (m,0) = (m,m-1) := by
  grind [pairStep, positiveStep]

theorem core_residue (B : Int) (hB : 30 ≤ B) (hr : B % 3 = 1) (p : Pair) (hp : core B p) :
    (p.1+p.2) % 3 = 0 ↔ ((pairStep B p).1+(pairStep B p).2) % 3 = 0 := by
  rw [core_pair_rule B hB p hp]
  unfold corePairStep
  unfold core margins at hp
  simp only [margins]
  split
  · dsimp; omega
  · split <;> dsimp <;> omega

#print axioms pair_valid
#print axioms core_pair_rule
#print axioms core_rule
#print axioms strip_bounds
#print axioms two_step_descent
#print axioms rare_predecessors
#print axioms rare_images
#print axioms core_residue
end FiveDigit
