import FiveDigit.Map

namespace FiveDigit

def box (p : Pair) : Prop := -2 ≤ p.1 ∧ p.1 ≤ 9 ∧ -2 ≤ p.2 ∧ p.2 ≤ 9

theorem box_invariance (p : Pair) (hp : box p) :
    box (marginStep (marginStep p)) ∧ height (marginStep p) ≤ 10 := by
  grind [marginStep, box, height]

def cycles : List (List Pair) :=
  [[(0,0)], [(-2,0),(1,3),(3,-2)],
   [(-1,-1),(-1,2),(2,2),(2,-1)], [(-1,0),(0,2),(2,0),(0,-1)]]

/-- Each row follows the map cyclically. Distinct entries across all rows
ensure that each displayed cycle has the displayed least period. -/
theorem listed_cycles :
    cycles.flatten.Nodup ∧ cycles.all (fun row =>
      !row.isEmpty && (row.map marginStep == row.drop 1 ++ row.take 1)) := by decide

/-- First index on a displayed cycle, with fuel+1 as the failure value. -/
def arrival : Nat → Pair → Nat
  | 0, p => if p ∈ cycles.flatten then 0 else 1
  | n+1, p => if p ∈ cycles.flatten then 0 else 1 + arrival n (marginStep p)

theorem arrival_first (n : Nat) (p : Pair) (h : arrival n p ≤ n) :
    iter marginStep (arrival n p) p ∈ cycles.flatten ∧
    ∀ j < arrival n p, iter marginStep j p ∉ cycles.flatten := by
  induction n generalizing p with
  | zero =>
    by_cases hc : p ∈ cycles.flatten <;> simp_all [arrival, iter]
  | succ n ih =>
    by_cases hc : p ∈ cycles.flatten
    · simp [arrival, hc, iter]
    · simp only [arrival, if_neg hc] at h ⊢
      obtain ⟨hi,hj⟩ := ih (marginStep p) (by omega)
      rw [Nat.add_comm 1, iter]
      refine ⟨hi, ?_⟩
      intro j hlt
      cases j with
      | zero => exact hc
      | succ j => exact hj j (by omega)

def scaledError (p : Pair) : Int := 3 * (arrival 10 p : Int) - 2 * height p

private theorem all_starts : ∀ i j : Fin 12,
    let p : Pair := ((i.val : Int)-2, (j.val : Int)-2)
    arrival 10 p ≤ 10 ∧ -9 ≤ scaledError p ∧ scaledError p ≤ 16 := by decide

theorem box_bounds (p : Pair) (hp : box p) :
    arrival 10 p ≤ 10 ∧ -9 ≤ scaledError p ∧ scaledError p ≤ 16 := by
  obtain ⟨ha,ha',hb,hb'⟩ := hp
  let i : Fin 12 := ⟨(p.1+2).toNat, by omega⟩
  let j : Fin 12 := ⟨(p.2+2).toNat, by omega⟩
  have hi : (i.val : Int)-2 = p.1 := by dsimp [i]; omega
  have hj : (j.val : Int)-2 = p.2 := by dsimp [j]; omega
  simpa [hi,hj] using all_starts i j

theorem box_extrema :
    box (5,7) ∧ arrival 10 (5,7) = 10 ∧ scaledError (5,7) = 16 ∧
    box (9,4) ∧ scaledError (9,4) = -9 := by
  unfold box
  decide

def periodic (p : Pair) : Prop := ∃ n, iter marginStep (n+1) p = p

private theorem iter_add (a b : Nat) (p : Pair) :
    iter marginStep (a+b) p = iter marginStep b (iter marginStep a p) := by
  induction a generalizing p with
  | zero => simp [iter]
  | succ a ih => simpa only [Nat.succ_add, iter] using ih (marginStep p)

private theorem return_multiple {p : Pair} {k : Nat}
    (ret : iter marginStep k p = p) (n : Nat) : iter marginStep (k*n) p = p := by
  induction n with
  | zero => rfl
  | succ n ih => rw [Nat.mul_succ, iter_add, ih, ret]

private theorem listed_orbits : ∀ p ∈ cycles.flatten,
    marginStep p ∈ cycles.flatten ∧ iter marginStep 12 p = p := by decide

private theorem listed_forward {p : Pair} (hp : p ∈ cycles.flatten) (n : Nat) :
    iter marginStep n p ∈ cycles.flatten := by
  induction n generalizing p with
  | zero => exact hp
  | succ n ih => exact ih (listed_orbits p hp).1

private theorem periodic_before_hit {p : Pair} (n : Nat)
    (hit : iter marginStep n p ∈ cycles.flatten) (hp : periodic p) : p ∈ cycles.flatten := by
  obtain ⟨k,ret⟩ := hp
  have h := listed_forward hit (k*n)
  rw [← iter_add] at h
  have he : n + k*n = (k+1)*n := by rw [Nat.succ_mul, Nat.add_comm]
  rwa [he, return_multiple ret] at h

/-- The time bounded in (22) is the least time at any periodic state of the
margin map. This excludes additional cycles reachable from any box start. -/
theorem box_transient (p : Pair) (hp : box p) :
    periodic (iter marginStep (arrival 10 p) p) ∧
    ∀ j < arrival 10 p, ¬ periodic (iter marginStep j p) := by
  obtain ⟨hit,early⟩ := arrival_first 10 p (box_bounds p hp).1
  refine ⟨⟨11,(listed_orbits _ hit).2⟩, ?_⟩
  intro j hj periodicAtJ
  apply early j hj
  apply periodic_before_hit (arrival 10 p - j) _ periodicAtJ
  rwa [← iter_add, Nat.add_sub_of_le (Nat.le_of_lt hj)]

#print axioms box_transient
#print axioms box_invariance
#print axioms listed_cycles
#print axioms arrival_first
#print axioms box_bounds
#print axioms box_extrema
end FiveDigit
