import Generic

namespace KaprekarBoundary

/-- Boundary walls in the two coordinates modulo the base.  Evaluating this
test at doubled coordinates gives the time to the first boundary visit. -/
def wall (B x y : Int) : Prop := B ∣ x ∨ B ∣ y ∨ B ∣ x+y ∨ B ∣ x-y

private theorem negdiff (B x y : Int) : B ∣ y-x ↔ B ∣ x-y := by
  rw [show y-x = -(x-y) by omega, Int.dvd_neg]

private theorem wall_swap (B x y : Int) : wall B y x ↔ wall B x y := by
  simp only [wall, Int.add_comm y x, negdiff B x y]
  grind

private theorem wall_neg (B x y : Int) : wall B (-x) y ↔ wall B x y := by
  simp only [wall, show -x+y = -(x-y) by omega, show -x-y = -(x+y) by omega, Int.dvd_neg]
  grind

private theorem wall_shift (B x y a : Int) : wall B (x+a*B) y ↔ wall B x y := by
  simp only [wall, show x+a*B+y = (x+y)+a*B by omega,
    show x+a*B-y = (x-y)+a*B by omega, Int.dvd_add_mul_self]

private theorem wall_fold_left (B d x y : Int) :
    wall B (d*absInt (2*x-B)) y ↔ wall B ((d*2)*x) y := by
  have eq : d*(2*x-B) = (d*2)*x + (-d)*B := by grind
  unfold absInt
  split
  · rw [Int.mul_neg, wall_neg, eq, wall_shift]
  · rw [eq, wall_shift]

private theorem wall_sort (B d x y : Int) :
    wall B (d*(sortPair x y).1) (d*(sortPair x y).2) ↔ wall B (d*x) (d*y) := by
  unfold sortPair
  split
  · rfl
  · exact wall_swap B _ _

theorem wall_generic (B d : Int) (p : Pair) :
    wall B (d*(genericMap B p).1) (d*(genericMap B p).2) ↔
      wall B ((d*2)*p.1) ((d*2)*p.2) := by
  unfold genericMap
  rw [wall_sort, wall_fold_left, wall_swap, wall_fold_left, wall_swap]

theorem wall_iter_generic (B d : Int) (p : Pair) (j : Nat) :
    wall B (d*(iter (genericMap B) j p).1) (d*(iter (genericMap B) j p).2) ↔
      wall B ((d*2^j)*p.1) ((d*2^j)*p.2) := by
  induction j generalizing d with
  | zero => simp [iter]
  | succ j ih =>
    simp only [iter, wall_generic, ih, Int.pow_succ]
    rw [show d*2*2^j = d*(2^j*2) by grind]

private theorem dvd_small (B x : Int) (hn : 0 ≤ x) (hu : x < B) : B ∣ x ↔ x = 0 := by
  constructor
  · intro hd
    by_cases hx : x = 0
    · exact hx
    · have := Int.le_of_dvd (show 0 < x by omega) hd
      omega
  · rintro rfl; simp

theorem wall_iff_boundary (B : Int) (p : Pair) (hp : validPair B p) :
    wall B p.1 p.2 ↔ boundary B p := by
  rcases p with ⟨u,v⟩
  obtain ⟨hu,hB,hv,hvu⟩ := hp
  simp only [wall, boundary, dvd_small B u (by omega) hB, dvd_small B v hv (by omega),
    dvd_small B (u-v) (by omega) (by omega)]
  have sum : B ∣ u+v ↔ u+v=B := by
    by_cases h : u+v < B
    · rw [dvd_small B _ (by omega) h]; omega
    · rw [show u+v = (u+v-B)+1*B by omega, Int.dvd_add_mul_self,
           dvd_small B _ (by omega) (by omega)]
      omega
  grind

theorem actual_generic_prefix (B : Int) (p : Pair) (k : Nat)
    (h : ∀ i, i < k → ¬ boundary B (iter (pairMap B) i p)) :
    iter (pairMap B) k p = iter (genericMap B) k p := by
  induction k with
  | zero => rfl
  | succ k ih =>
    simp only [iter]
    rw [map_generic _ _ (h k (by omega)), ih (fun i hi => h i (by omega))]

/-- The first boundary is reached at exactly the first nonnegative generic
time at which a divisibility test fires.  The start here is the immediate
output of a boundary or any other valid pair. -/
theorem boundary_clock (B : Int) (p : Pair) (hp : validPair B p) (j : Nat) :
    (boundary B (iter (pairMap B) j p) ∧
      ∀ i, i < j → ¬ boundary B (iter (pairMap B) i p)) ↔
    (wall B (2^j*p.1) (2^j*p.2) ∧ ∀ i : Nat, i < j → ¬ wall B (2^i*p.1) (2^i*p.2)) := by
  have test (k : Nat) (h : ∀ i, i < k → ¬ boundary B (iter (pairMap B) i p)) :
      boundary B (iter (pairMap B) k p) ↔ wall B (2^k*p.1) (2^k*p.2) := by
    rw [← wall_iff_boundary _ _ (iter_pair_valid B p hp k), actual_generic_prefix _ _ k h]
    simpa only [Int.one_mul] using wall_iter_generic B 1 p k
  constructor
  · rintro ⟨hj,h⟩
    exact ⟨(test j h).mp hj, fun i hi hwall =>
      h i hi ((test i (fun l hl => h l (by omega))).mpr hwall)⟩
  · rintro ⟨hj,h⟩
    have hfree (k : Nat) (hk : k < j) : ¬ boundary B (iter (pairMap B) k p) := by
      induction k using Nat.strongRecOn with
      | ind k ih =>
        exact fun hb => h k hk ((test k (fun l hl => ih l hl (by omega))).mp hb)
    exact ⟨(test j hfree).mpr hj, hfree⟩

/-- Positive-time first return from any boundary state. -/
theorem first_return_clock (B : Int) (p : Pair) (hp : validPair B p) (hb : boundary B p) (j : Nat) :
    (∃ q, FirstReturn (pairMap B) (boundary B) p q (j+1)) ↔
    (wall B (2^j*(pairMap B p).1) (2^j*(pairMap B p).2) ∧
      ∀ i : Nat, i < j → ¬ wall B (2^i*(pairMap B p).1) (2^i*(pairMap B p).2)) := by
  rw [← boundary_clock _ _ (pairMap_valid _ _ hp)]
  have advance (i : Nat) : iter (pairMap B) i (pairMap B p) = iter (pairMap B) (i+1) p := by
    change iter (pairMap B) i (iter (pairMap B) 1 p) = _
    rw [← iter_add, Nat.add_comm]
  constructor
  · rintro ⟨q, _, hq, _, he, hmin⟩
    rw [advance, he]
    exact ⟨hq, fun i hi => by rw [advance]; exact hmin _ (by omega) (by omega)⟩
  · rintro ⟨hj,hmin⟩
    refine ⟨_, hb, hj, by omega, (advance j).symm, ?_⟩
    rintro (_ | i) hpos hi
    · omega
    · rw [← advance]; exact hmin i (by omega)

private theorem wall_abs_left (B d x y : Int) :
    wall B (d*absInt x) y ↔ wall B (d*x) y := by
  unfold absInt
  split
  · rw [Int.mul_neg, wall_neg]
  · rfl

private theorem wall_shift_right (B x y a : Int) : wall B x (y+a*B) ↔ wall B x y := by
  rw [wall_swap, wall_shift, wall_swap]

theorem zero_wall (B d u : Int) :
    wall B (d*(pairMap B (u,0)).1) (d*(pairMap B (u,0)).2) ↔
    B ∣ d*(u-1) ∨ B ∣ d*u ∨ B ∣ d*(2*u-1) ∨ B ∣ d := by
  simp only [pairMap, if_true, wall_sort]
  rw [show d*(B-u) = -(d*u)+d*B by grind, wall_shift_right,
      wall_swap, wall_neg, wall_swap]
  simp only [wall, show d*(u-1)+d*u = d*(2*u-1) by grind,
    show d*(u-1)-d*u = -d by grind, Int.dvd_neg]

theorem diagonal_wall (B d u v : Int) (hv : v ≠ 0) (hb : u=v ∨ u+v=B) :
    wall B (d*(pairMap B (u,v)).1) (d*(pairMap B (u,v)).2) ↔
    B ∣ d*(2*u+1) ∨ B ∣ d*(2*u-1) ∨ B ∣ d*(4*u) ∨ B ∣ 2*d := by
  simp only [pairMap, if_neg hv, if_pos hb, wall_sort]
  rw [wall_abs_left, wall_swap, wall_abs_left, wall_swap]
  rw [show d*(2*u-B+1) = d*(2*u+1)+(-d)*B by grind, wall_shift,
      show d*(2*u-B-1) = d*(2*u-1)+(-d)*B by grind, wall_shift_right]
  simp only [wall, show d*(2*u+1)+d*(2*u-1) = d*(4*u) by grind,
    show d*(2*u+1)-d*(2*u-1) = 2*d by grind]

/-- Fold a residue modulo twice the base around its midpoint. -/
def triangle (B x : Int) : Int := absInt (x % (2*B) - B)

private theorem fold_twice (B : Int) (hB : 0 < B) (x : Int) :
    absInt (2*triangle B x-B) = triangle B (2*x) := by
  let r := x % (2*B)
  have lo : 0 ≤ r := Int.emod_nonneg _ (by omega)
  have hi : r < 2*B := Int.emod_lt_of_pos _ (by omega)
  have red : (2*x) % (2*B) = (2*r) % (2*B) := by
    simp [r, Int.mul_emod]
  unfold triangle
  rw [red]
  change absInt (2*absInt (r-B)-B) = absInt ((2*r) % (2*B)-B)
  by_cases hr : r < B
  · rw [Int.emod_eq_of_lt (by omega) (by omega)]
    grind [absInt]
  · have sub : (2*r)%(2*B) = 2*r-2*B := by
      rw [← Int.emod_eq_of_lt (a:=2*r-2*B) (b:=2*B) (by omega) (by omega)]
      simp
    rw [sub]
    grind [absInt]

private theorem generic_sort (B x y : Int) :
    genericMap B (sortPair x y) = sortPair (absInt (2*x-B)) (absInt (2*y-B)) := by
  grind [genericMap, sortPair]

theorem generic_closed_form (B : Int) (p : Pair) (hp : validPair B p) (j : Nat) :
    iter (genericMap B) (j+1) p =
      sortPair (triangle B (2^(j+1)*p.1)) (triangle B (2^(j+1)*p.2)) := by
  obtain ⟨hu,hub,hv,hvu⟩ := hp
  induction j with
  | zero =>
    simp only [iter, Nat.zero_add, Int.pow_one, genericMap, triangle]
    rw [Int.emod_eq_of_lt (a:=2*p.1) (b:=2*B) (by omega) (by omega),
        Int.emod_eq_of_lt (a:=2*p.2) (b:=2*B) (by omega) (by omega)]
  | succ j ih =>
    change genericMap B (iter (genericMap B) (j+1) p) = _
    rw [ih, generic_sort, fold_twice _ (by omega), fold_twice _ (by omega)]
    rw [Int.pow_succ (n:=j+1)]
    congr 2 <;> grind

/-- After the boundary branch, zero generic steps leave its output unchanged.
Positive generic times have a closed expression in the two output digits. -/
def returnTarget (B : Int) (p : Pair) (j : Nat) : Pair :=
  let q := pairMap B p
  match j with
  | 0 => q
  | k+1 => sortPair (triangle B (2^(k+1)*q.1)) (triangle B (2^(k+1)*q.2))

theorem first_return_iff (B : Int) (p q : Pair) (hp : validPair B p) (hb : boundary B p) (j : Nat) :
    FirstReturn (pairMap B) (boundary B) p q (j+1) ↔
    q = returnTarget B p j ∧ wall B (2^j*(pairMap B p).1) (2^j*(pairMap B p).2) ∧
      ∀ i : Nat, i < j → ¬ wall B (2^i*(pairMap B p).1) (2^i*(pairMap B p).2) := by
  have target (z : Pair) (hr : FirstReturn (pairMap B) (boundary B) p z (j+1)) :
      z = returnTarget B p j := by
    have advance (i : Nat) : iter (pairMap B) i (pairMap B p) = iter (pairMap B) (i+1) p := by
      change iter (pairMap B) i (iter (pairMap B) 1 p) = _
      rw [← iter_add, Nat.add_comm]
    have hfree : ∀ i, i < j → ¬ boundary B (iter (pairMap B) i (pairMap B p)) :=
      fun i hi => by rw [advance]; exact hr.2.2.2.2 _ (by omega) (by omega)
    rw [← hr.2.2.2.1, ← advance, actual_generic_prefix _ _ j hfree]
    cases j with
    | zero => rfl
    | succ k => exact generic_closed_form B _ (pairMap_valid B p hp) k
  constructor
  · intro hr
    exact ⟨target _ hr, (first_return_clock B p hp hb j).mp ⟨q, hr⟩⟩
  · rintro ⟨hq,hw⟩
    obtain ⟨z,hz⟩ := (first_return_clock B p hp hb j).mpr hw
    rw [hq, ← target z hz]
    exact hz

#print axioms first_return_iff
#print axioms generic_closed_form
#print axioms first_return_clock
#print axioms zero_wall
#print axioms diagonal_wall
#print axioms boundary_clock
end KaprekarBoundary
