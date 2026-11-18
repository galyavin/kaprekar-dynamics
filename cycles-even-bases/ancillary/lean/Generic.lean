import Digits

namespace KaprekarBoundary

def scale (s : Int) (p : Pair) : Pair := (s * p.1, s * p.2)

def away (B : Int) (p : Pair) : Prop :=
  ∀ j : Nat, ¬ boundary B (iter (genericMap B) j p)

/-- The scaled odd-core region of P1. -/
def oddCore (m : Int) (p : Pair) : Prop :=
  validPair m p ∧ 0 < p.2 ∧ p.2 < p.1 ∧ p.1 % 2 = 1 ∧ p.2 % 2 = 1

private theorem sort_scale (s x y : Int) (hs : 0 < s) :
    sortPair (s*x) (s*y) = scale s (sortPair x y) := by
  by_cases h : x ≥ y
  · have hh : s*x ≥ s*y := Int.mul_le_mul_of_nonneg_left h (by omega)
    simp [sortPair, scale, h, hh]
  · have hh : ¬ s*x ≥ s*y := by have := Int.mul_lt_mul_of_pos_left (show x < y by omega) hs; omega
    simp [sortPair, scale, h, hh]

private theorem fold_scale (s B x : Int) (hs : 0 < s) :
    absInt (2*(s*x)-s*B) = s*absInt (2*x-B) := by
  have id : 2*(s*x)-s*B = s*(2*x-B) := by grind
  rw [id]
  by_cases h : 2*x-B < 0
  · have h' := Int.mul_neg_of_pos_of_neg hs h
    simp [absInt, h, h', Int.mul_neg]
  · have h' := Int.mul_nonneg (show 0 ≤ s by omega) (show 0 ≤ 2*x-B by omega)
    simp [absInt, h, show ¬ s*(2*x-B) < 0 by omega]

theorem generic_scale (s B : Int) (hs : 0 < s) (p : Pair) :
    genericMap (s*B) (scale s p) = scale s (genericMap B p) := by
  simp only [genericMap, scale, fold_scale s B _ hs, sort_scale s _ _ hs]

private theorem iter_generic_scale (s B : Int) (hs : 0 < s) (p : Pair) (k : Nat) :
    iter (genericMap (s*B)) k (scale s p) = scale s (iter (genericMap B) k p) := by
  induction k with
  | zero => rfl
  | succ k ih => simp only [iter, ih, generic_scale s B hs]

private theorem scale_inj (s : Int) (hs : 0 < s) (p q : Pair) :
    scale s p = scale s q ↔ p = q := by
  rcases p with ⟨a,b⟩; rcases q with ⟨c,d⟩
  simp [scale, show s ≠ 0 by omega, Int.mul_eq_mul_left_iff]

private theorem boundary_scale (s B : Int) (hs : 0 < s) (p : Pair) :
    boundary (s*B) (scale s p) ↔ boundary B p := by
  simp [boundary, scale, Int.mul_eq_zero, ← Int.mul_add, show s ≠ 0 by omega, Int.mul_eq_mul_left_iff]

theorem generic_actual_iter (B : Int) (p : Pair) (ha : away B p) (k : Nat) :
    iter (pairMap B) k p = iter (genericMap B) k p := by
  induction k with
  | zero => rfl
  | succ k ih => simp only [iter, ih, map_generic _ _ (ha k)]

theorem away_of_actual (B : Int) (p : Pair)
    (ha : ∀ j : Nat, ¬ boundary B (iter (pairMap B) j p)) : away B p := by
  have id (k : Nat) : iter (genericMap B) k p = iter (pairMap B) k p := by
    induction k with
    | zero => rfl
    | succ k ih => simp only [iter, ih, map_generic _ _ (ha k)]
  intro j
  rw [id]
  exact ha j

theorem even_generic (B : Int) (p : Pair) :
    ∃ t, genericMap (2*B) p = scale 2 t := by
  let t := sortPair (absInt (p.1-B)) (absInt (p.2-B))
  refine ⟨t, ?_⟩
  simp only [t, genericMap, sortPair, scale, absInt]
  grind

private theorem odd_generic (B : Int) (hm : B % 2 = 1) (p : Pair) :
    (genericMap B p).1 % 2 = 1 ∧ (genericMap B p).2 % 2 = 1 := by
  grind [genericMap, sortPair, absInt]

private theorem descend (B : Int) (p : Pair) (k : Nat)
    (hp : validPair (2*B) p) (ha : away (2*B) p)
    (hk : iter (genericMap (2*B)) (k+1) p = p) :
    ∃ t, p = scale 2 t ∧ validPair B t ∧ away B t ∧ iter (genericMap B) (k+1) t = t := by
  obtain ⟨t, ht⟩ := even_generic B (iter (genericMap (2*B)) k p)
  have he : p = scale 2 t := hk.symm.trans ht
  subst p
  refine ⟨t, rfl, ?_, ?_, ?_⟩
  · grind [validPair, scale]
  · intro j h
    apply ha j
    rw [iter_generic_scale 2 B (by decide)]
    exact (boundary_scale 2 B (by decide) _).mpr h
  · rw [iter_generic_scale 2 B (by decide)] at hk
    exact (scale_inj 2 (by decide) _ _).mp hk

theorem generic_cycle_descends (m : Int) (hm : m % 2 = 1) (n : Nat)
    (p : Pair) (k : Nat) (hp : validPair (2^n*m) p) (ha : away (2^n*m) p)
    (hk : iter (genericMap (2^n*m)) (k+1) p = p) :
    ∃ t, p = scale (2^n) t ∧ oddCore m t ∧ iter (genericMap m) (k+1) t = t := by
  induction n generalizing p with
  | zero =>
    simp only [Int.pow_zero, Int.one_mul] at hp ha hk ⊢
    have hfree : ¬ boundary m p := ha 0
    have hpar := odd_generic m hm (iter (genericMap m) k p)
    change genericMap m (iter (genericMap m) k p) = p at hk
    rw [hk] at hpar
    refine ⟨p, by simp [scale], ?_, hk⟩
    grind [oddCore, validPair, boundary]
  | succ n ih =>
    have id : (2:Int)^(n+1)*m = 2*(2^n*m) := by rw [Int.pow_succ]; grind
    rw [id] at hp ha hk
    obtain ⟨q, eq, hq, aq, kq⟩ := descend _ p k hp ha hk
    obtain ⟨t, qt, ht, kt⟩ := ih q hq aq kq
    refine ⟨t, ?_, ht, kt⟩
    rw [eq, qt]
    simp [scale, Int.pow_succ]
    grind

theorem oddCore_away (m : Int) (hm : m % 2 = 1) (p : Pair) (hp : oddCore m p) :
    ¬ boundary m p := by
  grind [oddCore, validPair, boundary]

theorem oddCore_preserved (m : Int) (hm : m % 2 = 1) (p : Pair) (hp : oddCore m p) :
    oddCore m (genericMap m p) := by
  have hv := pairMap_valid m p hp.1
  rw [map_generic m p (oddCore_away m hm p hp)] at hv
  have hpar := odd_generic m hm p
  have hone := generic_not_diagonal m p.1 p.2 (by have := hp.2.2.1; omega)
    (by have h := oddCore_away m hm p hp; unfold boundary at h; omega)
  grind [oddCore, validPair]

theorem oddCore_orbit (m : Int) (hm : m % 2 = 1) (p : Pair) (hp : oddCore m p) (j : Nat) :
    oddCore m (iter (genericMap m) j p) := by
  induction j with
  | zero => exact hp
  | succ j ih => exact oddCore_preserved m hm _ ih

private theorem scale_away (s m : Int) (hs : 0 < s) (p : Pair) (ha : away m p) :
    away (s*m) (scale s p) := by
  intro j
  rw [iter_generic_scale s m hs, boundary_scale s m hs]
  exact ha j

/-- The non-boundary part of the classification, including every return time.
On the right, the smaller-base orbit is precisely in the odd core. -/
theorem generic_cycle_iff (m : Int) (hm : m % 2 = 1) (n : Nat) (p : Pair)
    (hp : validPair (2^n*m) p) (k : Nat) :
    ((∀ j, ¬ boundary (2^n*m) (iter (pairMap (2^n*m)) j p)) ∧
      iter (pairMap (2^n*m)) (k+1) p = p) ↔
    ∃ t, p = scale (2^n) t ∧ oddCore m t ∧ iter (pairMap m) (k+1) t = t := by
  have coreAway (t : Pair) (ht : oddCore m t) : away m t :=
    fun j => oddCore_away m hm _ (oddCore_orbit m hm t ht j)
  constructor
  · rintro ⟨ha,hk⟩
    have hfree := away_of_actual _ p ha
    rw [generic_actual_iter _ _ hfree] at hk
    obtain ⟨t, he, ht, kt⟩ := generic_cycle_descends m hm n p k hp hfree hk
    exact ⟨t, he, ht, (generic_actual_iter m t (coreAway t ht) _).trans kt⟩
  · rintro ⟨t, rfl, ht, kt⟩
    have hs : (0:Int) < 2^n := Int.pow_pos (by decide)
    have hfree := scale_away (2^n) m hs t (coreAway t ht)
    refine ⟨fun j => ?_, ?_⟩
    · rw [generic_actual_iter _ _ hfree]
      exact hfree j
    · rw [generic_actual_iter _ _ hfree, iter_generic_scale _ _ hs,
          ← generic_actual_iter _ _ (coreAway t ht), kt]

#print axioms generic_cycle_iff
end KaprekarBoundary
