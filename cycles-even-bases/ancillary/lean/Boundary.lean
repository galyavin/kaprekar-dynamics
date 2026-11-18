import Std

/-!
Self-contained, limited formalization of the four-digit boundary dependency.
The map below is the credited piecewise pair map. Its bridge to sorted digit
subtraction, the full six-row return table, and Theorem C are NOT proved here.
No result in this file assumes any of those three missing results.
-/
namespace KaprekarBoundary

abbrev Pair := Int × Int
def absInt (x : Int) : Int := if x < 0 then -x else x
def sortPair (a b : Int) : Pair := if a ≥ b then (a,b) else (b,a)
def boundary (B : Int) (p : Pair) : Prop :=
  p.2 = 0 ∨ p.1 = p.2 ∨ p.1 + p.2 = B
def genericMap (B : Int) (p : Pair) : Pair :=
  sortPair (absInt (2*p.1-B)) (absInt (2*p.2-B))
def pairMap (B : Int) (p : Pair) : Pair :=
  if p.2 = 0 then sortPair (p.1 - 1) (B - p.1)
  else if p.1 = p.2 ∨ p.1 + p.2 = B then
    sortPair (absInt (2*p.1-B+1)) (absInt (2*p.1-B-1))
  else genericMap B p

theorem zero_one (B : Int) (hB : B ≥ 4) :
    pairMap B (1,0) = (B-1,0) := by
  simp only [pairMap, sortPair]
  simp
  omega

theorem diagonal_midpoint (h : Int) (hh : h ≥ 2) :
    pairMap (2*h) (h,h) = (1,1) := by
  simp [pairMap, sortPair, absInt, show h ≠ 0 by omega]

theorem generic_not_diagonal (B u v : Int)
    (hne : u ≠ v) (hsum : u + v ≠ B) :
    (genericMap B (u,v)).1 ≠ (genericMap B (u,v)).2 := by
  simp only [genericMap, sortPair, absInt]
  split <;> split <;> split <;> simp_all <;> omega

theorem map_generic (B : Int) (p : Pair) (hp : ¬ boundary B p) :
    pairMap B p = genericMap B p := by
  simp only [boundary, not_or] at hp
  simp [pairMap, hp.1, hp.2.1, hp.2.2]

theorem diagonal_has_no_generic_predecessor (B : Int) (p : Pair)
    (hp : ¬ boundary B p) : (pairMap B p).1 ≠ (pairMap B p).2 := by
  rw [map_generic B p hp]
  simp only [boundary, not_or] at hp
  exact generic_not_diagonal B p.1 p.2 hp.2.1 hp.2.2

def iter {α : Type} (f : α → α) : Nat → α → α
  | 0, x => x
  | n+1, x => f (iter f n x)

theorem iter_add {α : Type} (f : α → α) (a b : Nat) (x : α) :
    iter f (a+b) x = iter f b (iter f a x) := by
  induction b with
  | zero => simp [iter]
  | succ b ih => simp [iter, ih]

def FirstReturn {α : Type} (f : α → α) (S : α → Prop)
    (x y : α) (n : Nat) : Prop :=
  S x ∧ S y ∧ 0 < n ∧ iter f n x = y ∧
    ∀ k, 0 < k → k < n → ¬ S (iter f k x)

theorem direct_first_return {α : Type} (f : α → α) (S : α → Prop)
    (x y : α) (hx : S x) (hy : S y) (hf : f x = y) :
    FirstReturn f S x y 1 := by
  refine ⟨hx, hy, by omega, ?_, ?_⟩
  · simpa [iter] using hf
  · intro k hk hk1
    omega

theorem zero_one_first_return (B : Int) (hB : B ≥ 4) :
    FirstReturn (pairMap B) (boundary B) (1,0) (B-1,0) 1 := by
  apply direct_first_return
  · exact Or.inl rfl
  · exact Or.inl rfl
  · exact zero_one B hB

theorem diagonal_midpoint_first_return (h : Int) (hh : h ≥ 2) :
    FirstReturn (pairMap (2*h)) (boundary (2*h)) (h,h) (1,1) 1 := by
  apply direct_first_return
  · exact Or.inr (Or.inl rfl)
  · exact Or.inr (Or.inl rfl)
  · exact diagonal_midpoint h hh

theorem diagonal_first_return_one (B : Int) (x y : Pair) (n : Nat)
    (h : FirstReturn (pairMap B) (boundary B) x y n) (hy : y.1 = y.2) :
    n = 1 := by
  by_cases hn : n = 1
  · exact hn
  have hk : 0 < n-1 := by have := h.2.2.1; omega
  have hkn : n-1 < n := by omega
  have hprev := h.2.2.2.2 (n-1) hk hkn
  have hne := diagonal_has_no_generic_predecessor B (iter (pairMap B) (n-1) x) hprev
  have he : pairMap B (iter (pairMap B) (n-1) x) = y := by
    have hsucc : n-1+1 = n := by omega
    change iter (pairMap B) (n-1+1) x = y
    rw [hsucc]
    exact h.2.2.2.1
  rw [he] at hne
  exact False.elim (hne hy)

theorem first_return_unique {α : Type} (f : α → α) (S : α → Prop)
    (x y z : α) (n m : Nat)
    (h : FirstReturn f S x y n) (g : FirstReturn f S x z m) :
    n = m ∧ y = z := by
  have he : n = m := by
    by_cases heq : n = m
    · exact heq
    have hlt : n < m ∨ m < n := by omega
    rcases hlt with hlt | hlt
    · exact False.elim (g.2.2.2.2 n h.2.2.1 hlt (by rw [h.2.2.2.1]; exact h.2.1))
    · exact False.elim (h.2.2.2.2 m g.2.2.1 hlt (by rw [g.2.2.2.1]; exact g.2.1))
  constructor
  · exact he
  · rw [← h.2.2.2.1, ← g.2.2.2.1, he]

theorem first_return_self_exact_period {α : Type} (f : α → α)
    (S : α → Prop) (x : α) (n : Nat) (h : FirstReturn f S x x n) :
    0 < n ∧ iter f n x = x ∧
      ∀ k, 0 < k → k < n → iter f k x ≠ x := by
  refine ⟨h.2.2.1, h.2.2.2.1, ?_⟩
  intro k hk hkn he
  exact h.2.2.2.2 k hk hkn (by rw [he]; exact h.1)

theorem two_returns_exact_period {α : Type} (f : α → α)
    (S : α → Prop) (x y : α) (n m : Nat)
    (hxy : x ≠ y) (h : FirstReturn f S x y n)
    (g : FirstReturn f S y x m) :
    iter f (n+m) x = x ∧
      ∀ k, 0 < k → k < n+m → iter f k x ≠ x := by
  constructor
  · rw [iter_add, h.2.2.2.1, g.2.2.2.1]
  · intro k hk hkn he
    by_cases hlt : k < n
    · exact h.2.2.2.2 k hk hlt (by rw [he]; exact h.1)
    by_cases heq : k = n
    · subst k
      exact hxy (he.symm.trans h.2.2.2.1)
    have hpos : 0 < k-n := by omega
    have hsmall : k-n < m := by omega
    have hadd : n+(k-n) = k := by omega
    have hiter : iter f (k-n) y = x := by
      rw [← h.2.2.2.1, ← iter_add, hadd, he]
    exact g.2.2.2.2 (k-n) hpos hsmall (by rw [hiter]; exact h.1)

-- The arithmetic exception behind the E <-> U boundary cycle.
theorem exceptional_core_five (m h : Int)
    (hh : 2*h = m+1) (hc : m-2 = h) : m = 5 := by
  omega

-- Pure powers of two use these involutions on their critical indices.
theorem zero_index_involution (n t : Int) : n-(n-t) = t := by omega
theorem anti_index_involution (n t : Int) : n-(n-t-2)-2 = t := by omega
theorem zero_two_return_clock (n t : Int) :
    (n-t+1) + (n-(n-t)+1) = n+2 := by omega
theorem anti_two_return_clock (n t : Int) :
    (n-t-1) + (n-(n-t-2)-1) = n := by omega

#print axioms zero_one
#print axioms diagonal_midpoint
#print axioms diagonal_has_no_generic_predecessor
#print axioms zero_one_first_return
#print axioms diagonal_midpoint_first_return
#print axioms diagonal_first_return_one
#print axioms first_return_unique
#print axioms first_return_self_exact_period
#print axioms two_returns_exact_period
#print axioms exceptional_core_five
#print axioms zero_index_involution
#print axioms anti_index_involution
#print axioms zero_two_return_clock
#print axioms anti_two_return_clock
end KaprekarBoundary
