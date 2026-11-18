import Critical

namespace KaprekarBoundary

private theorem five_power (k : Nat) : (2:Int)^(4*k+2)%10=4 := by
  induction k with
  | zero => decide
  | succ k ih =>
    rw [show 4*(k+1)+2 = (4*k+2)+4 by omega, Int.pow_add]
    change 2^(4*k+2)*16%10=4
    omega

/-- The second exceptional family in P1, with its exact period on digits. -/
theorem five_exception_digit_period (n : Nat) (hn : n%4=2) :
    let s : Int := 2^n
    let B := s*5
    let q := fromPair B (2*s,0)
    iter (step B) (2*(n+1)) q = q ∧ ∀ i, 0<i → i<2*(n+1) → iter (step B) i q ≠ q := by
  let s : Int := 2^n
  have hs : 0<s := Int.pow_pos (by decide)
  have he : s%10=4 := by
    change 2^n%10=4
    have eq : n=4*(n/4)+2 := by omega
    conv => lhs; rw [eq]
    exact five_power _
  have hp : validPair (s*5) (2*s,0) := by simp [validPair]; omega
  have hu : validPair (s*5) (3*s,2*s) := by simp [validPair]; omega
  have first := zero_critical_return 5 (2*s) (by decide) (by decide) n (by omega)
    (by omega) (by omega) hp (by omega)
  have second := diag_critical_return 5 (3*s) (2*s) (by decide) (by decide) n (by omega)
    (by omega) (Or.inr (show 3*s+2*s=s*5 by omega)) (by omega) hu (by omega)
  have a : triangle 5 (2*s-1)=2 := by unfold triangle absInt; split <;> omega
  have b : triangle 5 (2*s)=3 := by unfold triangle absInt; split <;> omega
  have c : triangle 5 (2*(3*s)+1)=0 := by unfold triangle absInt; split <;> omega
  have d : triangle 5 (2*(3*s)-1)=2 := by unfold triangle absInt; split <;> omega
  change FirstReturn _ _ _ (sortPair (s*_) (s*_)) _ at first second
  simp only [a,b,c,d,Int.mul_zero] at first second
  have eq1 : sortPair (s*2) (s*3) = (3*s,2*s) := by grind [sortPair]
  have eq2 : sortPair 0 (s*2) = (2*s,0) := by grind [sortPair]
  rw [eq1] at first
  rw [eq2] at second
  have per := two_returns_exact_period _ _ _ _ _ _ (by grind : (2*s,(0:Int)) ≠ (3*s,2*s)) first second
  rw [show n+1+(n+1)=2*(n+1) by omega] at per
  exact ⟨(digit_return_iff _ _ hp _).mpr per.1,
    fun i hi hj h => per.2 i hi hj ((digit_return_iff _ _ hp i).mp h)⟩

private theorem three_power (n : Nat) (hn : 0<n) :
    (2:Int)^n%6 = if n%2=0 then 4 else 2 := by
  cases n with
  | zero => omega
  | succ k =>
    induction k with
    | zero => decide
    | succ k ih =>
      rw [Int.pow_succ]
      split at ih <;> split <;> omega

private theorem three_returns {α : Type} (f : α → α) (b : α → Prop) (x y z : α) (k : Nat)
    (hxy : x≠y) (hxz : x≠z) (a : FirstReturn f b x y k)
    (c : FirstReturn f b y z k) (d : FirstReturn f b z x k) :
    iter f (3*k) x=x ∧ ∀ i, 0<i → i<3*k → iter f i x≠x := by
  have one := a.2.2.2.1
  have two : iter f (k+k) x=z := by rw [iter_add, one, c.2.2.2.1]
  constructor
  · rw [show 3*k=(k+k)+k by omega, iter_add, two, d.2.2.2.1]
  · intro i hi hsize hret
    by_cases hfirst : i<k
    · exact a.2.2.2.2 i hi hfirst (by rw [hret]; exact a.1)
    by_cases heq : i=k
    · exact hxy (hret.symm.trans (heq ▸ one))
    by_cases hsecond : i<k+k
    · have : iter f (i-k) y=x := by
        rw [← one, ← iter_add, show k+(i-k)=i by omega, hret]
      exact c.2.2.2.2 (i-k) (by omega) (by omega) (by rw [this]; exact a.1)
    by_cases heq2 : i=k+k
    · exact hxz (hret.symm.trans (heq2 ▸ two))
    have : iter f (i-(k+k)) z=x := by
      rw [← two, ← iter_add, show k+k+(i-(k+k))=i by omega, hret]
    exact d.2.2.2.2 (i-(k+k)) (by omega) (by omega) (by rw [this]; exact a.1)

/-- The critical graph of the odd core three, including the clock of each edge. -/
theorem three_critical_graph (n : Nat) (hn : 0<n) :
    let s : Int := 2^n
    let S : Pair := (s,0)
    let E : Pair := (2*s,0)
    let U : Pair := (2*s,s)
    let R := fun p q => FirstReturn (pairMap (s*3)) (boundary (s*3)) p q (n+1)
    R U E ∧ (if n%2=0 then R S S ∧ R E U else R S U ∧ R E S) := by
  let s : Int := 2^n
  have hs : 0<s := Int.pow_pos (by decide)
  have hres : s%6=if n%2=0 then 4 else 2 := three_power n hn
  have hsV : validPair (s*3) (s,0) := by simp [validPair]; omega
  have heV : validPair (s*3) (2*s,0) := by simp [validPair]; omega
  have huV : validPair (s*3) (2*s,s) := by simp [validPair]; omega
  have res : s%6=2 ∨ s%6=4 := by split at hres <;> omega
  have rs := zero_critical_return 3 s (by decide) (by decide) n hn (by omega) (by omega) hsV (by omega)
  have re := zero_critical_return 3 (2*s) (by decide) (by decide) n hn (by omega) (by omega) heV (by omega)
  have ru := diag_critical_return 3 (2*s) s (by decide) (by decide) n hn
    (by omega) (Or.inr (show 2*s+s=s*3 by omega)) (by omega) huV (by omega)
  change FirstReturn _ _ _ (sortPair (s*_) (s*_)) _ at rs re ru
  have tri0 (x : Int) (h : x%6=3) : triangle 3 x=0 := by unfold triangle absInt; omega
  have tri1 (x : Int) (h : x%6=2 ∨ x%6=4) : triangle 3 x=1 := by unfold triangle absInt; split <;> omega
  have tri2 (x : Int) (h : x%6=1 ∨ x%6=5) : triangle 3 x=2 := by unfold triangle absInt; split <;> omega
  dsimp only
  split at hres <;> rename_i hpar
  · simp only [if_pos hpar]
    rw [tri0 _ (by omega), tri1 _ (by omega)] at rs
    rw [tri2 _ (by omega), tri1 _ (by omega)] at re
    rw [tri2 _ (by omega), tri0 _ (by omega)] at ru
    simp only [Int.mul_zero, Int.mul_one, Int.mul_comm s 2] at rs re ru
    have e1 : sortPair 0 s = (s,0) := by grind [sortPair]
    have e2 : sortPair (2*s) s = (2*s,s) := by grind [sortPair]
    have e3 : sortPair (2*s) 0 = (2*s,0) := by grind [sortPair]
    have e4 : sortPair 0 (2*s) = (2*s,0) := by grind [sortPair]
    grind only
  · simp only [if_neg hpar]
    rw [tri2 _ (by omega), tri1 _ (by omega)] at rs
    rw [tri0 _ (by omega), tri1 _ (by omega)] at re
    rw [tri0 _ (by omega), tri2 _ (by omega)] at ru
    simp only [Int.mul_zero, Int.mul_one, Int.mul_comm s 2] at rs re ru
    have e1 : sortPair 0 s = (s,0) := by grind [sortPair]
    have e2 : sortPair (2*s) s = (2*s,s) := by grind [sortPair]
    have e3 : sortPair (2*s) 0 = (2*s,0) := by grind [sortPair]
    have e4 : sortPair 0 (2*s) = (2*s,0) := by grind [sortPair]
    grind only

/-- Exact positive period of an actual digit word. -/
def digitPeriod (B : Int) (p : Pair) (k : Nat) : Prop :=
  0<k ∧ iter (step B) k (fromPair B p)=fromPair B p ∧
    ∀ i, 0<i → i<k → iter (step B) i (fromPair B p)≠fromPair B p

theorem digitPeriod_of_pair (B : Int) (p : Pair) (hp : validPair B p) (k : Nat) (hk : 0<k)
    (h : iter (pairMap B) k p=p ∧ ∀ i, 0<i → i<k → iter (pairMap B) i p≠p) :
    digitPeriod B p k := ⟨hk, (digit_return_iff _ _ hp _).mpr h.1,
      fun i hi hj he => h.2 i hi hj ((digit_return_iff _ _ hp i).mp he)⟩

theorem three_digit_periods (n : Nat) (hn : 0<n) :
    let s : Int := 2^n
    let B := s*3
    if n%2=0 then digitPeriod B (s,0) (n+1) ∧ digitPeriod B (2*s,0) (2*(n+1))
    else digitPeriod B (s,0) (3*(n+1)) := by
  let s : Int := 2^n
  have hs : 0<s := Int.pow_pos (by decide)
  have hS : validPair (s*3) (s,0) := by simp [validPair]; omega
  have hE : validPair (s*3) (2*s,0) := by simp [validPair]; omega
  obtain ⟨ue, h⟩ := three_critical_graph n hn
  dsimp only
  split <;> rename_i hp
  · simp only [if_pos hp] at h
    obtain ⟨ss, eu⟩ := h
    have one := first_return_self_exact_period _ _ _ _ ss
    have two := two_returns_exact_period _ _ _ _ _ _
      (by grind : (2*s,(0:Int)) ≠ (2*s,s)) eu ue
    rw [show n+1+(n+1)=2*(n+1) by omega] at two
    exact ⟨digitPeriod_of_pair _ _ hS _ one.1 one.2, digitPeriod_of_pair _ _ hE _ (by omega) two⟩
  · simp only [if_neg hp] at h
    obtain ⟨su, es⟩ := h
    have three := three_returns _ _ _ _ _ _
      (by grind : (s,(0:Int)) ≠ (2*s,s)) (by grind : (s,(0:Int)) ≠ (2*s,0)) su ue es
    exact digitPeriod_of_pair _ _ hS _ (by omega) three

#print axioms three_critical_graph
#print axioms three_digit_periods
#print axioms five_exception_digit_period
end KaprekarBoundary
