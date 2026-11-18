import Clock

namespace KaprekarBoundary

theorem odd_cancel_two (m x : Int) (hm : m%2=1) : m ∣ 2*x ↔ m ∣ x := by
  constructor
  · rintro ⟨q,hq⟩
    have he : 2*((m+1)/2) = m+1 := by omega
    refine ⟨q*((m+1)/2)-x, ?_⟩
    calc
      x = (2*x)*((m+1)/2)-m*x := by grind
      _ = m*(q*((m+1)/2)-x) := by rw [hq]; grind
  · exact Int.dvd_mul_of_dvd_right

theorem odd_cancel_pow (m x : Int) (j : Nat) (hm : m%2=1) : m ∣ 2^j*x ↔ m ∣ x := by
  induction j with
  | zero => simp
  | succ j ih =>
    rw [Int.pow_succ, show 2^j*2*x=2*(2^j*x) by grind, odd_cancel_two m _ hm, ih]

theorem odd_no_early (m x : Int) (n i : Nat) (hi : i < n) (hx : x%2=1) :
    ¬ 2^n*m ∣ 2^i*x := by
  have he : n = i + (n-i-1) + 1 := by omega
  have factor : 2^n*m = 2^i*(2*(2^(n-i-1)*m)) := by
    conv => lhs; rw [he, Int.pow_succ, Int.pow_add]
    grind
  rw [factor, Int.mul_dvd_mul_iff_left (Int.ne_of_gt (Int.pow_pos (by decide : (0:Int)<2)))]
  rintro ⟨k,hk⟩
  have heven : x = 2*((2^(n-i-1)*m)*k) := by rw [hk]; grind
  omega

theorem odd_core_divides (m x : Int) (n i : Nat) (hm : m%2=1)
    (h : 2^n*m ∣ 2^i*x) : m ∣ x :=
  (odd_cancel_pow m x i hm).mp (Int.dvd_trans (Int.dvd_mul_left _ _) h)

theorem clock_bound (m : Int) (hm : m%2=1) (n j : Nat) (p q : Pair)
    (hp : validPair (2^n*m) p) (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q (j+1)) :
    j ≤ n := by
  have h := (first_return_clock _ p hp hr.1 j).mp ⟨q,hr⟩
  by_cases hle : j ≤ n
  · exact hle
  · apply False.elim
    apply h.2 n (by omega)
    obtain ⟨hit,_⟩ := h
    unfold wall at *
    rw [← Int.mul_add, ← Int.mul_sub] at hit ⊢
    simp only [Int.mul_dvd_mul_iff_left (Int.ne_of_gt (Int.pow_pos (by decide : (0:Int)<2)))]
    rcases hit with h|h|h|h
    · exact Or.inl (odd_core_divides m _ n j hm h)
    · exact Or.inr (Or.inl (odd_core_divides m _ n j hm h))
    · exact Or.inr (Or.inr (Or.inl (odd_core_divides m _ n j hm h)))
    · exact Or.inr (Or.inr (Or.inr (odd_core_divides m _ n j hm h)))

/-- For even zero-coordinate parameters prime to the odd core in the weak
sense m ∤ u, the only possible return depth is n. -/
theorem zero_critical_clock (m u : Int) (hm : m%2=1) (n : Nat)
    (hu : u%2=0) (hnon : ¬m ∣ u) (hp : validPair (2^n*m) (u,0)) :
    (∃ q, FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) (u,0) q (n+1)) ↔
    m ∣ u-1 ∨ m ∣ 2*u-1 ∨ m ∣ 1 := by
  rw [first_return_clock _ _ hp (Or.inl rfl), zero_wall]
  have cancel (x : Int) : 2^n*m ∣ 2^n*x ↔ m ∣ x :=
    Int.mul_dvd_mul_iff_left (Int.ne_of_gt (Int.pow_pos (by decide : (0:Int)<2)))
  have last : 2^n*m ∣ 2^n ↔ m ∣ 1 := by simpa using cancel 1
  simp only [cancel, last, hnon, false_or]
  constructor
  · exact And.left
  · intro h
    refine ⟨h, fun i hi => ?_⟩
    rw [zero_wall]
    have hone := odd_no_early m 1 n i hi (by decide)
    simp only [Int.mul_one] at hone
    exact fun h => h.elim (odd_no_early m _ n i hi (by omega)) (fun h =>
      h.elim (fun hu => hnon (odd_core_divides m _ n i hm hu)) (fun h =>
        h.elim (odd_no_early m _ n i hi (by omega)) hone))

private theorem tri_scale (s m x : Int) (hs : 0 < s) :
    triangle (s*m) (s*x) = s*triangle m x := by
  simp only [triangle, show 2*(s*m) = s*(2*m) by grind, Int.mul_emod_mul_of_pos _ _ hs,
    ← Int.mul_sub]
  unfold absInt
  by_cases h : x % (2*m) - m < 0
  · have hh := Int.mul_neg_of_pos_of_neg hs h
    simp [h, hh, Int.mul_neg]
  · have hh := Int.mul_nonneg (show 0 ≤ s by omega) (show 0 ≤ x%(2*m)-m by omega)
    simp [h, show ¬ s*(x%(2*m)-m)<0 by omega]

theorem triangle_at_center (m k e : Int) (hl : -m ≤ e) (hh : e < m) :
    triangle m (m*(2*k+1)+e) = absInt e := by
  simp only [triangle, show m*(2*k+1)+e = (m+e)+k*(2*m) by grind,
    Int.add_mul_emod_self_right,
    Int.emod_eq_of_lt (a:=m+e) (b:=2*m) (by omega) (by omega)]
  congr 1; omega

theorem sortPair_map (f : Int → Int) (x y : Int) :
    sortPair (f (sortPair x y).1) (f (sortPair x y).2) = sortPair (f x) (f y) := by
  grind [sortPair]

private theorem zero_target (m u : Int) (n : Nat) (hn : 0 < n) :
    returnTarget (2^n*m) (u,0) n =
    sortPair (2^n*triangle m (u-1)) (2^n*triangle m (2^n*m-u)) := by
  cases n with
  | zero => omega
  | succ n =>
    simp only [returnTarget, pairMap, if_true, tri_scale _ _ _ (Int.pow_pos (by decide : (0:Int)<2))]
    exact sortPair_map (fun x => 2^(n+1)*triangle m x) _ _

private theorem exists_odd_factor (m x : Int) (hm : m%2=1) (hx : x%2=1) (hd : m∣x) :
    ∃ k, x = m*(2*k+1) := by
  obtain ⟨w,hw⟩ := hd
  have hh := congrArg (· % 2) hw
  simp only [Int.mul_emod m w, hm, Int.one_mul, Int.emod_emod] at hh
  exact ⟨w/2, by rw [hw, show 2*(w/2)+1=w by omega]⟩

theorem zero_residue_one_return (m u : Int) (hm : m%2=1) (hm3 : 3 ≤ m) (n : Nat)
    (hn : 0 < n) (hp : validPair (2^n*m) (u,0)) (hu : u%2=0) (hd : m ∣ u-1) :
    FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) (u,0) (2^n,0) (n+1) := by
  have hnon : ¬m ∣ u := by
    intro h
    have h1 := Int.dvd_sub h hd
    rw [show u-(u-1) = 1 by omega] at h1
    have := Int.le_of_dvd (by decide : (0:Int)<1) h1
    omega
  obtain ⟨q,hq⟩ := (zero_critical_clock m u hm n hu hnon hp).mpr (Or.inl hd)
  have target := ((first_return_iff _ _ _ hp (Or.inl rfl) n).mp hq).1
  have hs : (0:Int)<2^n := Int.pow_pos (by decide)
  obtain ⟨t,he⟩ := exists_odd_factor m (u-1) hm (by omega) hd
  obtain ⟨k, sk⟩ : ∃ k : Int, 2^n=2*k := by
    cases n with
    | zero => omega
    | succ n => exact ⟨2^n, by rw [Int.pow_succ]; grind⟩
  have ea : triangle m (u-1) = 0 := by
    rw [he, ← Int.add_zero (m*(2*t+1)), triangle_at_center _ _ _ (by omega) (by omega)]
    decide
  have eb : triangle m (2^n*m-u) = 1 := by
    have id : 2^n*m-u = m*(2*(k-t-1)+1)+(-1) := by rw [sk]; grind
    rw [id, triangle_at_center _ _ _ (by omega) (by omega)]
    decide
  rw [zero_target _ _ n hn, ea, eb] at target
  simp [sortPair, show ¬(0:Int) ≥ 2^n by omega] at target
  rwa [target] at hq

/-- The first of P1's two exceptional cycles, as a digit orbit with its least
positive return time. -/
theorem residue_one_digit_period (m : Int) (hm : m%2=1) (hm3 : 3 ≤ m) (n : Nat)
    (hn : 0<n) (hres : m ∣ 2^n-1) :
    let B := 2^n*m
    let q := fromPair B (2^n,0)
    iter (step B) (n+1) q = q ∧ ∀ i, 0<i → i<n+1 → iter (step B) i q ≠ q := by
  have hs : (0:Int)<2^n := Int.pow_pos (by decide)
  have hp : validPair (2^n*m) (2^n,0) := by
    have hh := Int.mul_lt_mul_of_pos_left (show 1<m by omega) hs
    simp only [Int.mul_one] at hh
    exact ⟨by omega, hh, Int.le_refl 0, by omega⟩
  have ev : (2:Int)^n%2=0 := by
    cases n with
    | zero => omega
    | succ n => simp [Int.pow_succ]
  have hr := zero_residue_one_return m (2^n) hm hm3 n hn hp ev hres
  have per := first_return_self_exact_period _ _ _ _ hr
  exact ⟨(digit_return_iff _ _ hp _).mpr per.2.1,
    fun i hi hj h => per.2.2 i hi hj ((digit_return_iff _ _ hp i).mp h)⟩

theorem diagonal_critical_clock (m u v : Int) (hm : m%2=1) (n : Nat)
    (hv : v≠0) (hb : u=v ∨ u+v=2^n*m) (hnon : ¬m ∣ 4*u) (h2 : ¬m ∣ 2)
    (hp : validPair (2^n*m) (u,v)) :
    (∃ q, FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) (u,v) q (n+1)) ↔
    m ∣ 2*u+1 ∨ m ∣ 2*u-1 := by
  rw [first_return_clock _ _ hp (Or.inr hb), diagonal_wall _ _ _ _ hv hb]
  have cancel (x : Int) : 2^n*m ∣ 2^n*x ↔ m ∣ x :=
    Int.mul_dvd_mul_iff_left (Int.ne_of_gt (Int.pow_pos (by decide : (0:Int)<2)))
  simp only [Int.mul_comm 2 (2^n), cancel, hnon, h2, or_false]
  constructor
  · exact And.left
  · intro h
    refine ⟨h, fun i hi => ?_⟩
    rw [diagonal_wall _ _ _ _ hv hb, Int.mul_comm 2 (2^i)]
    exact fun h => h.elim (odd_no_early m _ n i hi (by omega)) (fun h =>
      h.elim (odd_no_early m _ n i hi (by omega)) (fun h =>
        h.elim (fun hu => hnon (odd_core_divides m _ n i hm hu))
          (fun hu => h2 (odd_core_divides m _ n i hm hu))))

private theorem tri_neg (m x : Int) (hm : 0<m) : triangle m (-x) = triangle m x := by
  simp only [triangle, Int.neg_emod]
  split
  · rename_i h
    rw [Int.emod_eq_zero_of_dvd h]
  · have hh : ((2*m).natAbs : Int) = 2*m := Int.natAbs_of_nonneg (by omega)
    rw [hh]
    grind [absInt]

private theorem tri_abs (m x : Int) (hm : 0<m) : triangle m (absInt x) = triangle m x := by
  unfold absInt
  split
  · exact tri_neg _ _ hm
  · rfl

private theorem zero_target_mod (m u : Int) (n : Nat) (hn : 0<n) (hm : 0<m) :
    returnTarget (2^n*m) (u,0) n =
      sortPair (2^n*triangle m (u-1)) (2^n*triangle m u) := by
  rw [zero_target m u n hn]
  have eq : triangle m (2^n*m-u) = triangle m (-u) := by
    cases n with
    | zero => omega
    | succ k =>
      simp only [triangle, Int.pow_succ, show 2^k*2*m-u = (-u)+2^k*(2*m) by grind,
        Int.add_mul_emod_self_right]
  rw [eq, tri_neg m u hm]

private theorem diag_target_mod (m u v : Int) (n : Nat) (hn : 0<n) (hm : 0<m)
    (hv : v≠0) (hb : u=v ∨ u+v=2^n*m) :
    returnTarget (2^n*m) (u,v) n =
      sortPair (2^n*triangle m (2*u+1)) (2^n*triangle m (2*u-1)) := by
  cases n with
  | zero => omega
  | succ k =>
    simp only [returnTarget, pairMap, if_neg hv, if_pos hb,
      tri_scale _ _ _ (Int.pow_pos (by decide : (0:Int)<2))]
    rw [sortPair_map (fun x => 2^(k+1)*triangle m x), tri_abs _ _ hm, tri_abs _ _ hm]
    have eq (a : Int) : triangle m (2*u-2^(k+1)*m+a) = triangle m (2*u+a) := by
      simp only [triangle, Int.pow_succ, show 2*u-2^k*2*m+a=(2*u+a)+(-2^k)*(2*m) by grind,
        Int.add_mul_emod_self_right]
    rw [eq 1, show 2*u-2^(k+1)*m-1=2*u-2^(k+1)*m+(-1) by omega, eq (-1)]
    rfl

theorem zero_critical_return (m u : Int) (hm : m%2=1) (hm3 : 3≤m) (n : Nat)
    (hn : 0<n) (hu : u%2=0) (hnon : ¬m ∣ u) (hp : validPair (2^n*m) (u,0))
    (hd : m ∣ u-1 ∨ m ∣ 2*u-1) :
    FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) (u,0)
      (sortPair (2^n*triangle m (u-1)) (2^n*triangle m u)) (n+1) := by
  obtain ⟨q,hq⟩ := (zero_critical_clock m u hm n hu hnon hp).mpr (by grind)
  have eq := ((first_return_iff _ _ _ hp (Or.inl rfl) n).mp hq).1
  rwa [eq, zero_target_mod m u n hn (by omega)] at hq

theorem diag_critical_return (m u v : Int) (hm : m%2=1) (hm3 : 3≤m) (n : Nat)
    (hn : 0<n) (hv : v≠0) (hb : u=v ∨ u+v=2^n*m) (hnon : ¬m ∣ 4*u) (hp : validPair (2^n*m) (u,v))
    (hd : m ∣ 2*u+1 ∨ m ∣ 2*u-1) :
    FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) (u,v)
      (sortPair (2^n*triangle m (2*u+1)) (2^n*triangle m (2*u-1))) (n+1) := by
  have h2 : ¬m∣2 := fun h => by have := Int.le_of_dvd (by decide : (0:Int)<2) h; omega
  obtain ⟨q,hq⟩ := (diagonal_critical_clock m u v hm n hv hb hnon h2 hp).mpr hd
  have eq := ((first_return_iff _ _ _ hp (Or.inr hb) n).mp hq).1
  rwa [eq, diag_target_mod m u v n hn (by omega) hv hb] at hq

#print axioms residue_one_digit_period
#print axioms clock_bound
#print axioms zero_critical_clock
end KaprekarBoundary
