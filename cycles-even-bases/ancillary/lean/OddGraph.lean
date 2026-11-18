import Cycles

namespace KaprekarBoundary

private theorem no_dvd_one (m : Int) (h : 3≤m) : ¬m∣1 :=
  fun hd => by have := Int.le_of_dvd (by decide : (0:Int)<1) hd; omega

private theorem no_dvd_power (m : Int) (hm : m%2=1) (h : 3≤m) (n : Nat) : ¬m∣2^n :=
  fun hd => no_dvd_one m h ((odd_cancel_pow m 1 n hm).mp (by simpa using hd))

private theorem tri_zero_dvd (m x : Int) (h : triangle m x=0) : m∣x := by
  have he : x%(2*m)=m := by unfold triangle absInt at h; split at h <;> omega
  have id := Int.emod_add_mul_ediv x (2*m)
  exact ⟨1+2*(x/(2*m)), by grind⟩

private theorem zero_successor (m : Int) (hm : m%2=1) (h3 : 3≤m) (n : Nat) (hn : 0<n) (u : Int)
    (hev : u%2=0) (hnd : ¬m∣u) (hv : validPair (2^n*m) (u,0))
    (hp : periodic (pairMap (2^n*m)) (u,0)) :
    let s : Int := 2^n
    (m∣u-1 ∧ reaches (pairMap (s*m)) (u,0) (s,0) ∧ periodic (pairMap (s*m)) (s,0)) ∨
    (m∣2*u-1 ∧ ∃ z, z.1+z.2=s*m ∧ z.1-z.2=s ∧
      reaches (pairMap (s*m)) (u,0) z ∧ periodic (pairMap (s*m)) z) := by
  let s : Int := 2^n
  let B := s*m
  let f := pairMap B
  have hs : 0<s := Int.pow_pos (by decide)
  obtain ⟨q,j,hr,_,_,_⟩ := periodic_first_return (pairMap B) (boundary B) (u,0) hp (Or.inl rfl)
  have tests := ((first_return_clock _ _ hv hr.1 j).mp ⟨q,hr⟩).1
  rw [zero_wall] at tests
  by_cases h : m∣u-1
  · have ret := zero_residue_one_return m u hm h3 n hn hv hev h
    have reach : reaches f (u,0) (s,0) := ⟨n+1,ret.2.2.2.1⟩
    exact Or.inl ⟨h,reach,(periodic_reachable _ _ _ hp reach).1⟩
  · have one : ¬m∣1 := no_dvd_one m h3
    have cond : m∣2*u-1 := by
      rcases tests with h0|h1|h2|h3
      · exact False.elim (h (odd_core_divides m _ n j hm h0))
      · exact False.elim (hnd (odd_core_divides m _ n j hm h1))
      · exact odd_core_divides m _ n j hm h2
      · exact False.elim (one (odd_core_divides m 1 n j hm (by simpa using h3)))
    have ret := zero_critical_return m u hm h3 n hn hev hnd hv (Or.inr cond)
    let a := triangle m (u-1)
    let b := triangle m u
    let z := sortPair (s*a) (s*b)
    have neA : a≠0 := fun h0 => h (tri_zero_dvd _ _ h0)
    have neB : b≠0 := fun h0 => hnd (tri_zero_dvd _ _ h0)
    have nz : z.2≠0 := by simp only [z,sortPair];split <;> simp [Int.mul_eq_zero,neA,neB,Int.ne_of_gt hs]
    have shape := full_depth_target m hm n _ z hn hv ret
    have reach : reaches f (u,0) z := ⟨n+1,ret.2.2.2.1⟩
    obtain ⟨sm,df⟩ : z.1+z.2=B ∧ z.1-z.2=s := shape.resolve_left (by grind) |>.resolve_left (by grind)
    exact Or.inr ⟨cond,z,sm,df,reach,(periodic_reachable _ _ _ hp reach).1⟩

private theorem tri_neighbors (m x : Int) (hm : m%2=1) (hg : 3≤m) (hx : x%2=1) (hd : m∣x) :
    triangle m x=0 ∧ triangle m (x+2)=2 ∧ triangle m (x-2)=2 := by
  obtain ⟨w,hw⟩ := hd
  have h := congrArg (·%2) hw
  simp only [Int.mul_emod m w, hm,Int.one_mul,Int.emod_emod] at h
  have he : x=m*(2*(w/2)+1) := by rw [hw,show 2*(w/2)+1=w by omega]
  rw [he]
  constructor
  · rw [← Int.add_zero (m*(2*(w/2)+1)),triangle_at_center _ _ _ (by omega) (by omega)];decide
  constructor
  · rw [triangle_at_center _ _ _ (by omega) (by omega)];decide
  · rw [show m*(2*(w/2)+1)-2=m*(2*(w/2)+1)+(-2) by omega,triangle_at_center _ _ _ (by omega) (by omega)];decide

private theorem anti_successor (m : Int) (hm : m%2=1) (h3 : 3≤m) (n : Nat) (hn : 0<n) (p : Pair)
    (hV : validPair (2^n*m) p) (ha : p.1+p.2=2^n*m) (hg : p.1-p.2=2^n)
    (hp : periodic (pairMap (2^n*m)) p) :
    let s : Int := 2^n
    (m∣s+1 ∨ m∣s-1) ∧ reaches (pairMap (s*m)) p (2*s,0) ∧ periodic (pairMap (s*m)) (2*s,0) := by
  let s : Int := 2^n
  let B := s*m
  have hs : 0<s := Int.pow_pos (by decide)
  change p.1+p.2=B at ha
  change p.1-p.2=s at hg
  have ph : 2*p.1 = B+s := by omega
  have nod : ¬m∣4*p.1 := by
    intro h
    have hh := (odd_cancel_two m (2*p.1) hm).mp (by simpa [show (4:Int)*p.1=2*(2*p.1) by omega] using h)
    rw [ph, Int.dvd_add_right (Int.dvd_mul_left s m)] at hh
    exact no_dvd_power m hm h3 n hh
  have hv : p.2≠0 := by have := hV.2.1; change p.1<B at this; omega
  obtain ⟨q,j,hr,_,_,_⟩ := periodic_first_return (pairMap B) (boundary B) p hp (Or.inr (Or.inr ha))
  have tests := ((first_return_clock _ _ hV hr.1 j).mp ⟨q,hr⟩).1
  rw [show p=(p.1,p.2) from rfl, diagonal_wall _ _ _ _ hv (Or.inr ha)] at tests
  have cond : m∣2*p.1+1 ∨ m∣2*p.1-1 := by
    rcases tests with h|h|h|h
    · exact Or.inl (odd_core_divides m _ n j hm h)
    · exact Or.inr (odd_core_divides m _ n j hm h)
    · exact False.elim (nod (odd_core_divides m _ n j hm h))
    · rw [Int.mul_comm] at h
      have := Int.le_of_dvd (by decide : (0:Int)<2) (odd_core_divides m 2 n j hm h)
      omega
  have ret := diag_critical_return m p.1 p.2 hm h3 n hn hv (Or.inr ha) nod hV cond
  have tar : sortPair (s*triangle m (2*p.1+1)) (s*triangle m (2*p.1-1))=(2*s,0) := by
    rcases cond with h|h
    · have tri := tri_neighbors m (2*p.1+1) hm h3 (by omega) h
      rw [show 2*p.1+1-2=2*p.1-1 by omega] at tri
      rw [tri.1,tri.2.2]
      grind [sortPair]
    · have tri := tri_neighbors m (2*p.1-1) hm h3 (by omega) h
      rw [show 2*p.1-1+2=2*p.1+1 by omega] at tri
      rw [tri.1,tri.2.1]
      grind [sortPair]
  have reach : reaches (pairMap B) p (2*s,0) := ⟨n+1, ret.2.2.2.1.trans tar⟩
  refine ⟨?_,reach,(periodic_reachable _ _ _ hp reach).1⟩
  rw [ph,show B+s+1=B+(s+1) by omega,show B+s-1=B+(s-1) by omega,
    Int.dvd_add_right (Int.dvd_mul_left s m),Int.dvd_add_right (Int.dvd_mul_left s m)] at cond
  exact cond

private theorem valid_reached (B : Int) (p q : Pair) (hp : validPair B p) (h : reaches (pairMap B) p q) :
    validPair B q := by obtain ⟨i,rfl⟩ := h; exact iter_pair_valid _ _ hp i

private theorem div_combine (m x y a b : Int) (hx : m∣x) (hy : m∣y) : m∣a*x+b*y :=
  Int.dvd_add (Int.dvd_mul_of_dvd_right hx) (Int.dvd_mul_of_dvd_right hy)

private theorem large_S_periodic (m : Int) (hm : m%2=1) (h5 : 5≤m) (n : Nat) (hn : 0<n)
    (hper : periodic (pairMap (2^n*m)) (2^n,0)) : m∣2^n-1 := by
  let s : Int := 2^n
  let B := s*m
  have hs : 0<s := Int.pow_pos (by decide)
  have hV : validPair B (s,0) := by
    have hh := Int.mul_lt_mul_of_pos_left (show 1<m by omega) hs
    exact ⟨by omega,by simpa [B] using hh,Int.le_refl 0,by omega⟩
  have hev : s%2=0 := by
    cases n with
    | zero => omega
    | succ n => simp [s,Int.pow_succ]
  rcases zero_successor m hm (by omega) n hn s hev (no_dvd_power m hm (by omega) n) hV hper with
    ⟨hd,_,_⟩|⟨hd,z,ha,hg,hreach,hz⟩
  · exact hd
  · obtain ⟨h,_,_⟩ := anti_successor m hm (by omega) n hn z (valid_reached _ _ _ hV hreach) ha hg hz
    change m∣s+1 ∨ m∣s-1 at h
    rcases h with h|h
    · have hc := div_combine m _ _ (-1) 2 hd h
      have : m∣3 := by rwa [show -(1:Int)*(2*s-1)+2*(s+1)=3 by omega] at hc
      have := Int.le_of_dvd (by decide : (0:Int)<3) this
      omega
    · exact h

theorem five_power_mod (n : Nat) : (2:Int)^n%5 =
    (if n%4=0 then 1 else if n%4=1 then 2 else if n%4=2 then 4 else 3) := by
  induction n with
  | zero => decide
  | succ n ih =>
    rw [Int.pow_succ]
    repeat' split at ih
    all_goals repeat' split
    all_goals omega

private theorem large_E_periodic (m : Int) (hm : m%2=1) (h5 : 5≤m) (n : Nat) (hn : 0<n)
    (hper : periodic (pairMap (2^n*m)) (2*2^n,0)) : m=5 ∧ n%4=2 := by
  let s : Int := 2^n
  let B := s*m
  have hs : 0<s := Int.pow_pos (by decide)
  have hV : validPair B (2*s,0) := by
    have hh := Int.mul_lt_mul_of_pos_left (show 2<m by omega) hs
    exact ⟨by omega,by dsimp [B]; grind,Int.le_refl 0,by omega⟩
  have nd : ¬m∣2*s := fun h => no_dvd_power m hm (by omega) n ((odd_cancel_two _ _ hm).mp h)
  rcases zero_successor m hm (by omega) n hn (2*s) (by omega) nd hV hper with
    ⟨hd,_,hS⟩|⟨hd,z,ha,hg,hreach,hz⟩
  · have h := large_S_periodic m hm h5 n hn hS
    have hc := div_combine m _ _ 1 (-2) hd h
    change m∣1*(2*s-1)+(-2)*(s-1) at hc
    rw [show 1*(2*s-1)+(-2)*(s-1)=1 by omega] at hc
    exact False.elim (no_dvd_one _ (by omega) hc)
  · obtain ⟨h,_,_⟩ := anti_successor m hm (by omega) n hn z (valid_reached _ _ _ hV hreach) ha hg hz
    change m∣s+1 ∨ m∣s-1 at h
    rcases h with h|h
    · have hc := div_combine m _ _ (-1) 4 hd h
      rw [show -(1:Int)*(2*(2*s)-1)+4*(s+1)=5 by omega] at hc
      have hm5 : m=5 := by have := Int.le_of_dvd (by decide : (0:Int)<5) hc; omega
      subst m
      have hh := five_power_mod n
      change s%5=_ at hh
      refine ⟨rfl,?_⟩
      repeat' split at hh
      all_goals omega
    · have hc := div_combine m _ _ 1 (-4) hd h
      rw [show 1*(2*(2*s)-1)+(-4)*(s-1)=3 by omega] at hc
      have := Int.le_of_dvd (by decide : (0:Int)<3) hc
      omega

/-- Exhaustive classification on pairs for every odd core at least five. -/
theorem large_pair_classification (m : Int) (hm : m%2=1) (h5 : 5≤m) (n : Nat) (hn : 0<n)
    (p : Pair) (hp : validPair (2^n*m) p) :
    let s : Int := 2^n
    let f := pairMap (s*m)
    periodic f p ↔ (∃ t, p=scale s t ∧ oddCore m t ∧ periodic (pairMap m) t) ∨
      (m∣s-1 ∧ reaches f (s,0) p) ∨ (m=5 ∧ n%4=2 ∧ reaches f (2*s,0) p) := by
  let s : Int := 2^n
  let B := s*m
  let f := pairMap B
  have hs : 0<s := Int.pow_pos (by decide)
  have vS : validPair B (s,0) := by
    have h := Int.mul_lt_mul_of_pos_left (show 1<m by omega) hs
    exact ⟨by omega,by simpa [B] using h,Int.le_refl 0,by omega⟩
  have vE : validPair B (2*s,0) := by
    have h := Int.mul_lt_mul_of_pos_left (show 2<m by omega) hs
    exact ⟨by omega,by dsimp [B];grind,Int.le_refl 0,by omega⟩
  constructor
  · intro h
    rcases cycle_core_or_graph m hm (by omega) n hn p hp h with core|⟨z,zp,zper,zg⟩
    · exact Or.inl core
    · rcases zg with eq|eq|⟨ha,hg⟩
      · subst z; exact Or.inr (Or.inl ⟨large_S_periodic m hm h5 n hn zper,zp⟩)
      · subst z; obtain ⟨h5,he⟩ := large_E_periodic m hm h5 n hn zper
        exact Or.inr (Or.inr ⟨h5,he,zp⟩)
      · have hz := valid_reached _ _ _ hp (periodic_reachable _ z p zper zp).2
        obtain ⟨_,ze,ep⟩ := anti_successor m hm (by omega) n hn z hz ha hg zper
        obtain ⟨m5,n4⟩ := large_E_periodic m hm h5 n hn ep
        have ez := (periodic_reachable _ _ _ zper ze).2
        exact Or.inr (Or.inr ⟨m5,n4,reaches_trans _ _ _ _ ez zp⟩)
  · rintro (⟨t,eq,ht,k,hk⟩|⟨hd,hsr⟩|⟨rfl,n4,her⟩)
    · exact ⟨k,((generic_cycle_iff m hm n p hp k).mpr ⟨t,eq,ht,hk⟩).2⟩
    · have h := residue_one_digit_period m hm (by omega) n hn hd
      have per : periodic f (s,0) := ⟨n,(digit_return_iff _ _ vS _).mp h.1⟩
      exact (periodic_reachable f _ _ per hsr).1
    · have h := five_exception_digit_period n n4
      have per : periodic (pairMap (s*5)) (2*s,0) := by
        refine ⟨2*n+1,?_⟩
        rw [show 2*n+1+1=2*(n+1) by omega]
        exact (digit_return_iff _ _ vE _).mp h.1
      exact (periodic_reachable _ _ _ per her).1

theorem large_digit_classification (m : Int) (hm : m%2=1) (h5 : 5≤m) (n : Nat) (hn : 0<n)
    (q : Digits) (hq : q.valid (2^n*m)) (hne : (diff q).1≠0) :
    let s : Int := 2^n
    let B := s*m
    periodic (step B) q ↔ ∃ p, validPair B p ∧ q=fromPair B p ∧
      ((∃ t, p=scale s t ∧ oddCore m t ∧ periodic (pairMap m) t) ∨
       (m∣s-1 ∧ reaches (pairMap B) (s,0) p) ∨
       (m=5 ∧ n%4=2 ∧ reaches (pairMap B) (2*s,0) p)) := by
  dsimp only
  rw [word_periodic_iff _ _ hq hne]
  constructor
  · rintro ⟨p,hp,eq,h⟩
    exact ⟨p,hp,eq,(large_pair_classification m hm h5 n hn p hp).mp h⟩
  · rintro ⟨p,hp,eq,h⟩
    exact ⟨p,hp,eq,(large_pair_classification m hm h5 n hn p hp).mpr h⟩

#print axioms large_digit_classification
end KaprekarBoundary
