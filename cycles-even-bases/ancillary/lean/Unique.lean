import CoreCounts

namespace KaprekarBoundary

/-- Exactly one nonconstant cycle: there is a periodic word and its forward
orbit contains every nonconstant periodic word. -/
def uniqueDigitCycle (B : Int) : Prop := ∃ q : Digits, q.valid B ∧ q.nonconstant ∧
  periodic (step B) q ∧ ∀ z : Digits, z.valid B → z.nonconstant → periodic (step B) z → reaches (step B) q z

def uniquePairCycle (B : Int) : Prop := ∃ p, validPair B p ∧ periodic (pairMap B) p ∧
  ∀ z, validPair B z → periodic (pairMap B) z → reaches (pairMap B) p z

private theorem word_ne (B : Int) (p : Pair) (hp : validPair B p) : (diff (fromPair B p)).1≠0 := by
  rw [pair_of_subtraction B p hp]
  have := (pairMap_valid B p hp).1
  omega

private theorem lift_reach (B : Int) (p q : Pair) (hp : validPair B p) (h : reaches (pairMap B) p q) :
    reaches (step B) (fromPair B p) (fromPair B q) := by
  obtain ⟨k,hk⟩ := h
  exact ⟨k,by rw [iter_fromPair B p hp,hk]⟩

theorem unique_digit_iff_pair (B : Int) : uniqueDigitCycle B ↔ uniquePairCycle B := by
  constructor
  · rintro ⟨q,hV,hne,hqp,hall⟩
    obtain ⟨p,hp,rfl,hper⟩ := (word_periodic_iff B q hV ((nonconstant_iff q).mp hne)).mp hqp
    refine ⟨p,hp,hper,?_⟩
    intro z hz zper
    have wper := (word_periodic_iff B (fromPair B z) (fromPair_valid B z hz) (word_ne B z hz)).mpr
      ⟨z,hz,rfl,zper⟩
    obtain ⟨i,hi⟩ := hall _ (fromPair_valid B z hz) ((nonconstant_iff _).mpr (word_ne B z hz)) wper
    rw [iter_fromPair B p hp i] at hi
    exact ⟨i,fromPair_injective _ _ _ (iter_pair_valid B p hp i) hz hi⟩
  · rintro ⟨p,hp,hper,hall⟩
    have hw := (word_periodic_iff B (fromPair B p) (fromPair_valid B p hp) (word_ne B p hp)).mpr
      ⟨p,hp,rfl,hper⟩
    refine ⟨_,fromPair_valid B p hp,(nonconstant_iff _).mpr (word_ne B p hp),hw,?_⟩
    intro q hq hn hper
    obtain ⟨t,ht,rfl,kt⟩ := (word_periodic_iff B q hq ((nonconstant_iff q).mp hn)).mp hper
    exact lift_reach _ _ _ hp (hall t ht kt)

private theorem joined (B : Int) (hu : uniquePairCycle B) (x y : Pair)
    (vx : validPair B x) (vy : validPair B y) (hx : periodic (pairMap B) x) (hy : periodic (pairMap B) y) :
    reaches (pairMap B) x y := by
  obtain ⟨p,_,hp,hall⟩ := hu
  exact reaches_trans _ x p y (periodic_reachable _ p x hp (hall x vx hx)).2 (hall y vy hy)

private theorem return_on_reach {α : Type} (f : α → α) (p q : α) (k : Nat) (hp : iter f k p=p)
    (h : reaches f p q) : iter f k q=q := by
  obtain ⟨i,rfl⟩ := h
  rw [iter_commute,hp]

private theorem pure_not_unique (n : Nat) (hn : 0<n) : ¬uniquePairCycle (2^n) := by
  intro hu
  by_cases he : n=1
  · subst n
    have f0 : pairMap 2 (1,0)=(1,0) := by decide
    have f1 : pairMap 2 (1,1)=(1,1) := by decide
    have vals : validPair 2 (1,0) ∧ validPair 2 (1,1) := by dsimp [validPair];omega
    have h := joined 2 (by simpa using hu) _ _ vals.1 vals.2 ⟨0,f0⟩ ⟨0,f1⟩
    obtain ⟨k,hk⟩ := h
    have fixed (i : Nat) : iter (pairMap 2) i (1,0)=(1,0) := by
      induction i with
      | zero => rfl
      | succ i ih => simp only [iter,ih,f0]
    rw [fixed k] at hk; contradiction
  · have ex : n=0+(n-2)+2 := by omega
    have re := power_rep_returns 0 (n-2)
    rw [← ex] at re
    have h0 : validPair (2^n) (powZ 0) ∧ validPair (2^n) (powA n 0) := by
      have h := pow_valid 0 (n-2)
      rwa [← ex] at h
    have per (p : Pair) (v : validPair (2^n) p) (e : p=powZ 0 ∨ p=powA n 0) : periodic (pairMap (2^n)) p := by
      apply (power_pair_classification n (by omega) p v).mpr
      refine ⟨0,n-2,ex,by omega,?_⟩
      rcases e with rfl|rfl
      · exact Or.inl ⟨0,rfl⟩
      · exact Or.inr ⟨0,rfl⟩
    have hit := joined _ hu _ _ h0.1 h0.2 (per _ h0.1 (Or.inl rfl)) (per _ h0.2 (Or.inr rfl))
    have bad := (power_reps_disjoint n 0 (n-2) 0 (n-2) ex ex (by omega) (by omega) false true).mp hit
    simp at bad

private theorem three_unique (n : Nat) (hn : 0<n) :
    uniquePairCycle (2^n*3) ↔ n%2=1 := by
  let s : Int := 2^n
  let B := s*3
  have hs : 0<s := Int.pow_pos (by decide)
  have vS : validPair B (s,0) := by dsimp [validPair,B];omega
  have vE : validPair B (2*s,0) := by dsimp [validPair,B];omega
  have perS := (three_pair_classification n hn (s,0) vS).mpr (Or.inl ⟨0,rfl⟩)
  constructor
  · intro hu
    by_cases h : n%2=0
    · have perE := (three_pair_classification n hn (2*s,0) vE).mpr (Or.inr ⟨h,0,rfl⟩)
      have eq := three_digit_periods n hn
      simp only [if_pos h] at eq
      have hh := return_on_reach _ _ _ _ ((digit_return_iff _ _ vS _).mp eq.1.2.1) (joined B hu _ _ vS vE perS perE)
      exact False.elim (eq.2.2.2 (n+1) (by omega) (by omega) ((digit_return_iff _ _ vE _).mpr hh))
    · omega
  · intro h
    refine ⟨_,vS,perS,?_⟩
    intro z vz hp
    exact ((three_pair_classification n hn z vz).mp hp).resolve_right (fun hh => by have := hh.1;omega)

private theorem large_not_unique (m : Int) (hm : m%2=1) (h9 : 9≤m) (n : Nat) :
    ¬uniquePairCycle (2^n*m) := by
  have c1 : oddCore m (m-2,1) := by dsimp [oddCore,validPair];omega
  have c2 : oddCore m (3,1) := by dsimp [oddCore,validPair];omega
  have v1 := scaled_core_valid m (2^n) (Int.pow_pos (by decide)) _ c1
  have v2 := scaled_core_valid m (2^n) (Int.pow_pos (by decide)) _ c2
  obtain ⟨p1,p2,sep⟩ := two_large_core_cycles m hm h9 n
  exact fun hu => sep (joined _ hu _ _ v1 v2 p1 p2)

private theorem seven_power (n : Nat) : (2:Int)^n%7 = (if n%3=0 then 1 else if n%3=1 then 2 else 4) := by
  induction n with
  | zero => decide
  | succ n ih =>
    rw [Int.pow_succ]
    repeat' split at ih
    all_goals repeat' split
    all_goals omega

private theorem small_core_reach (m : Int) (hm : m=5 ∨ m=7) (p : Pair) (hp : oddCore m p) :
    reaches (pairMap m) (3,1) p := by
  rcases p with ⟨a,b⟩
  rcases hm with rfl|rfl
  · have ha : a=3 ∧ b=1 := by grind [oddCore,validPair]
    exact ⟨0,by rw [ha.1,ha.2];rfl⟩
  · have ha : (a=3 ∧ b=1) ∨ (a=5 ∧ b=1) ∨ (a=5 ∧ b=3) := by grind [oddCore,validPair]
    rcases ha with ⟨rfl,rfl⟩|⟨rfl,rfl⟩|⟨rfl,rfl⟩
    · exact ⟨0,rfl⟩
    · exact ⟨1,by decide⟩
    · exact ⟨2,by decide⟩

private theorem small_large_unique (m : Int) (hm : m=5 ∨ m=7) (n : Nat) (hn : 0<n) :
    uniquePairCycle (2^n*m) ↔ (¬m∣2^n-1 ∧ ¬(m=5 ∧ n%4=2)) := by
  let s : Int := 2^n
  let B := s*m
  let f := pairMap B
  have hmO : m%2=1 := by rcases hm with rfl|rfl <;> decide
  have hm5 : 5≤m := by omega
  have ht : oddCore m (3,1) := by rcases hm with rfl|rfl <;> dsimp [oddCore,validPair] <;> omega
  have hs : 0<s := Int.pow_pos (by decide)
  have hu := scaled_core_valid m s hs (3,1) ht
  have per := scaled_periodic m hmO n _ ht
  have vS : validPair B (s,0) := by
    have hh := Int.mul_lt_mul_of_pos_left (show 1<m by omega) hs
    exact ⟨by omega,by simpa [B] using hh,Int.le_refl 0,by omega⟩
  have vE : validPair B (2*s,0) := by
    have hh := Int.mul_lt_mul_of_pos_left (show 2<m by omega) hs
    exact ⟨by omega,by dsimp [B];grind,Int.le_refl 0,by omega⟩
  constructor
  · intro h
    have no (p : Pair) (vp : validPair B p) (bp : boundary B p) : ¬periodic f p := by
      intro kp
      obtain ⟨i,hi⟩ := joined B h _ p hu vp per kp
      exact scaled_avoids_boundary m hmO n _ ht i (by rw [hi];exact bp)
    constructor
    · intro hd
      exact no _ vS (Or.inl rfl) ((large_pair_classification m hmO hm5 n hn _ vS).mpr (Or.inr (Or.inl ⟨hd,0,rfl⟩)))
    · rintro ⟨m5,n2⟩
      exact no _ vE (Or.inl rfl) ((large_pair_classification m hmO hm5 n hn _ vE).mpr (Or.inr (Or.inr ⟨m5,n2,0,rfl⟩)))
  · rintro ⟨hS,hE⟩
    refine ⟨_,hu,per,?_⟩
    intro z vz zp
    rcases (large_pair_classification m hmO hm5 n hn z vz).mp zp with ⟨t,rfl,ct,_⟩|⟨h,_⟩|⟨m5,n2,_⟩
    · obtain ⟨i,hi⟩ := small_core_reach m hm t ct
      exact ⟨i,by rw [iter_scaled_core m hmO n _ ht,hi]⟩
    · exact False.elim (hS h)
    · exact False.elim (hE ⟨m5,n2⟩)

/-- P1's unique-cycle corollary. The odd core and positive exponent are the
canonical factorization of an even base. -/
theorem unique_cycle (m : Int) (hm : m%2=1) (hpos : 0<m) (n : Nat) (hn : 0<n) :
    uniqueDigitCycle (2^n*m) ↔ (m=3 ∧ n%2=1) ∨ (m=5 ∧ n%2=1) ∨ (m=7 ∧ n%3≠0) := by
  rw [unique_digit_iff_pair]
  rcases (show m=1 ∨ m=3 ∨ m=5 ∨ m=7 ∨ 9≤m by omega) with rfl|rfl|rfl|rfl|h9
  · simp [pure_not_unique n hn]
  · simpa using three_unique n hn
  · rw [small_large_unique 5 (Or.inl rfl) n hn]
    have h := five_power_mod n
    repeat' split at h
    all_goals omega
  · rw [small_large_unique 7 (Or.inr rfl) n hn]
    have h := seven_power n
    repeat' split at h
    all_goals omega
  · have h := large_not_unique m hm h9 n
    grind

private theorem scale_period (m : Int) (hm : m%2=1) (n : Nat) (t : Pair) (ht : oddCore m t) (k : Nat)
    (hk : 0<k) (h : iter (pairMap m) k t=t ∧ ∀ i, 0<i → i<k → iter (pairMap m) i t≠t) :
    digitPeriod (2^n*m) (scale (2^n) t) k := by
  have hs : (0:Int)<2^n := Int.pow_pos (by decide)
  apply digitPeriod_of_pair _ _ (scaled_core_valid m _ hs t ht) _ hk
  refine ⟨by rw [iter_scaled_core m hm n t ht,h.1],?_⟩
  intro i hi hj he
  rw [iter_scaled_core m hm n t ht] at he
  apply h.2 i hi hj
  apply Prod.ext <;> apply (Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)).mp
  · exact congrArg Prod.fst he
  · exact congrArg Prod.snd he

/-- Lengths of the unique digit orbit in the three congruence families. -/
theorem unique_cycle_periods (n : Nat) (hn : 0<n) :
    (n%2=1 → digitPeriod (2^n*3) (2^n,0) (6*(n/2+1))) ∧
    (n%2=1 → digitPeriod (2^n*5) (scale (2^n) (3,1)) 1) ∧
    (n%3≠0 → digitPeriod (2^n*7) (scale (2^n) (3,1)) 3) := by
  constructor
  · intro h
    have per := three_digit_periods n hn
    simpa only [if_neg (show ¬n%2=0 by omega),show 3*(n+1)=6*(n/2+1) by omega] using per
  constructor
  · intro _
    apply scale_period 5 (by decide) n (3,1) (by dsimp [oddCore,validPair];omega) 1 (by decide)
    exact ⟨by decide,fun i hi hj => by omega⟩
  · intro _
    apply scale_period 7 (by decide) n (3,1) (by dsimp [oddCore,validPair];omega) 3 (by decide)
    refine ⟨by decide,?_⟩
    intro i hi hj
    rcases (show i=1 ∨ i=2 by omega) with rfl|rfl <;> decide

#print axioms unique_cycle
#print axioms unique_cycle_periods
end KaprekarBoundary
