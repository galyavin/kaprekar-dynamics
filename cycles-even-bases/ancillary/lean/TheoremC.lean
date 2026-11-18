import Unique

namespace KaprekarBoundary

theorem word_reach_iff (B : Int) (p q : Pair) (hp : validPair B p) (hq : validPair B q) :
    reaches (step B) (fromPair B p) (fromPair B q) ↔ reaches (pairMap B) p q := by
  constructor
  · rintro ⟨i,he⟩
    rw [iter_fromPair _ _ hp] at he
    exact ⟨i,fromPair_injective _ _ _ (iter_pair_valid B p hp i) hq he⟩
  · rintro ⟨i,he⟩;exact ⟨i,by rw [iter_fromPair _ _ hp,he]⟩

/-- Classification for an odd core at least five, including exact periods,
and separation of each additional cycle from every lifted cycle. -/
structure LargePart (m : Int) (n : Nat) : Prop where
  words : ∀ q, q.valid (2^n*m) → q.nonconstant →
    (periodic (step (2^n*m)) q ↔
      (∃ t, oddCore m t ∧ periodic (pairMap m) t ∧ q=fromPair (2^n*m) (scale (2^n) t)) ∨
      (m ∣ 2^n-1 ∧ reaches (step (2^n*m)) (fromPair (2^n*m) (2^n,0)) q) ∨
      (m=5 ∧ n%4=2 ∧ reaches (step (2^n*m)) (fromPair (2^n*m) (2*2^n,0)) q))
  lift_clock : ∀ t, oddCore m t → ∀ i,
    (iter (step (2^n*m)) i (fromPair (2^n*m) (scale (2^n) t))=fromPair (2^n*m) (scale (2^n) t) ↔
     iter (pairMap m) i t=t)
  residue_one : m ∣ 2^n-1 → digitPeriod (2^n*m) (2^n,0) (n+1)
  exceptional_five : m=5 → n%4=2 → digitPeriod (2^n*m) (2*2^n,0) (2*(n+1))
  exclusive : ¬ (m ∣ 2^n-1 ∧ m=5 ∧ n%4=2)
  separate : ∀ t, oddCore m t → ∀ p, validPair (2^n*m) p → boundary (2^n*m) p →
    ¬reaches (step (2^n*m)) (fromPair (2^n*m) (scale (2^n) t)) (fromPair (2^n*m) p)

private theorem large_part (m : Int) (hm : m%2=1) (h5 : 5≤m) (n : Nat) (hn : 0<n) : LargePart m n := by
  let s : Int := 2^n
  let B := s*m
  have hs : 0<s := Int.pow_pos (by decide)
  have vs : validPair B (s,0) := by
    have hi := Int.mul_lt_mul_of_pos_left (show 1<m by omega) hs
    exact ⟨by omega,by simpa [B] using hi,Int.le_refl 0,by omega⟩
  have ve : validPair B (2*s,0) := by
    have hi := Int.mul_lt_mul_of_pos_left (show 2<m by omega) hs
    exact ⟨by omega,by dsimp [B];grind,Int.le_refl 0,by omega⟩
  constructor
  · intro q vq hne
    constructor
    · intro per
      replace hne := (nonconstant_iff q).mp hne
      obtain ⟨p,vp,rfl,h⟩ := (large_digit_classification m hm h5 n hn q vq hne).mp per
      rcases h with ⟨t,rfl,ht,pt⟩|⟨hr,reach⟩|⟨m5,h4,reach⟩
      · exact Or.inl ⟨t,ht,pt,rfl⟩
      · exact Or.inr (Or.inl ⟨hr,(word_reach_iff _ _ _ vs vp).mpr reach⟩)
      · exact Or.inr (Or.inr ⟨m5,h4,(word_reach_iff _ _ _ ve vp).mpr reach⟩)
    · rintro (⟨t,ht,kt,rfl⟩|⟨hr,reach⟩|⟨m5,h4,reach⟩)
      · exact (word_periodic_iff _ _ vq ((nonconstant_iff _).mp hne)).mpr
          ⟨_,scaled_core_valid m s hs t ht,rfl,scaled_periodic m hm n t ht⟩
      · have per := residue_one_digit_period m hm (by omega) n hn hr
        exact (periodic_reachable _ _ _ ⟨n,per.1⟩ reach).1
      · subst m
        have per := five_exception_digit_period n h4
        exact (periodic_reachable _ _ _ ⟨2*n+1,by rw [show 2*n+1+1=2*(n+1) by omega];exact per.1⟩ reach).1
  · intro t ht i
    rw [digit_return_iff _ _ (scaled_core_valid m s hs t ht),iter_scaled_core m hm n t ht]
    constructor
    · intro eq
      apply Prod.ext <;> apply (Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)).mp
      · exact congrArg Prod.fst eq
      · exact congrArg Prod.snd eq
    · exact congrArg (scale s)
  · intro h;exact ⟨by omega,residue_one_digit_period m hm (by omega) n hn h⟩
  · rintro rfl h;exact ⟨by omega,five_exception_digit_period n h⟩
  · rintro ⟨hd,rfl,h4⟩
    have h := five_power_mod n
    repeat' split at h
    all_goals omega
  · intro t ht p vp bp reach
    obtain ⟨i,hi⟩ := (word_reach_iff B _ p (scaled_core_valid m s hs t ht) vp).mp reach
    exact scaled_avoids_boundary m hm n t ht i (by rw [hi];exact bp)

structure ThreePart (n : Nat) : Prop where
  words : ∀ q : Digits, q.valid (2^n*3) → q.nonconstant →
    (periodic (step (2^n*3)) q ↔ reaches (step (2^n*3)) (fromPair (2^n*3) (2^n,0)) q ∨
      n%2=0 ∧ reaches (step (2^n*3)) (fromPair (2^n*3) (2*2^n,0)) q)
  periods : if n%2=0 then digitPeriod (2^n*3) (2^n,0) (n+1) ∧ digitPeriod (2^n*3) (2*2^n,0) (2*(n+1))
    else digitPeriod (2^n*3) (2^n,0) (3*(n+1))

/-- Exactly two families indexed by `Fin (n/2)`. In each family the indices
below `(n-1)/2` give the full periods and the last index for even n gives
the half period. The disjointness statement is on actual words. -/
structure PowerPart (n : Nat) : Prop where
  words : ∀ q : Digits, q.valid (2^n) → q.nonconstant →
    (periodic (step (2^n)) q ↔ ∃ (i : Fin (n/2)) (anti : Bool),
      reaches (step (2^n)) (fromPair (2^n) (powerRep n i.val anti)) q)
  disjoint : ∀ (i j : Fin (n/2)) (b c : Bool),
    (reaches (step (2^n)) (fromPair (2^n) (powerRep n i.val b)) (fromPair (2^n) (powerRep n j.val c)) ↔
      i=j ∧ b=c)
  periods : ∀ i : Fin (n/2),
    digitPeriod (2^n) (powZ i.val) (if i.val<(n-1)/2 then n+2 else n/2+1) ∧
    digitPeriod (2^n) (powA n i.val) (if i.val<(n-1)/2 then n else n/2)
  odd_counts : n%2=1 → n/2=(n-1)/2
  even_counts : n%2=0 → (n-1)/2=(n-2)/2 ∧ n/2=(n-2)/2+1

private theorem power_part (n : Nat) (hn : 2≤n) : PowerPart n := by
  let B : Int := 2^n
  have valid (i : Fin (n/2)) (b : Bool) : validPair B (powerRep n i.val b) := by
    have h := pow_valid i.val (n-i.val-2)
    have ht := i.isLt
    rw [show i.val+(n-i.val-2)+2=n by omega] at h
    cases b <;> first | exact h.1 | exact h.2
  constructor
  · intro q hV hne
    constructor
    · intro per
      obtain ⟨p,vp,rfl,t,r,eq,ht,hr⟩ := (power_digit_classification n hn q hV ((nonconstant_iff q).mp hne)).mp per
      let i : Fin (n/2) := ⟨t,by omega⟩
      rcases hr with h|h
      · exact ⟨i,false,(word_reach_iff B _ p (valid i false) vp).mpr h⟩
      · exact ⟨i,true,(word_reach_iff B _ p (valid i true) vp).mpr h⟩
    · rintro ⟨i,b,h⟩
      have hi := i.isLt
      have eq : n=i.val+(n-i.val-2)+2 := by omega
      have v := valid i b
      have per : periodic (pairMap B) (powerRep n i.val b) := by
        apply (power_pair_classification n hn _ v).mpr
        refine ⟨i.val,n-i.val-2,eq,by omega,?_⟩
        cases b
        · exact Or.inl ⟨0,rfl⟩
        · exact Or.inr ⟨0,rfl⟩
      obtain ⟨j,hj⟩ := per
      exact (periodic_reachable _ _ _ ⟨j,(digit_return_iff B _ v _).mpr hj⟩ h).1
  · intro i j b c
    have hi := i.isLt
    have hj := j.isLt
    rw [word_reach_iff _ _ _ (valid i b) (valid j c),
      power_reps_disjoint n i.val (n-i.val-2) j.val (n-j.val-2) (by omega) (by omega) (by omega) (by omega)]
    exact ⟨fun h => ⟨Fin.ext h.1,h.2⟩,fun h => ⟨congrArg Fin.val h.1,h.2⟩⟩
  · intro i;exact (power_period_counts n i.val hn i.isLt).2.2
  · omega
  · omega

/-- P1 Theorem C and its unique-cycle corollary, for the positional four-slot
subtraction map. The outer factorization ranges over every even base. -/
structure TheoremCStatement (m : Int) (n : Nat) : Prop where
  sorting : ∀ q : Digits, (sort4 q).a ≥ (sort4 q).b ∧ (sort4 q).b ≥ (sort4 q).c ∧
    (sort4 q).c ≥ (sort4 q).d ∧ ∀ x, (sort4 q).count x=q.count x
  subtraction : ∀ q : Digits, q.valid (2^n*m) → (step (2^n*m) q).valid (2^n*m) ∧
    (step (2^n*m) q).value (2^n*m) = (sort4 q).value (2^n*m) - (sort4 q).rev.value (2^n*m)
  odd_core : 5≤m → LargePart m n
  core_three : m=3 → ThreePart n
  binary : m=1 → n=1 → ∀ q : Digits, q.valid 2 → q.nonconstant →
      (periodic (step 2) q ↔ q=⟨0,1,1,1⟩ ∨ q=⟨1,0,0,1⟩)
  binary_fixed : m=1 → n=1 →
      step 2 ⟨0,1,1,1⟩=⟨0,1,1,1⟩ ∧ step 2 ⟨1,0,0,1⟩=⟨1,0,0,1⟩
  pure_power : m=1 → 2≤n → PowerPart n
  uniqueness : uniqueDigitCycle (2^n*m) ↔ (m=3 ∧ n%2=1) ∨ (m=5 ∧ n%2=1) ∨ (m=7 ∧ n%3≠0)
  unique_lengths :
    (m=3 → n%2=1 → digitPeriod (2^n*m) (2^n,0) (6*(n/2+1))) ∧
    (m=5 → n%2=1 → digitPeriod (2^n*m) (scale (2^n) (3,1)) 1) ∧
    (m=7 → n%3≠0 → digitPeriod (2^n*m) (scale (2^n) (3,1)) 3)

theorem theoremC (m : Int) (n : Nat) (hm : 0<m) (hodd : m%2=1) (hn : 0<n) : TheoremCStatement m n := by
  have hs : (2:Int)≤2^n := by
    have lo := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) hn
    simp only [Int.pow_zero] at lo
    omega
  have hB : 2≤2^n*m := by
    have h := Int.mul_le_mul_of_nonneg_left (show 1≤m by omega) (show 0≤(2:Int)^n by omega)
    simp only [Int.mul_one] at h
    omega
  refine ⟨fun q => ⟨(sort4_ordered q).1,(sort4_ordered q).2.1,(sort4_ordered q).2.2,sort4_count q⟩,
    fun q hq => ⟨step_valid _ hB q hq,step_value _ q⟩,fun h5 => large_part m hodd h5 n hn, fun _ => ⟨fun q hq ne => three_digit_classification n hn q hq ((nonconstant_iff q).mp ne),three_digit_periods n hn⟩,
    fun _ _ q hq ne => base_two_cycles q hq ((nonconstant_iff q).mp ne), fun _ _ => by decide, fun _ h => power_part n h,
    unique_cycle m hodd hm n hn,?_⟩
  have per := unique_cycle_periods n hn
  exact ⟨fun eq => eq ▸ per.1,fun eq => eq ▸ per.2.1,fun eq => eq ▸ per.2.2⟩

private theorem odd_part_factorization (x : Nat) (hx : 0<x) :
    ∃ (m : Int) (n : Nat), 0<m ∧ m%2=1 ∧ (x:Int)=2^n*m := by
  induction x using Nat.strongRecOn with
  | ind x ih =>
    by_cases ho : x%2=1
    · exact ⟨x,0,by omega,by omega,by simp⟩
    · obtain ⟨m,n,hm,hodd,he⟩ := ih (x/2) (by omega) (by omega)
      refine ⟨m,n+1,hm,hodd,?_⟩
      have twice : (x:Int)=2*((x/2:Nat):Int) := by omega
      rw [twice,he,Int.pow_succ]
      grind

/-- Equivalent entry point quantified over the base itself. -/
theorem every_even_base (B : Int) (hB : 2≤B) (heven : B%2=0) :
    ∃ (m : Int) (n : Nat), 0<m ∧ m%2=1 ∧ 0<n ∧ B=2^n*m ∧ TheoremCStatement m n := by
  obtain ⟨m,n,hm,hodd,he⟩ := odd_part_factorization B.toNat (by omega)
  rw [Int.toNat_of_nonneg (show 0≤B by omega)] at he
  have hn : 0<n := by
    cases n with
    | zero => simp only [Int.pow_zero,Int.one_mul] at he;omega
    | succ k => omega
  exact ⟨m,n,hm,hodd,hn,he,theoremC m n hm hodd hn⟩

#print axioms every_even_base
#print axioms theoremC
end KaprekarBoundary
