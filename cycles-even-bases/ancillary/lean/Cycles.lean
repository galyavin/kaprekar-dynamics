import Exhaustion

namespace KaprekarBoundary

def periodic {α : Type} (f : α → α) (x : α) : Prop := ∃ k : Nat, iter f (k+1) x=x

def reaches {α : Type} (f : α → α) (x y : α) : Prop := ∃ j : Nat, iter f j x=y

theorem reaches_trans {α : Type} (f : α → α) (x y z : α)
    (h : reaches f x y) (g : reaches f y z) : reaches f x z := by
  obtain ⟨i,hi⟩ := h
  obtain ⟨j,hj⟩ := g
  exact ⟨i+j,by rw [iter_add,hi,hj]⟩

theorem iter_commute {α : Type} (f : α → α) (i j : Nat) (x : α) :
    iter f i (iter f j x)=iter f j (iter f i x) := by
  rw [← iter_add,Nat.add_comm,iter_add]

theorem periodic_reachable {α : Type} (f : α → α) (x y : α) (hx : periodic f x) (h : reaches f x y) :
    periodic f y ∧ reaches f y x := by
  obtain ⟨k,hk⟩ := hx
  obtain ⟨i,rfl⟩ := h
  have forward (j : Nat) : iter f (k+1) (iter f j x)=iter f j x := by
    rw [iter_commute,hk]
  refine ⟨⟨k,forward i⟩,?_⟩
  induction i with
  | zero => exact ⟨0,rfl⟩
  | succ i ih =>
    apply reaches_trans f _ (iter f i x) _ _ ih
    refine ⟨k,?_⟩
    change iter f k (iter f 1 (iter f i x))=_
    rw [← iter_add,Nat.add_comm,forward]

private theorem first_hit {P : Nat → Prop} (k : Nat) (hk : P k) :
    ∃ i, P i ∧ ∀ j, j<i → ¬P j := by
  induction k using Nat.strongRecOn with
  | ind k ih =>
    by_cases h : ∃ j, j<k ∧ P j
    · obtain ⟨j,hj,hp⟩ := h
      exact ih j hj hp
    · exact ⟨k,hk,fun j hj hh => h ⟨j,hj,hh⟩⟩

theorem periodic_first_return {α : Type} (f : α → α) (b : α → Prop) (p : α)
    (hp : periodic f p) (hb : b p) :
    ∃ q j, FirstReturn f b p q (j+1) ∧ periodic f q ∧ reaches f p q ∧ reaches f q p := by
  obtain ⟨k,hk⟩ := hp
  obtain ⟨j,hj,hmin⟩ := first_hit (P:=fun j => b (iter f (j+1) p)) k (by change b (iter f (k+1) p); rw [hk];exact hb)
  have reach : reaches f p (iter f (j+1) p) := ⟨j+1,rfl⟩
  have hh := periodic_reachable f p _ ⟨k,hk⟩ reach
  refine ⟨_,j,⟨hb,hj,by omega,rfl,?_⟩,hh.1,reach,hh.2⟩
  rintro (_ | i) hi hj
  · omega
  · exact hmin i (by omega)

theorem periodic_meets_graph (m : Int) (hm : m%2=1) (h3 : 3≤m) (n : Nat) (hn : 0<n)
    (p : Pair) (hp : validPair (2^n*m) p) (hper : periodic (pairMap (2^n*m)) p)
    (hb : ∃ i, boundary (2^n*m) (iter (pairMap (2^n*m)) i p)) :
    let s : Int := 2^n
    ∃ z, reaches (pairMap (s*m)) p z ∧ reaches (pairMap (s*m)) z p ∧ periodic (pairMap (s*m)) z ∧
      (z=(s,0) ∨ z=(2*s,0) ∨ (z.1+z.2=s*m ∧ z.1-z.2=s)) := by
  obtain ⟨i,hib⟩ := hb
  let B := 2^n*m
  let a := iter (pairMap B) i p
  have pa : reaches (pairMap B) p a := ⟨i,rfl⟩
  have ha := periodic_reachable _ _ _ hper pa
  obtain ⟨q,j,hr,qper,aq,_⟩ := periodic_first_return _ _ a ha.1 hib
  obtain ⟨z,k,ht,zper,qz,_⟩ := periodic_first_return _ _ q qper hr.2.1
  have hgraph := two_returns_to_graph m hm h3 n j k hn a q z (iter_pair_valid _ p hp i) hr ht
  have pz := reaches_trans _ p a z pa (reaches_trans _ a q z aq qz)
  exact ⟨z,pz,(periodic_reachable _ p z hper pz).2,zper,hgraph.2⟩

/-- Every periodic pair is either in the lifted core, or on a cycle containing
a periodic critical vertex. The latter point itself reaches the input. -/
theorem cycle_core_or_graph (m : Int) (hm : m%2=1) (h3 : 3≤m) (n : Nat) (hn : 0<n)
    (p : Pair) (hp : validPair (2^n*m) p) (hper : periodic (pairMap (2^n*m)) p) :
    let s : Int := 2^n
    (∃ t, p=scale s t ∧ oddCore m t ∧ periodic (pairMap m) t) ∨
    ∃ z, reaches (pairMap (s*m)) z p ∧ periodic (pairMap (s*m)) z ∧
      (z=(s,0) ∨ z=(2*s,0) ∨ (z.1+z.2=s*m ∧ z.1-z.2=s)) := by
  by_cases h : ∃ i, boundary (2^n*m) (iter (pairMap (2^n*m)) i p)
  · obtain ⟨z,_,zp,zper,graph⟩ := periodic_meets_graph m hm h3 n hn p hp hper h
    exact Or.inr ⟨z,zp,zper,graph⟩
  · obtain ⟨k,hk⟩ := hper
    obtain ⟨t,eq,ht,kt⟩ := (generic_cycle_iff m hm n p hp k).mp ⟨fun i hi => h ⟨i,hi⟩,hk⟩
    exact Or.inl ⟨t,eq,ht,k,kt⟩

private theorem pair_period_of_digit (B : Int) (p : Pair) (hp : validPair B p) (k : Nat)
    (h : digitPeriod B p k) : periodic (pairMap B) p := by
  cases k with
  | zero => have := h.1; omega
  | succ k => exact ⟨k,(digit_return_iff _ _ hp _).mp h.2.1⟩

theorem three_pair_classification (n : Nat) (hn : 0<n) (p : Pair) (hp : validPair (2^n*3) p) :
    let s : Int := 2^n
    let f := pairMap (s*3)
    periodic f p ↔ reaches f (s,0) p ∨ (n%2=0 ∧ reaches f (2*s,0) p) := by
  let s : Int := 2^n
  let f := pairMap (s*3)
  have hs : 0<s := Int.pow_pos (by decide)
  have hS : validPair (s*3) (s,0) := by simp [validPair]; omega
  have hE : validPair (s*3) (2*s,0) := by simp [validPair]; omega
  have graph := three_critical_graph n hn
  have periods := three_digit_periods n hn
  change periodic f p ↔ reaches f (s,0) p ∨ (n%2=0 ∧ reaches f (2*s,0) p)
  constructor
  · intro hper
    rcases cycle_core_or_graph 3 (by decide) (by decide) n hn p hp hper with ⟨t,_,ht,_⟩|⟨z,zp,_,hz⟩
    · have := ht.1
      grind [oddCore,validPair]
    · rcases hz with he|he|⟨su,di⟩
      · exact Or.inl (he ▸ zp)
      · by_cases hpar : n%2=0
        · exact Or.inr ⟨hpar,he ▸ zp⟩
        · simp only [if_neg hpar] at graph
          have su : reaches f (s,0) (2*s,s) := ⟨n+1,graph.2.1.2.2.2.1⟩
          have ue : reaches f (2*s,s) (2*s,0) := ⟨n+1,graph.1.2.2.2.1⟩
          exact Or.inl (reaches_trans f _ _ _ (reaches_trans f _ _ _ su ue) (he ▸ zp))
      · have he : z=(2*s,s) := by apply Prod.ext <;> dsimp <;> change z.1+z.2=s*3 at su <;> change z.1-z.2=s at di <;> omega
        by_cases hpar : n%2=0
        · simp only [if_pos hpar] at graph
          have eu : reaches f (2*s,0) (2*s,s) := ⟨n+1,graph.2.2.2.2.2.1⟩
          exact Or.inr ⟨hpar,reaches_trans f _ _ _ eu (he ▸ zp)⟩
        · simp only [if_neg hpar] at graph
          have su : reaches f (s,0) (2*s,s) := ⟨n+1,graph.2.1.2.2.2.1⟩
          exact Or.inl (reaches_trans f _ _ _ su (he ▸ zp))
  · rintro (h|⟨hpar,h⟩)
    · by_cases hpar : n%2=0
      · simp only [if_pos hpar] at periods
        exact (periodic_reachable f _ _ (pair_period_of_digit _ _ hS _ periods.1) h).1
      · simp only [if_neg hpar] at periods
        exact (periodic_reachable f _ _ (pair_period_of_digit _ _ hS _ periods) h).1
    · simp only [if_pos hpar] at periods
      exact (periodic_reachable f _ _ (pair_period_of_digit _ _ hE _ periods.2) h).1

theorem word_periodic_iff (B : Int) (q : Digits) (hq : q.valid B) (hne : (diff q).1≠0) :
    periodic (step B) q ↔ ∃ p, validPair B p ∧ q=fromPair B p ∧ periodic (pairMap B) p := by
  constructor
  · rintro ⟨k,hk⟩
    obtain ⟨p,hp,hq,_⟩ := periodic_digit_has_pair B q hq hne k hk
    refine ⟨p,hp,hq.symm,k,(digit_return_iff _ _ hp _).mp ?_⟩
    rwa [hq]
  · rintro ⟨p,hp,rfl,k,hk⟩
    exact ⟨k,(digit_return_iff _ _ hp _).mpr hk⟩

theorem three_digit_classification (n : Nat) (hn : 0<n) (q : Digits)
    (hq : q.valid (2^n*3)) (hne : (diff q).1≠0) :
    let s : Int := 2^n
    let B := s*3
    periodic (step B) q ↔ (reaches (step B) (fromPair B (s,0)) q ∨
      n%2=0 ∧ reaches (step B) (fromPair B (2*s,0)) q) := by
  let s : Int := 2^n
  let B := s*3
  have hs : 0<s := Int.pow_pos (by decide)
  have hS : validPair B (s,0) := by dsimp [validPair,B]; omega
  have hE : validPair B (2*s,0) := by dsimp [validPair,B]; omega
  have lift (a : Pair) (ha : validPair B a) (p : Pair) (h : reaches (pairMap B) a p) :
      reaches (step B) (fromPair B a) (fromPair B p) := by
    obtain ⟨i,hi⟩ := h
    exact ⟨i,by rw [iter_fromPair B a ha,hi]⟩
  constructor
  · intro h
    obtain ⟨p,hp,rfl,hper⟩ := (word_periodic_iff B q hq hne).mp h
    rcases (three_pair_classification n hn p hp).mp hper with h|⟨hpar,h⟩
    · exact Or.inl (lift _ hS p h)
    · exact Or.inr ⟨hpar,lift _ hE p h⟩
  · intro h
    rcases h with ⟨i,hi⟩|⟨hpar,i,hi⟩
    · have h := (three_pair_classification n hn (s,0) hS).mpr (Or.inl ⟨0,rfl⟩)
      obtain ⟨k,hk⟩ := h
      exact (periodic_reachable _ _ _ ⟨k,(digit_return_iff _ _ hS _).mpr hk⟩ ⟨i,hi⟩).1
    · have h := (three_pair_classification n hn (2*s,0) hE).mpr (Or.inr ⟨hpar,0,rfl⟩)
      obtain ⟨k,hk⟩ := h
      exact (periodic_reachable _ _ _ ⟨k,(digit_return_iff _ _ hE _).mpr hk⟩ ⟨i,hi⟩).1

#print axioms three_digit_classification
#print axioms periodic_meets_graph
end KaprekarBoundary
