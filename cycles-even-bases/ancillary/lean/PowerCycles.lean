import OddGraph

namespace KaprekarBoundary

def powZ (t : Nat) : Pair := (2^(t+1),0)
def powA (n t : Nat) : Pair := (2^(n-1)+2^t,2^(n-1)-2^t)

private theorem pow_index (a n : Nat) (h : (2:Int)^a<2^n) : a<n := by
  by_cases he : a=n
  · subst n; omega
  · by_cases hl : a<n
    · exact hl
    · have := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) (show n<a by omega)
      omega

private theorem half_pow (n : Nat) (hn : 0<n) : (2:Int)^n=2*2^(n-1) := by
  cases n with
  | zero => omega
  | succ k => simp [Int.pow_succ,Int.mul_comm]

private theorem power_first_shape (n : Nat) (hn : 2≤n) (p q : Pair) (hp : validPair (2^n) p) (j : Nat)
    (hr : FirstReturn (pairMap (2^n)) (boundary (2^n)) p q (j+1)) :
    (q=(1,1) ∧ p=(2^(n-1),2^(n-1))) ∨ (q=(2^n-1,0) ∧ p=(1,0)) ∨
    ∃ t r, n=t+r+2 ∧ (q=powZ t ∨ q=powA n t) := by
  have hv : validPair (2^n) q := hr.2.2.2.1 ▸ iter_pair_valid _ _ hp _
  have fromGap (t : Nat) (hgap : q.1-q.2=2^(t+1)) (hshape : q.2=0 ∨ q.1+q.2=2^n) :
      ∃ a b, n=a+b+2 ∧ (q=powZ a ∨ q=powA n a) := by
    have ht := pow_index (t+1) n (by have := hv.2.1; have := hv.2.2.1;omega)
    refine ⟨t,n-t-2,by omega,?_⟩
    rcases hshape with h|h
    · left; apply Prod.ext <;> dsimp [powZ] <;> omega
    · right
      have half := half_pow n (by omega)
      rw [Int.pow_succ] at hgap
      apply Prod.ext <;> dsimp [powA] <;> omega
  cases j with
  | zero =>
    have half := half_pow n (by omega)
    rw [half] at hp hr
    rcases direct_return_shape _ _ _ hp hr with h|h|h
    · exact Or.inl h
    · exact Or.inr (Or.inl (by simpa only [← half] using h))
    · exact Or.inr (Or.inr (fromGap 0 h.2 (Or.inr (by simpa only [← half] using h.1))))
  | succ j =>
    have h := positive_return_shape 1 (by decide) n j p q (by simpa using hp) (by simpa using hr)
    simp only [Int.mul_one] at h
    right;right
    rcases h.2 with hd|hd
    · exact fromGap j hd h.1
    · have e : (2:Int)*2^(j+1)=2^((j+1)+1) := by rw [Int.pow_succ (n:=j+1)];grind
      exact fromGap (j+1) (hd.trans e) h.1

private theorem power_meets (n : Nat) (hn : 2≤n) (p : Pair)
    (hp : validPair (2^n) p) (hper : periodic (pairMap (2^n)) p) :
    ∃ t r, n=t+r+2 ∧ (reaches (pairMap (2^n)) (powZ t) p ∨
      reaches (pairMap (2^n)) (powA n t) p) := by
  let B : Int := 2^n
  have hb : ∃ i, boundary B (iter (pairMap B) i p) := by
    apply Classical.byContradiction
    intro hehe
    obtain ⟨k,hk⟩ := hper
    obtain ⟨t,_,ht,_⟩ := (generic_cycle_iff 1 (by decide) n p (by simpa using hp) k).mp
      ⟨by simpa using (show ∀ i, ¬boundary B (iter (pairMap B) i p) from fun i h => hehe ⟨i,h⟩),
       by simpa using hk⟩
    grind [oddCore,validPair]
  obtain ⟨i,hi⟩ := hb
  let a := iter (pairMap B) i p
  have pa : reaches (pairMap B) p a := ⟨i,rfl⟩
  have aPer := (periodic_reachable _ p a hper pa).1
  obtain ⟨q,j,hr,qp,aq,_⟩ := periodic_first_return _ _ a aPer hi
  have pq := reaches_trans _ _ _ _ pa aq
  have finish (z : Pair) (pz : reaches (pairMap B) p z) (hz : ∃ t r, n=t+r+2 ∧ (z=powZ t ∨ z=powA n t)) :
      ∃ t r, n=t+r+2 ∧ (reaches (pairMap B) (powZ t) p ∨ reaches (pairMap B) (powA n t) p) := by
    obtain ⟨t,r,hn,he|he⟩ := hz
    · exact ⟨t,r,hn,Or.inl (he ▸ (periodic_reachable _ p z hper pz).2)⟩
    · exact ⟨t,r,hn,Or.inr (he ▸ (periodic_reachable _ p z hper pz).2)⟩
  have hh : (1:Int)<2^(n-1) := by
    have := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) (show 0<n-1 by omega)
    simpa using this
  have hB := half_pow n (by omega)
  have next (bad : q=(1,1) ∨ q=(B-1,0)) :
      ∃ t r, n=t+r+2 ∧ (reaches (pairMap B) (powZ t) p ∨ reaches (pairMap B) (powA n t) p) := by
    obtain ⟨z,k,ht,_,qz,_⟩ := periodic_first_return _ _ q qp hr.2.1
    have hq : validPair B q := hr.2.2.2.1 ▸ iter_pair_valid B a (iter_pair_valid B p hp i) _
    rcases power_first_shape n hn q z hq k ht with h|h|h
    · have hhq := h.2
      rcases bad with eq|eq <;> grind
    · have hhq := h.2
      rcases bad with eq|eq <;> dsimp [B] at * <;> grind
    · exact finish z (reaches_trans _ _ _ _ pq qz) h
  rcases power_first_shape n hn a q (iter_pair_valid _ _ hp i) j hr with h|h|h
  · exact next (Or.inl h.1)
  · exact next (Or.inr h.1)
  · exact finish q pq h

private theorem pow_data (t r : Nat) :
    (2:Int)^(t+r+2)=2^(t+1)*2^(r+1) ∧ (2:Int)^((t+r+2)-1)=2*2^t*2^r := by
  rw [show t+r+2=t+1+(r+1) by omega,Int.pow_add]
  rw [show t+1+(r+1)-1=(t+r)+1 by omega,Int.pow_succ (n:=t+r),Int.pow_add]
  grind

theorem power_rep_returns (t r : Nat) :
    let n := t+r+2
    let f := pairMap (2^n)
    let b := boundary (2^n)
    FirstReturn f b (powZ t) (powZ r) (r+2) ∧ FirstReturn f b (powA n t) (powA n r) (r+1) := by
  have ha := power_anti_return t r
  have hz := power_zero_return t r
  have data := pow_data t r
  dsimp only
  constructor
  · rw [data.1]; exact hz
  · have bEq : (2:Int)^(t+r+2)=2*(2*2^t*2^r) := by
      rw [data.1,Int.pow_succ (n:=t),Int.pow_succ (n:=r)]
      grind
    simp only [bEq,powA,data.2]
    exact ha

theorem power_rep_periods (t r : Nat) :
    let n := t+r+2
    digitPeriod (2^n) (powZ t) (if t=r then r+2 else n+2) ∧
    digitPeriod (2^n) (powA n t) (if t=r then r+1 else n) := by
  have ha := power_anti_period t r
  have hz := power_zero_period t r
  have data := pow_data t r
  dsimp only
  constructor
  · rw [data.1]
    simpa only [powZ,show t+1+(r+1)+2=t+r+2+2 by omega] using hz
  · have bEq : (2:Int)^(t+r+2)=2*(2*2^t*2^r) := by
      rw [data.1,Int.pow_succ (n:=t),Int.pow_succ (n:=r)]
      grind
    simp only [bEq,powA,data.2]
    exact ha

theorem pow_valid (t r : Nat) :
    let n := t+r+2
    validPair (2^n) (powZ t) ∧ validPair (2^n) (powA n t) := by
  have hx : (0:Int)<2^t := Int.pow_pos (by decide)
  have hy : (0:Int)<2^r := Int.pow_pos (by decide)
  have pos := Int.mul_lt_mul_of_pos_left (show 1<2*(2:Int)^r by omega) hx
  have h := pow_data t r
  have b : (2:Int)^(t+r+2)=2*(2*2^t*2^r) := by rw [h.1,Int.pow_succ (n:=t),Int.pow_succ (n:=r)];grind
  dsimp only
  simp only [validPair,powZ,powA,h.2,b,Int.pow_succ (n:=t)]
  grind

theorem power_pair_classification (n : Nat) (hn : 2≤n) (p : Pair) (hp : validPair (2^n) p) :
    periodic (pairMap (2^n)) p ↔ ∃ t r, n=t+r+2 ∧ t≤r ∧
      (reaches (pairMap (2^n)) (powZ t) p ∨ reaches (pairMap (2^n)) (powA n t) p) := by
  constructor
  · intro h
    obtain ⟨t,r,rfl,hr⟩ := power_meets n hn p hp h
    by_cases ht : t≤r
    · exact ⟨t,r,rfl,ht,hr⟩
    · refine ⟨r,t,by omega,by omega,?_⟩
      have rev := power_rep_returns r t
      rw [show r+t+2=t+r+2 by omega] at rev
      rcases hr with h|h
      · exact Or.inl (reaches_trans _ _ _ _ ⟨_,rev.1.2.2.2.1⟩ h)
      · exact Or.inr (reaches_trans _ _ _ _ ⟨_,rev.2.2.2.2.1⟩ h)
  · rintro ⟨t,r,rfl,_,h⟩
    have periods := power_rep_periods t r
    have hV := pow_valid t r
    have fromPeriod (z : Pair) (v : validPair (2^(t+r+2)) z) (k : Nat)
        (hp : digitPeriod (2^(t+r+2)) z k) : periodic (pairMap (2^(t+r+2))) z := by
      cases k with
      | zero => have := hp.1;omega
      | succ k => exact ⟨k,(digit_return_iff _ _ v _).mp hp.2.1⟩
    rcases h with h|h
    · exact (periodic_reachable _ _ _ (fromPeriod _ hV.1 _ periods.1) h).1
    · exact (periodic_reachable _ _ _ (fromPeriod _ hV.2 _ periods.2) h).1

theorem power_digit_classification (n : Nat) (hn : 2≤n) (q : Digits) (hq : q.valid (2^n))
    (hne : (diff q).1≠0) :
    periodic (step (2^n)) q ↔ ∃ p, validPair (2^n) p ∧ q=fromPair (2^n) p ∧
      ∃ t r, n=t+r+2 ∧ t≤r ∧
      (reaches (pairMap (2^n)) (powZ t) p ∨ reaches (pairMap (2^n)) (powA n t) p) := by
  rw [word_periodic_iff _ _ hq hne]
  constructor
  · rintro ⟨p,hp,eq,h⟩
    exact ⟨p,hp,eq,(power_pair_classification n hn p hp).mp h⟩
  · rintro ⟨p,hp,eq,h⟩
    exact ⟨p,hp,eq,(power_pair_classification n hn p hp).mpr h⟩

private theorem two_return_members {α : Type} (f : α → α) (b : α → Prop) (x y z : α) (a c : Nat)
    (h : FirstReturn f b x y a) (g : FirstReturn f b y x c) (hz : b z) (hr : reaches f x z) :
    z=x ∨ z=y := by
  obtain ⟨i,hi⟩ := hr
  have period : iter f (a+c) x=x := by rw [iter_add,h.2.2.2.1,g.2.2.2.1]
  induction i using Nat.strongRecOn with
  | ind i ih =>
    by_cases iz : i=0
    · subst i; exact Or.inl hi.symm
    by_cases ea : i=a
    · subst i; exact Or.inr (hi.symm.trans h.2.2.2.1)
    by_cases low : i<a
    · exact False.elim (h.2.2.2.2 i (by omega) low (hi ▸ hz))
    by_cases mid : i<a+c
    · have he : iter f (i-a) y=z := by rw [← h.2.2.2.1,← iter_add,show a+(i-a)=i by omega,hi]
      exact False.elim (g.2.2.2.2 (i-a) (by omega) (by omega) (he ▸ hz))
    · have ha := h.2.2.1
      apply ih (i-(a+c)) (by omega)
      rw [← period,← iter_add,show a+c+(i-(a+c))=i by omega,hi]

def powerRep (n t : Nat) (anti : Bool) : Pair := if anti then powA n t else powZ t

private theorem rep_eq (n t r j k : Nat) (hn : n=t+r+2) (hn' : n=j+k+2)
    (b c : Bool) (h : powerRep n t b = powerRep n j c) : t=j ∧ b=c := by
  have v := pow_valid t r
  have w := pow_valid j k
  have half := half_pow n (by omega)
  have va : (powA n t).1+(powA n t).2=2^n := by dsimp [powA];omega
  have wa : (powA n j).1+(powA n j).2=2^n := by dsimp [powA];omega
  rw [← hn] at v
  rw [← hn'] at w
  have ppos : 0<(powA n t).2 := by have := v.2.2.1;omega
  have qpos : 0<(powA n j).2 := by have := w.2.2.1;omega
  have unique (x y : Nat) (hh : (2:Int)^x=2^y) : x=y := by
    by_cases eq : x=y
    · exact eq
    · rcases (show x<y ∨ y<x by omega) with lt|lt <;>
        have := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) lt <;> omega
  cases b <;> cases c <;> simp only [powerRep,Bool.false_eq_true,if_false,if_true] at h
  · have eq : (2:Int)^(t+1)=2^(j+1) := congrArg Prod.fst h
    exact ⟨by have := unique _ _ eq; omega,rfl⟩
  · have eq := congrArg Prod.snd h
    change 0=(powA n j).2 at eq
    omega
  · have eq := congrArg Prod.snd h
    change (powA n t).2=0 at eq
    omega
  · have eq := congrArg Prod.fst h
    change 2^(n-1)+(2:Int)^t=2^(n-1)+2^j at eq
    exact ⟨unique _ _ (by omega),rfl⟩

/-- Distinct canonical indices and distinct families give disjoint cycles.
Together with exhaustion this identifies the pair of families, each indexed
by precisely the integers `t < n/2`. -/
theorem power_reps_disjoint (n t r j k : Nat) (ht : n=t+r+2) (hj : n=j+k+2)
    (htr : t≤r) (hjk : j≤k) (b c : Bool) :
    reaches (pairMap (2^n)) (powerRep n t b) (powerRep n j c) ↔ t=j ∧ b=c := by
  constructor
  · intro h
    have f := power_rep_returns t r
    have g := power_rep_returns r t
    have dest := power_rep_returns j k
    rw [← ht] at f
    rw [show r+t+2=n by omega] at g
    rw [← hj] at dest
    have bound : boundary (2^n) (powerRep n j c) := by cases c <;> first | exact dest.1.1 | exact dest.2.1
    have eq : powerRep n j c=powerRep n t b ∨ powerRep n j c=powerRep n r b := by
      cases b
      · exact two_return_members _ _ _ _ _ _ _ f.1 g.1 bound h
      · exact two_return_members _ _ _ _ _ _ _ f.2 g.2 bound h
    rcases eq with eq|eq
    · exact rep_eq n t r j k ht hj b c eq.symm
    · obtain ⟨hjr,hcb⟩ := rep_eq n j k r t hj (by omega) c b eq
      exact ⟨by omega,hcb.symm⟩
  · rintro ⟨rfl,rfl⟩;exact ⟨0,rfl⟩

theorem power_index_count (n t : Nat) :
    (∃ r, n=t+r+2 ∧ t≤r) ↔ t<n/2 := by
  constructor
  · rintro ⟨r,hr,ht⟩;omega
  · intro h; exact ⟨n-t-2,by omega,by omega⟩

/-- The doubled cycles occupy precisely the first `(n-1)/2` indices.
Even exponents contribute the remaining fixed index in each family. -/
theorem power_period_counts (n t : Nat) (hn : 2≤n) (ht : t<n/2) :
    let r := n-t-2
    (t<r ↔ t<(n-1)/2) ∧ (t=r ↔ n%2=0 ∧ t+1=n/2) ∧
    digitPeriod (2^n) (powZ t) (if t<(n-1)/2 then n+2 else n/2+1) ∧
    digitPeriod (2^n) (powA n t) (if t<(n-1)/2 then n else n/2) := by
  have eq : n=t+(n-t-2)+2 := by omega
  have per := power_rep_periods t (n-t-2)
  rw [← eq] at per
  dsimp only
  refine ⟨by omega,by omega,?_⟩
  by_cases hs : t<(n-1)/2
  · simpa only [if_pos hs,if_neg (show t≠n-t-2 by omega)] using per
  · simpa only [if_neg hs,if_pos (show t=n-t-2 by omega),show n-t-2+2=n/2+1 by omega,
      show n-t-2+1=n/2 by omega] using per

#print axioms power_reps_disjoint
#print axioms power_period_counts
#print axioms power_digit_classification
#print axioms power_rep_periods
end KaprekarBoundary
