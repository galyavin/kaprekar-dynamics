import PowerRows

namespace KaprekarBoundary

private theorem generic_gap (B : Int) (p : Pair) (h : p.2 ≤ p.1) :
    (genericMap B p).1-(genericMap B p).2 ≤ 2*(p.1-p.2) := by
  grind [genericMap,absInt,sortPair]

private theorem generic_order (B : Int) (p : Pair) : (genericMap B p).2≤(genericMap B p).1 := by
  grind [genericMap,sortPair]

private theorem generic_iter_order (B : Int) (p : Pair) (h : p.2≤p.1) (k : Nat) :
    (iter (genericMap B) k p).2≤(iter (genericMap B) k p).1 := by
  cases k with
  | zero => exact h
  | succ k => exact generic_order B _

private theorem generic_iter_gap (B : Int) (p : Pair) (h : p.2≤p.1) (k : Nat) :
    (iter (genericMap B) k p).1-(iter (genericMap B) k p).2 ≤ 2^k*(p.1-p.2) := by
  induction k with
  | zero => simp [iter]
  | succ k ih =>
    have bound := generic_gap B _ (generic_iter_order B p h k)
    simp only [iter,Int.pow_succ]
    grind

private theorem boundary_first_fold_gap (B : Int) (p : Pair) (hb : boundary B p) :
    let q := genericMap B (pairMap B p)
    q.1-q.2 ≤ 4 := by
  by_cases hz : p.2=0
  · simp only [pairMap,if_pos hz,genericMap]
    grind [absInt,sortPair]
  · have hd : (pairMap B p).1-(pairMap B p).2 ≤ 2 := by
      simp only [pairMap,if_neg hz,if_pos (hb.resolve_left hz)]
      grind [absInt,sortPair]
    have ho : (pairMap B p).2 ≤ (pairMap B p).1 := by
      simp only [pairMap,if_neg hz,if_pos (hb.resolve_left hz)]
      grind [sortPair]
    exact Int.le_trans (generic_gap B _ ho) (by omega)

/-- During a first return the folded pair has separation at most twice the
power of two at its final generic depth. -/
theorem return_gap (B : Int) (p q : Pair) (j : Nat)
    (hr : FirstReturn (pairMap B) (boundary B) p q ((j+1)+1)) :
    q.1-q.2 ≤ 2*2^(j+1) := by
  let z := pairMap B p
  have step (i : Nat) : iter (pairMap B) i z = iter (pairMap B) (i+1) p := by
    change iter (pairMap B) i (iter (pairMap B) 1 p) = _
    rw [← iter_add,Nat.add_comm]
  have hfree : ∀ i, i<j+1 → ¬boundary B (iter (pairMap B) i z) :=
    fun i hi => by rw [step]; exact hr.2.2.2.2 (i+1) (by omega) (by omega)
  have eq : iter (genericMap B) j (genericMap B z) = q := by
    change iter (genericMap B) j (iter (genericMap B) 1 z) = q
    rw [← iter_add,Nat.add_comm,← actual_generic_prefix _ _ _ hfree,step,hr.2.2.2.1]
  have bound := generic_iter_gap B (genericMap B z) (generic_order B z) j
  rw [eq] at bound
  have hg := boundary_first_fold_gap B p hr.1
  change (genericMap B z).1-(genericMap B z).2 ≤ 4 at hg
  have pos : (0:Int)<2^j := Int.pow_pos (by decide)
  have hh := Int.mul_le_mul_of_nonneg_left hg (Int.le_of_lt pos)
  rw [Int.pow_succ]
  grind

private theorem generic_factor (j : Nat) (L : Int) (p : Pair) :
    ∃ t, iter (genericMap (2^j*L)) j p = scale (2^j) t := by
  induction j generalizing p L with
  | zero => exact ⟨p,by simp [iter,scale]⟩
  | succ j ih =>
    obtain ⟨t,ht⟩ := ih (2*L) p
    rw [show 2^j*(2*L)=2^(j+1)*L by rw [Int.pow_succ];grind] at ht
    change ∃ t, genericMap _ (iter _ j p) = scale (2^(j+1)) t
    rw [ht, show 2^(j+1)*L=2^j*(2*L) by rw [Int.pow_succ];grind,
      generic_scale _ _ (Int.pow_pos (by decide : (0:Int)<2))]
    obtain ⟨w,he⟩ := even_generic L t
    exact ⟨w,by rw [he]; simp [scale,Int.pow_succ,Int.mul_assoc]⟩

theorem return_scale (B : Int) (p q : Pair) (j : Nat) (L : Int) (hdiv : B=2^j*L)
    (hr : FirstReturn (pairMap B) (boundary B) p q (j+1)) : ∃ w, q=scale (2^j) w := by
  let z := pairMap B p
  have step (i : Nat) : iter (pairMap B) i z = iter (pairMap B) (i+1) p := by
    change iter (pairMap B) i (iter (pairMap B) 1 p) = _
    rw [← iter_add,Nat.add_comm]
  have hfree : ∀ i, i<j → ¬boundary B (iter (pairMap B) i z) :=
    fun i hi => by rw [step]; exact hr.2.2.2.2 (i+1) (by omega) (by omega)
  have endeq : iter (genericMap B) (j) z=q := by
    rw [← actual_generic_prefix _ _ _ hfree, step,hr.2.2.2.1]
  rw [← endeq,hdiv]
  exact generic_factor (j) _ z

/-- Every target of a positive-depth boundary return has separation exactly
one or two units of the current power of two. A zero coordinate or an
anti-diagonal then specifies the target. -/
theorem positive_return_shape (m : Int) (hm : m%2=1) (n j : Nat) (p q : Pair)
    (hp : validPair (2^n*m) p)
    (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q ((j+1)+1)) :
    (q.2=0 ∨ q.1+q.2=2^n*m) ∧ (q.1-q.2=2^(j+1) ∨ q.1-q.2=2*2^(j+1)) := by
  let B := 2^n*m
  let d : Int := 2^(j+1)
  have hj := clock_bound m hm n (j+1) p q hp hr
  have hdiv : B=d*(2^(n-(j+1))*m) := by
    dsimp [B,d]
    conv => lhs; rw [show n=(j+1)+(n-(j+1)) by omega, Int.pow_add]
    grind
  have qv : validPair B q := hr.2.2.2.1 ▸ iter_pair_valid B p hp _
  have hd : ¬q.1=q.2 := fun he => by have := diagonal_first_return_one _ _ _ _ hr he; omega
  have hbound := hr.2.1
  change boundary B q at hbound
  have sep := return_gap B p q j hr
  have fact := return_scale B p q (j+1) (2^(n-(j+1))*m) hdiv hr
  obtain ⟨w,hw⟩ := fact
  have unit : 0<d := Int.pow_pos (by decide)
  have e : q.1-q.2=d*(w.1-w.2) := by rw [hw]; simp [scale, d, Int.mul_sub]
  have hpos : 0<q.1-q.2 := by have := qv.2.2.2; omega
  change q.1-q.2≤2*d at sep
  have upper : w.1-w.2≤2 := by
    apply Classical.byContradiction
    intro h
    have := Int.mul_lt_mul_of_pos_left (show 2<w.1-w.2 by omega) unit
    grind
  have lower : 0<w.1-w.2 := by
    apply Classical.byContradiction
    intro h
    have := Int.mul_le_mul_of_nonneg_left (show w.1-w.2≤0 by omega) (Int.le_of_lt unit)
    grind
  change (q.2=0 ∨ q.1+q.2=B) ∧ (q.1-q.2=d ∨ q.1-q.2=2*d)
  constructor
  · grind [boundary]
  · have : w.1-w.2=1 ∨ w.1-w.2=2 := by omega
    grind

/-- If no generic step intervenes, the only possible targets are the two
exceptional singletons or an anti-diagonal pair with gap two. -/
theorem direct_return_shape (H : Int) (p q : Pair) (hp : validPair (2*H) p)
    (hr : FirstReturn (pairMap (2*H)) (boundary (2*H)) p q 1) :
    (q=(1,1) ∧ p=(H,H)) ∨ (q=(2*H-1,0) ∧ p=(1,0)) ∨
      (q.1+q.2=2*H ∧ q.1-q.2=2) := by
  obtain ⟨hb,hq,_,he,_⟩ := hr
  change pairMap (2*H) p=q at he
  subst q
  rcases p with ⟨u,v⟩
  grind (splits := 64) [validPair, boundary, pairMap, sortPair, absInt]

/-- A boundary visit after depth `n` hits one of the three scaled graph
vertices. The anti-diagonal gap must be one scaled unit since m is odd. -/
theorem full_depth_target (m : Int) (hm : m%2=1) (n : Nat) (p q : Pair)
    (hn : 0<n) (hp : validPair (2^n*m) p)
    (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q (n+1)) :
    let s : Int := 2^n
    q=(s,0) ∨ q=(2*s,0) ∨ (q.1+q.2=s*m ∧ q.1-q.2=s) := by
  cases n with
  | zero => omega
  | succ k =>
    obtain ⟨hb, hsep⟩ := positive_return_shape m hm (k+1) k p q hp hr
    obtain ⟨w,hw⟩ := return_scale _ _ _ (k+1) m rfl hr
    let s : Int := 2^(k+1)
    have hs : 0<s := Int.pow_pos (by decide)
    have qe : q.1=s*w.1 ∧ q.2=s*w.2 := by rw [hw]; exact ⟨rfl,rfl⟩
    change q.2=0 ∨ q.1+q.2=s*m at hb
    change q.1-q.2=s ∨ q.1-q.2=2*s at hsep
    change q=(s,0) ∨ q=(2*s,0) ∨ _
    rcases hb with hz|ha
    · rcases hsep with e|e <;> grind
    · right; right
      refine ⟨ha, hsep.resolve_right ?_⟩
      intro bad
      have sum : w.1+w.2=m := by
        apply (Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)).mp
        rw [Int.mul_add]; grind
      have gap : w.1-w.2=2 := by
        apply (Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)).mp
        rw [Int.mul_sub]; grind
      omega

private def ready (m B : Int) (q : Pair) : Prop :=
  (q.2=0 ∧ ((q.1%2=0 ∧ ¬m∣q.1) ∨ (q.1%2=1 ∧ ¬m∣q.1-1))) ∨
  (q.2≠0 ∧ (q.1=q.2 ∨ q.1+q.2=B) ∧ ¬m∣4*q.1)

private theorem return_ready (m : Int) (hm : m%2=1) (h3 : 3≤m) (n j : Nat) (hn : 0<n)
    (p q : Pair) (hp : validPair (2^n*m) p)
    (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q (j+1)) : ready m (2^n*m) q := by
  let B : Int := 2^n*m
  have hmB : m∣B := Int.dvd_mul_left _ _
  have heven : ∃ H : Int, B=2*H := by
    cases n with
    | zero => omega
    | succ k => exact ⟨2^k*m, by dsimp [B]; rw [Int.pow_succ]; grind⟩
  obtain ⟨H,hH⟩ := heven
  have one : ¬m∣1 := fun h => by have := Int.le_of_dvd (by decide : (0:Int)<1) h; omega
  have pw (t : Nat) : ¬m∣(2:Int)^t := by
    simpa using (show ¬m∣2^t*(1:Int) by rw [odd_cancel_pow m 1 t hm]; exact one)
  have hq : validPair B q := hr.2.2.2.1 ▸ iter_pair_valid B p hp _
  have anti (a : Nat) (hs : q.1+q.2=B) (hd : q.1-q.2=2^a) : ready m B q := by
    right
    have hv : q.2≠0 := by have := hq.2.1; omega
    refine ⟨hv,Or.inr hs,fun h => ?_⟩
    have hh : m∣2*(q.1-q.2) := by
      rw [show 2*(q.1-q.2)=4*q.1-2*B by omega]
      exact Int.dvd_sub h (Int.dvd_mul_of_dvd_right hmB)
    exact pw a (hd ▸ (odd_cancel_two _ _ hm).mp hh)
  cases j with
  | zero =>
    change FirstReturn (pairMap B) (boundary B) p q 1 at hr
    change validPair B p at hp
    rw [hH] at hr hp
    rcases direct_return_shape H p q hp hr with ⟨ha,_⟩|⟨ha,_⟩|⟨hs,hd⟩
    · subst q
      right
      exact ⟨by decide,Or.inl rfl,by simpa using pw 2⟩
    · subst q
      left
      refine ⟨rfl,Or.inr ⟨by dsimp; omega,?_⟩⟩
      intro h
      have he : m∣B-(2*H-1-1) := Int.dvd_sub hmB h
      rw [show B-(2*H-1-1)=2 by omega] at he
      exact pw 1 he
    · exact anti 1 (by omega) hd
  | succ j =>
    obtain ⟨hs,hd⟩ := positive_return_shape m hm n j p q hp hr
    have zero (a : Nat) (ha : 0<a) (hz : q.2=0) (hd : q.1-q.2=2^a) : ready m B q := by
      left
      refine ⟨hz,Or.inl ⟨?_,?_⟩⟩
      · rw [show q.1=2^a by omega]
        cases a with
        | zero => omega
        | succ a => simp [Int.pow_succ]
      · rw [show q.1=2^a by omega]; exact pw a
    have hshift : (2:Int)*2^(j+1)=2^(j+2) := by rw [Int.pow_succ (n:=j+1)]; grind
    rcases hs with hz|hs <;> rcases hd with hd|hd
    · exact zero _ (by omega) hz hd
    · exact zero _ (by omega) hz (hd.trans hshift)
    · exact anti _ hs hd
    · exact anti _ hs (hd.trans hshift)

private theorem ready_depth (m : Int) (hm : m%2=1) (h3 : 3≤m) (n j : Nat)
    (p q : Pair) (hp : validPair (2^n*m) p) (ha : ready m (2^n*m) p)
    (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q (j+1)) : j=n := by
  have hbound := clock_bound m hm n j p q hp hr
  by_cases he : j=n
  · exact he
  have hj : j<n := by omega
  have hit := ((first_return_clock _ _ hp hr.1 j).mp ⟨q,hr⟩).1
  have odd (x : Int) (hx : x%2=1) := odd_no_early m x n j hj hx
  have core (x : Int) (hx : ¬m∣x) : ¬2^n*m∣2^j*x :=
    fun h => hx (odd_core_divides m x n j hm h)
  obtain ⟨u,v⟩ := p
  rcases ha with ⟨hz,heven|hodd⟩|⟨hv,hb,h4⟩
  · change v=0 at hz; subst v
    rw [zero_wall] at hit
    have h1 := odd 1 (by decide)
    simp only [Int.mul_one] at h1
    obtain ⟨hu,hn⟩ := heven
    change u%2=0 at hu
    exact False.elim (hit.elim (odd _ (by omega)) (fun h => h.elim (core _ hn)
      (fun g => g.elim (odd _ (by omega)) h1)))
  · change v=0 at hz; subst v
    rw [zero_wall] at hit
    have h1 := odd 1 (by decide)
    simp only [Int.mul_one] at h1
    obtain ⟨hu,hn⟩ := hodd
    change u%2=1 at hu
    exact False.elim (hit.elim (core _ hn) (fun h => h.elim (odd _ hu)
      (fun g => g.elim (odd _ (by omega)) h1)))
  · rw [diagonal_wall _ _ _ _ hv hb,Int.mul_comm 2 (2^j)] at hit
    have h2 : ¬m∣2 := fun h => by have := Int.le_of_dvd (by decide : (0:Int)<2) h; omega
    exact False.elim (hit.elim (odd _ (by omega)) (fun h => h.elim (odd _ (by omega))
      (fun g => g.elim (core _ h4) (core _ h2))))

/-- Two successive boundary returns at odd core at least three necessarily
end in the three-vertex graph. -/
theorem two_returns_to_graph (m : Int) (hm : m%2=1) (h3 : 3≤m) (n j k : Nat) (hn : 0<n)
    (p q z : Pair) (hp : validPair (2^n*m) p)
    (hr : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) p q (j+1))
    (ht : FirstReturn (pairMap (2^n*m)) (boundary (2^n*m)) q z (k+1)) :
    let s : Int := 2^n
    k=n ∧ (z=(s,0) ∨ z=(2*s,0) ∨ (z.1+z.2=s*m ∧ z.1-z.2=s)) := by
  have hq := hr.2.2.2.1 ▸ iter_pair_valid _ _ hp _
  have ha := return_ready m hm h3 n j hn p q hp hr
  have eq := ready_depth m hm h3 n k q z hq ha ht
  subst k
  exact ⟨rfl, full_depth_target m hm n q z hn hq ht⟩

#print axioms two_returns_to_graph
#print axioms direct_return_shape
#print axioms full_depth_target
#print axioms positive_return_shape
end KaprekarBoundary
