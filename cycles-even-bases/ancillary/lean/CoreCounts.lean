import PowerCycles

namespace KaprekarBoundary

private theorem stream_repeats {α : Type} (xs : List α) (a : Nat → α) (ha : ∀ i, a i ∈ xs) :
    ∃ i j, i<j ∧ a i=a j := by
  induction xs generalizing a with
  | nil => have := ha 0; simp at this
  | cons x xs ih =>
    by_cases h : ∃ i, a i=x
    · obtain ⟨i,hi⟩ := h
      by_cases hnext : ∃ j, i<j ∧ a j=x
      · obtain ⟨j,lt,hj⟩ := hnext
        exact ⟨i,j,lt,hi.trans hj.symm⟩
      · have tail (k : Nat) : a (i+1+k)∈xs := by
          have hh := ha (i+1+k)
          rw [List.mem_cons] at hh
          exact hh.resolve_left (fun eq => hnext ⟨_,by omega,eq⟩)
        obtain ⟨u,v,lt,eq⟩ := ih (fun k => a (i+1+k)) tail
        exact ⟨_,_,by omega,eq⟩
    · apply ih a
      intro i
      exact (List.mem_cons.mp (ha i)).resolve_left (fun eq => h ⟨i,eq⟩)

private theorem core_inj (m : Int) (hm : m%2=1) (p q : Pair) (hp : oddCore m p) (hq : oddCore m q)
    (eq : genericMap m p=genericMap m q) : p=q := by
  rcases p with ⟨u,v⟩;rcases q with ⟨a,b⟩
  grind (splits := 60) [oddCore,validPair,genericMap,absInt,sortPair]

theorem oddCore_periodic (m : Int) (hm : m%2=1) (p : Pair) (hp : oddCore m p) : periodic (pairMap m) p := by
  let f := genericMap m
  let ls : List Int := (List.range m.toNat).map Int.ofNat
  let pairs : List Pair := ls.flatMap (fun x => ls.map (fun y => (x,y)))
  have mem (x : Int) (lo : 0≤x) (hi : x<m) : x∈ls := by
    apply List.mem_map.mpr
    exact ⟨x.toNat,List.mem_range.mpr (by omega),Int.toNat_of_nonneg lo⟩
  have allmem (i : Nat) : iter f i p ∈ pairs := by
    have hc := (oddCore_orbit m hm p hp i).1
    let t := iter f i p
    change validPair m t at hc
    apply List.mem_flatMap.mpr
    refine ⟨t.1,mem _ (by have := hc.1; omega) hc.2.1,List.mem_map.mpr ?_⟩
    exact ⟨t.2,mem _ hc.2.2.1 (by have := hc.2.2.2; have := hc.2.1;omega), rfl⟩
  obtain ⟨i,j,hlt,eq⟩ := stream_repeats pairs (fun k => iter f k p) allmem
  have cancel (k : Nat) (x y : Pair) (hx : oddCore m x) (hy : oddCore m y)
      (eq : iter f k x=iter f k y) : x=y := by
    induction k with
    | zero => exact eq
    | succ k ih => exact ih (core_inj m hm _ _ (oddCore_orbit m hm x hx k) (oddCore_orbit m hm y hy k) eq)
  have he : p=iter f (j-i) p := by
    apply cancel i p _ hp (oddCore_orbit m hm p hp _)
    rw [← iter_add,show j-i+i=j by omega];exact eq
  refine ⟨j-i-1,?_⟩
  rw [show j-i-1+1=j-i by omega,generic_actual_iter m p (fun k => oddCore_away m hm _ (oddCore_orbit m hm p hp k))]
  exact he.symm

theorem scaled_core_valid (m s : Int) (hs : 0<s) (p : Pair) (hp : oddCore m p) :
    validPair (s*m) (scale s p) := by
  have h := hp.1
  have lo := Int.mul_pos hs (show 0<p.1 by have := h.1;omega)
  exact ⟨by change 1≤s*p.1;omega,Int.mul_lt_mul_of_pos_left h.2.1 hs,
    Int.mul_nonneg (Int.le_of_lt hs) h.2.2.1,Int.mul_le_mul_of_nonneg_left h.2.2.2 (Int.le_of_lt hs)⟩

theorem iter_scaled_core (m : Int) (hm : m%2=1) (n : Nat) (p : Pair) (hp : oddCore m p) (i : Nat) :
    iter (pairMap (2^n*m)) i (scale (2^n) p)=scale (2^n) (iter (pairMap m) i p) := by
  let s : Int := 2^n
  have hs : 0<s := Int.pow_pos (by decide)
  have hcore (k : Nat) : oddCore m (iter (pairMap m) k p) := by
    rw [generic_actual_iter m p (fun k => oddCore_away m hm _ (oddCore_orbit m hm p hp k))]
    exact oddCore_orbit m hm p hp k
  induction i with
  | zero => rfl
  | succ k ih =>
    let z := iter (pairMap m) k p
    have hz := hcore k
    have away := oddCore_away m hm z hz
    have sc : ¬boundary (s*m) (scale s z) := by
      simp only [boundary,scale,← Int.mul_add,Int.mul_eq_zero,Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)]
      grind [boundary]
    change pairMap (s*m) (iter (pairMap (s*m)) k (scale s p))=scale s (pairMap m z)
    rw [ih,map_generic _ _ sc,map_generic _ _ away,generic_scale s m hs z]

theorem scaled_avoids_boundary (m : Int) (hm : m%2=1) (n : Nat) (p : Pair) (hp : oddCore m p) (i : Nat) :
    ¬boundary (2^n*m) (iter (pairMap (2^n*m)) i (scale (2^n) p)) := by
  obtain ⟨k,hk⟩ := oddCore_periodic m hm p hp
  exact ((generic_cycle_iff m hm n _ (scaled_core_valid m (2^n) (Int.pow_pos (by decide)) p hp) k).mpr
    ⟨p,rfl,hp,hk⟩).1 i

theorem scaled_periodic (m : Int) (hm : m%2=1) (n : Nat) (p : Pair) (hp : oddCore m p) :
    periodic (pairMap (2^n*m)) (scale (2^n) p) := by
  obtain ⟨k,hk⟩ := oddCore_periodic m hm p hp
  exact ⟨k,by rw [iter_scaled_core m hm n p hp,hk]⟩

private def adjacent (m : Int) (p : Pair) : Prop :=
  p.1=absInt (2*p.2-m) ∨ p.2=absInt (2*p.1-m)

private theorem adjacent_next (m : Int) (p : Pair) (ha : adjacent m p) : adjacent m (genericMap m p) := by
  grind [genericMap,sortPair,adjacent]

theorem two_large_core_cycles (m : Int) (hm : m%2=1) (h9 : 9≤m) (n : Nat) :
    let x := scale (2^n) (m-2,1)
    let y := scale (2^n) (3,1)
    periodic (pairMap (2^n*m)) x ∧ periodic (pairMap (2^n*m)) y ∧ ¬reaches (pairMap (2^n*m)) x y := by
  have hx : oddCore m (m-2,1) := by dsimp [oddCore,validPair];omega
  have hy : oddCore m (3,1) := by dsimp [oddCore,validPair];omega
  refine ⟨scaled_periodic m hm n _ hx,scaled_periodic m hm n _ hy,?_⟩
  rintro ⟨i,hi⟩
  rw [iter_scaled_core m hm n _ hx] at hi
  have hs : (0:Int)<2^n := Int.pow_pos (by decide)
  have eq : iter (pairMap m) i (m-2,1)=(3,1) := by
    apply Prod.ext <;> apply (Int.mul_eq_mul_left_iff (Int.ne_of_gt hs)).mp
    · exact congrArg Prod.fst hi
    · exact congrArg Prod.snd hi
  have free : away m (m-2,1) := fun k => oddCore_away m hm _ (oddCore_orbit m hm _ hx k)
  rw [generic_actual_iter _ _ free] at eq
  have hadj (j : Nat) : adjacent m (iter (genericMap m) j (m-2,1)) := by
    induction j with
    | zero => simp only [iter,adjacent];grind [absInt]
    | succ j ih => exact adjacent_next m _ ih
  have bad := hadj i
  rw [eq] at bad
  grind [adjacent,absInt]

#print axioms two_large_core_cycles
end KaprekarBoundary
