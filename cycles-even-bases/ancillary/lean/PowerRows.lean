import SmallCores

namespace KaprekarBoundary

private theorem higher_wall (B x y : Int) (a b : Nat) (h : a≤b) :
    wall B (2^a*x) (2^a*y) → wall B (2^b*x) (2^b*y) := by
  unfold wall
  rw [← Int.mul_add, ← Int.mul_sub, ← Int.mul_add, ← Int.mul_sub]
  have lift (z : Int) (hz : B∣2^a*z) : B∣2^b*z := by
    have id : (2:Int)^b*z = 2^(b-a)*(2^a*z) := by
      conv => lhs; rw [show b=(b-a)+a by omega, Int.pow_add]
      grind
    rw [id]
    exact Int.dvd_mul_of_dvd_right hz
  grind

/-- For a positive generic depth, only the immediately preceding wall test is
needed to establish minimality. -/
theorem first_return_next (B : Int) (p q : Pair) (hp : validPair B p) (hb : boundary B p) (k : Nat) :
    FirstReturn (pairMap B) (boundary B) p q ((k+1)+1) ↔
    q=returnTarget B p (k+1) ∧
    wall B (2^(k+1)*(pairMap B p).1) (2^(k+1)*(pairMap B p).2) ∧
    ¬wall B (2^k*(pairMap B p).1) (2^k*(pairMap B p).2) := by
  rw [first_return_iff B p q hp hb (k+1)]
  exact ⟨fun h => ⟨h.1,h.2.1,h.2.2 k (by omega)⟩, fun h =>
    ⟨h.1,h.2.1, fun i hi hw => h.2.2 (higher_wall _ _ _ i k (by omega) hw)⟩⟩

/-- Row 3 of the P1 table. All legal zero-boundary parameters in either
specified family first return to the displayed power of two. -/
theorem zero_block_return (k : Nat) (L a e : Int) (hL : 2≤L) (he : e=0 ∨ e=1)
    (u : Int) (hu : u=L*(2*a+1)+e)
    (hp : validPair (2^(k+1)*L) (u,0)) :
    FirstReturn (pairMap (2^(k+1)*L)) (boundary (2^(k+1)*L)) (u,0) (2^(k+1),0) ((k+1)+1) := by
  let d : Int := 2^(k+1)
  let t : Int := 2^k
  let B := d*L
  have ht : 0<t := Int.pow_pos (by decide)
  have hdt : d=2*t := by simp [d,t,Int.pow_succ,Int.mul_comm]
  have hd : 0<d := by omega
  have hbd : d<B := by
    have h := Int.mul_lt_mul_of_pos_left (show 1<L by omega) hd
    simpa [B] using h
  change FirstReturn (pairMap B) (boundary B) (u,0) (d,0) ((k+1)+1)
  rw [first_return_next _ _ _ hp (Or.inl rfl) k, zero_wall, zero_wall]
  change (d,0)=returnTarget B (u,0) (k+1) ∧
    (B∣d*(u-1) ∨ B∣d*u ∨ B∣d*(2*u-1) ∨ B∣d) ∧
    ¬ (B∣t*(u-1) ∨ B∣t*u ∨ B∣t*(2*u-1) ∨ B∣t)
  have ca (x : Int) : B ∣ d*x ↔ L ∣ x := Int.mul_dvd_mul_iff_left (Int.ne_of_gt hd)
  have cb (x : Int) : B ∣ t*x ↔ 2*L ∣ x := by
    change d*L ∣ t*x ↔ _
    rw [hdt, show 2*t*L=t*(2*L) by grind, Int.mul_dvd_mul_iff_left (Int.ne_of_gt ht)]
  have cb1 : B∣t ↔ 2*L∣1 := by simpa using cb 1
  simp only [ca,cb,cb1]
  have eqA : d*(u-1)=B*(2*a+1)+d*(e-1) := by dsimp [B]; grind
  have eqB : d*(B-u)=B*(2*(t-a-1)+1)+(-d*e) := by dsimp [B]; grind
  constructor
  · simp only [returnTarget, pairMap, if_true]
    rw [sortPair_map (fun x => triangle B (d*x)), eqA, eqB]
    rcases he with rfl|rfl
    · simp only [Int.zero_sub, Int.mul_neg, Int.mul_one, Int.mul_zero]
      rw [triangle_at_center _ _ _ (by omega) (by omega), triangle_at_center _ _ _ (by omega) (by omega)]
      grind [absInt,sortPair]
    · simp only [Int.sub_self, Int.mul_zero, Int.mul_one]
      rw [triangle_at_center _ _ _ (by omega) (by omega), triangle_at_center _ _ _ (by omega) (by omega)]
      grind [absInt,sortPair]
  constructor
  · rcases he with rfl|rfl
    · exact Or.inr (Or.inl ⟨2*a+1, by omega⟩)
    · exact Or.inl ⟨2*a+1, by omega⟩
  · have rem (r z : Int) (hz : z=(2*L)*a+r) (hr : 0<r ∧ r<2*L) : ¬2*L∣z := by
      rw [hz, Int.dvd_add_right (Int.dvd_mul_right _ _)]
      intro h
      have := Int.le_of_dvd hr.1 h
      omega
    have one : ¬2*L∣1 := fun h => by have := Int.le_of_dvd (by decide : (0:Int)<1) h; omega
    have hneg : ¬2*L∣u-1 := rem (L+e-1) _ (by grind) (by omega)
    have hnear : ¬2*L∣u := rem (L+e) _ (by grind) (by omega)
    have hdiff : ¬2*L∣2*u-1 := by
      rw [show 2*u-1=(2*L)*(2*a+1)+(2*e-1) by grind, Int.dvd_add_right (Int.dvd_mul_right _ _)]
      rcases he with rfl|rfl
      · exact fun h => one (Int.dvd_neg.mp (by simpa using h))
      · simpa using one
    grind

/-- The zero-boundary involution for powers of two. The row source index is
`t+1` and the return index is `r+1`; their sum is the base exponent. -/
theorem power_zero_return (t r : Nat) :
    let x : Int := 2^(t+1)
    let y : Int := 2^(r+1)
    FirstReturn (pairMap (x*y)) (boundary (x*y)) (x,0) (y,0) ((r+1)+1) := by
  let x : Int := 2^(t+1)
  let y : Int := 2^(r+1)
  have hx : 2≤x := by have : (0:Int)<2^t := Int.pow_pos (by decide); simp [x, Int.pow_succ]; omega
  have hy : 2≤y := by have : (0:Int)<2^r := Int.pow_pos (by decide); simp [y, Int.pow_succ]; omega
  have hb : x < y*x := by
    have := Int.mul_lt_mul_of_pos_right (show 1<y by omega) (show 0<x by omega)
    simpa using this
  have hp : validPair (y*x) (x,0) := ⟨by omega,hb,Int.le_refl 0,by omega⟩
  have hr := zero_block_return r x 0 0 hx (Or.inl rfl) x (by simp) hp
  change FirstReturn (pairMap (x*y)) (boundary (x*y)) (x,0) (y,0) _
  rw [Int.mul_comm x y]; exact hr

theorem power_zero_period (t r : Nat) :
    let x : Int := 2^(t+1)
    let y : Int := 2^(r+1)
    digitPeriod (x*y) (x,0) (if t=r then r+2 else (t+1)+(r+1)+2) := by
  let x : Int := 2^(t+1)
  let y : Int := 2^(r+1)
  have f := power_zero_return t r
  have g := power_zero_return r t
  change FirstReturn (pairMap (y*x)) (boundary (y*x)) _ _ _ at g
  rw [Int.mul_comm y x] at g
  change FirstReturn (pairMap (x*y)) (boundary (x*y)) (x,0) (y,0) _ at f
  have hv : validPair (x*y) (x,0) := by
    have hx : 0<x := Int.pow_pos (by decide)
    have hy : 1<y := by have : (0:Int)<2^r := Int.pow_pos (by decide); simp [y, Int.pow_succ]; omega
    have h := Int.mul_lt_mul_of_pos_left hy hx
    exact ⟨by omega,by simpa using h,Int.le_refl 0,by omega⟩
  dsimp only
  split <;> rename_i he
  · subst r
    have per := first_return_self_exact_period _ _ _ _ f
    exact digitPeriod_of_pair _ _ hv _ per.1 per.2
  · have hne : x≠y := by
      change (2:Int)^(t+1) ≠ 2^(r+1)
      intro h
      rcases (show t+1<r+1 ∨ r+1<t+1 by omega) with hlt|hlt <;>
        have := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) hlt <;> omega
    have per := two_returns_exact_period _ _ _ _ _ _ (show (x,(0:Int))≠(y,0) by grind) f g
    rw [show r+1+1+(t+1+1) = (t+1)+(r+1)+2 by omega] at per
    exact digitPeriod_of_pair _ _ hv _ (by omega) per

/-- The anti-diagonal power-of-two return, with exponent `t+r+2` and
the index involution `t ↦ r`. -/
theorem power_anti_return (t r : Nat) :
    let x : Int := 2^t
    let y : Int := 2^r
    let H := 2*x*y
    FirstReturn (pairMap (2*H)) (boundary (2*H)) (H+x,H-x) (H+y,H-y) (r+1) := by
  let x : Int := 2^t
  have hx : 0<x := Int.pow_pos (by decide)
  cases r with
  | zero =>
    simp only [Int.pow_zero, Int.mul_one]
    change FirstReturn _ _ (2*x+x,2*x-x) (2*x+1,2*x-1) 1
    apply direct_first_return
    · right; right; omega
    · right; right; omega
    · grind [pairMap, sortPair, absInt]
  | succ k =>
    let z : Int := 2^k
    let y : Int := 2^(k+1)
    let H := 2*x*y
    let B := 2*H
    let u := H+x
    let v := H-x
    have hz : 0<z := Int.pow_pos (by decide)
    have hyz : y=2*z := by simp [y,z,Int.pow_succ,Int.mul_comm]
    have hy : 0<y := by omega
    have hHx : x<H := by
      have h := Int.mul_lt_mul_of_pos_left (show 1<2*y by omega) hx
      change x<2*x*y
      grind
    have hHy : y<H := by
      have h := Int.mul_lt_mul_of_pos_left (show 1<2*x by omega) hy
      change y<2*x*y
      grind
    have hv : v≠0 := by dsimp [v]; omega
    have hab : u+v=B := by dsimp [u,v,B]; omega
    have hp : validPair B (u,v) := by dsimp [validPair,u,v,B]; omega
    change FirstReturn (pairMap B) (boundary B) (u,v) (H+y,H-y) ((k+1)+1)
    rw [first_return_next _ _ _ hp (Or.inr (Or.inr hab)) k, diagonal_wall _ _ _ _ hv (Or.inr hab),
      diagonal_wall _ _ _ _ hv (Or.inr hab)]
    constructor
    · have branch : pairMap B (u,v) = (2*x+1,2*x-1) := by
        dsimp [pairMap,u,v,B]
        grind [absInt,sortPair]
      simp only [returnTarget, branch]
      have eqA : y*(2*x+1)=H+y := by dsimp [H]; grind
      have eqB : y*(2*x-1)=H-y := by dsimp [H]; grind
      change (H+y,H-y)=sortPair (triangle B (y*(2*x+1))) (triangle B (y*(2*x-1)))
      rw [eqA,eqB]
      simp only [triangle, Int.emod_eq_of_lt (a:=H+y) (b:=2*B) (by omega) (by dsimp [B];omega),
        Int.emod_eq_of_lt (a:=H-y) (b:=2*B) (by omega) (by dsimp [B];omega)]
      dsimp [B]; grind [absInt,sortPair]
    · constructor
      · right; right; left
        exact ⟨2*y+1, by change y*(4*u)=B*(2*y+1); dsimp [u,B,H]; grind⟩
      · have fact : B=z*(8*x) := by dsimp [B,H]; grind
        have cancel (a : Int) : B ∣ z*a ↔ 8*x ∣ a := by rw [fact, Int.mul_dvd_mul_iff_left (by omega)]
        change ¬ (B∣z*(2*u+1) ∨ B∣z*(2*u-1) ∨ B∣z*(4*u) ∨ B∣2*z)
        rw [Int.mul_comm 2 z]
        simp only [cancel]
        have even (a : Int) (hh : 8*x∣a) : 2∣a :=
          Int.dvd_trans (show (2:Int)∣8*x from ⟨4*x,by grind⟩) hh
        have h4 : ¬8*x∣4*u := by
          rw [show 8*x=(4*x)*2 by grind, show 4*u=(4*x)*(2*y+1) by dsimp [u,H]; grind,
            Int.mul_dvd_mul_iff_left (by omega)]
          omega
        have h2 : ¬8*x∣2 := fun h => by have := Int.le_of_dvd (by decide : (0:Int)<2) h; omega
        intro h
        rcases h with h|h|h|h
        · have := even _ h; omega
        · have := even _ h; omega
        · exact h4 h
        · exact h2 h

theorem power_anti_period (t r : Nat) :
    let x : Int := 2^t
    let y : Int := 2^r
    let H := 2*x*y
    digitPeriod (2*H) (H+x,H-x) (if t=r then r+1 else t+r+2) := by
  let x : Int := 2^t
  let y : Int := 2^r
  let H := 2*x*y
  have hx : 0<x := Int.pow_pos (by decide)
  have hy : 0<y := Int.pow_pos (by decide)
  have hH : x<H := by
    have := Int.mul_lt_mul_of_pos_left (show 1<2*y by omega) hx
    dsimp [H]; grind
  have hv : validPair (2*H) (H+x,H-x) := by dsimp [validPair]; omega
  have f := power_anti_return t r
  have g := power_anti_return r t
  change FirstReturn _ _ (H+x,H-x) (H+y,H-y) _ at f
  change FirstReturn (pairMap (2*(2*y*x))) (boundary (2*(2*y*x)))
    (2*y*x+y,2*y*x-y) (2*y*x+x,2*y*x-x) _ at g
  rw [show 2*y*x=H by dsimp [H]; grind] at g
  dsimp only
  split <;> rename_i he
  · subst t
    have per := first_return_self_exact_period _ _ _ _ f
    exact digitPeriod_of_pair _ _ hv _ per.1 per.2
  · have hxy : x≠y := by
      intro h
      rcases (show t<r ∨ r<t by omega) with hlt|hlt <;>
        have := Int.pow_lt_pow_of_lt (by decide : (1:Int)<2) hlt <;> change 2^t=2^r at h <;> omega
    have per := two_returns_exact_period _ _ _ _ _ _ (show (H+x,H-x)≠(H+y,H-y) by grind) f g
    rw [show r+1+(t+1)=t+r+2 by omega] at per
    exact digitPeriod_of_pair _ _ hv _ (by omega) per

#print axioms power_anti_period
#print axioms power_anti_return
#print axioms power_zero_period
#print axioms zero_block_return
end KaprekarBoundary
