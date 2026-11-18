import Boundary

namespace KaprekarBoundary

/-- Four positional digits, most significant first; leading zeroes are allowed. -/
structure Digits where
  a : Int
  b : Int
  c : Int
  d : Int
  deriving DecidableEq, Repr

def Digits.valid (B : Int) (q : Digits) : Prop :=
  (0 ≤ q.a ∧ q.a < B) ∧ (0 ≤ q.b ∧ q.b < B) ∧
  (0 ≤ q.c ∧ q.c < B) ∧ (0 ≤ q.d ∧ q.d < B)

def Digits.value (B : Int) (q : Digits) : Int :=
  ((q.a * B + q.b) * B + q.c) * B + q.d

def Digits.rev (q : Digits) : Digits := ⟨q.d,q.c,q.b,q.a⟩

/-- Sort two pairs, then compare their upper and lower entries. The middle two
entries are the smaller upper entry and the larger lower entry, in sorted order. -/
def sort4 (q : Digits) : Digits :=
  let (i,j) := sortPair q.a q.b
  let (k,l) := sortPair q.c q.d
  let (a,x) := sortPair i k
  let (y,d) := sortPair j l
  let (b,c) := sortPair x y
  ⟨a,b,c,d⟩

def diff (q : Digits) : Pair := ((sort4 q).a - (sort4 q).d, (sort4 q).b - (sort4 q).c)

def Digits.nonconstant (q : Digits) : Prop := ¬(q.a=q.b ∧ q.b=q.c ∧ q.c=q.d)

theorem nonconstant_iff (q : Digits) : q.nonconstant ↔ (diff q).1≠0 := by
  grind (splits := 60) [Digits.nonconstant, diff, sort4, sortPair]

def validPair (B : Int) (p : Pair) : Prop := 1 ≤ p.1 ∧ p.1 < B ∧ 0 ≤ p.2 ∧ p.2 ≤ p.1

/-- Four digits produced by borrowing in descending minus ascending order. -/
def fromPair (B : Int) (p : Pair) : Digits :=
  if p.2 = 0 then ⟨p.1-1,B-1,B-1,B-p.1⟩
  else ⟨p.1,p.2-1,B-1-p.2,B-p.1⟩

def step (B : Int) (q : Digits) : Digits :=
  if (diff q).1 = 0 then ⟨0,0,0,0⟩ else fromPair B (diff q)

theorem sort4_ordered (q : Digits) :
    (sort4 q).a ≥ (sort4 q).b ∧ (sort4 q).b ≥ (sort4 q).c ∧
    (sort4 q).c ≥ (sort4 q).d := by
  grind [sort4, sortPair]

theorem pair_of_subtraction (B : Int) (p : Pair) (hp : validPair B p) :
    diff (fromPair B p) = pairMap B p := by
  rcases p with ⟨u,v⟩
  obtain ⟨hu,hB,hv,hvu⟩ := hp
  by_cases hz : v = 0
  · subst v
    simp only [diff, fromPair, pairMap, if_pos, sort4]
    grind (splits := 60) [sortPair]
  · simp only [diff, fromPair, pairMap, if_neg hz, sort4, genericMap]
    grind (splits := 60) [sortPair, absInt]


theorem sort4_valid (B : Int) (q : Digits) (hq : q.valid B) : (sort4 q).valid B := by
  grind [Digits.valid, sort4, sortPair]

/-- Characterize a digit multiset by its count of each possible entry. -/
def Digits.count (q : Digits) (x : Int) : Nat :=
  (if q.a = x then 1 else 0) + (if q.b = x then 1 else 0) +
  (if q.c = x then 1 else 0) + (if q.d = x then 1 else 0)

theorem sort4_count (q : Digits) (x : Int) : (sort4 q).count x = q.count x := by
  unfold Digits.count
  generalize (fun a : Int => if a = x then 1 else 0 : Int → Nat) = f
  grind (splits := 60) [sort4, sortPair]

theorem diff_valid (B : Int) (q : Digits) (hq : q.valid B) (hne : (diff q).1 ≠ 0) :
    validPair B (diff q) := by
  have ho := sort4_ordered q
  have hv := sort4_valid B q hq
  unfold Digits.valid at hv
  unfold diff at *
  unfold validPair
  dsimp
  omega

theorem fromPair_valid (B : Int) (p : Pair) (hp : validPair B p) : (fromPair B p).valid B := by
  grind [fromPair, validPair, Digits.valid]

theorem step_value (B : Int) (q : Digits) :
    (step B q).value B = (sort4 q).value B - (sort4 q).rev.value B := by
  have ho := sort4_ordered q
  simp only [step, diff, fromPair, Digits.value, Digits.rev]
  generalize sort4 q = z at *
  rcases z with ⟨a,b,c,d⟩
  dsimp at *
  split
  · have : a = d := by omega
    have : b = d := by omega
    have : c = d := by omega
    grind
  · split <;> grind

theorem step_valid (B : Int) (hB : 2 ≤ B) (q : Digits) (hq : q.valid B) :
    (step B q).valid B := by
  unfold step
  split
  · simp [Digits.valid]; omega
  · apply fromPair_valid; exact diff_valid B q hq (by assumption)

theorem pairMap_valid (B : Int) (p : Pair) (hp : validPair B p) :
    validPair B (pairMap B p) := by
  rcases p with ⟨u,v⟩
  obtain ⟨hu, hB, hv, hvu⟩ := hp
  grind (splits := 60) [pairMap, validPair, sortPair, genericMap, absInt]

theorem fromPair_injective (B : Int) (p t : Pair) (hp : validPair B p)
    (ht : validPair B t) (he : fromPair B p = fromPair B t) : p = t := by
  rcases p with ⟨u,v⟩
  rcases t with ⟨w,z⟩
  have hlast := congrArg Digits.d he
  have hsecond := congrArg Digits.b he
  grind [validPair, fromPair]

theorem step_fromPair (B : Int) (p : Pair) (hp : validPair B p) :
    step B (fromPair B p) = fromPair B (pairMap B p) := by
  have hvalid := pairMap_valid B p hp
  simp only [step, pair_of_subtraction B p hp]
  rw [if_neg (show (pairMap B p).1 ≠ 0 by have := hvalid.1; omega)]

theorem iter_pair_valid (B : Int) (p : Pair) (hp : validPair B p) (k : Nat) :
    validPair B (iter (pairMap B) k p) := by
  induction k with
  | zero => exact hp
  | succ k ih => exact pairMap_valid B _ ih

theorem iter_fromPair (B : Int) (p : Pair) (hp : validPair B p) (k : Nat) :
    iter (step B) k (fromPair B p) = fromPair B (iter (pairMap B) k p) := by
  induction k with
  | zero => rfl
  | succ k ih => exact (congrArg (step B) ih).trans (step_fromPair B _ (iter_pair_valid B p hp k))

theorem digit_return_iff (B : Int) (p : Pair) (hp : validPair B p) (k : Nat) :
    iter (step B) k (fromPair B p) = fromPair B p ↔ iter (pairMap B) k p = p := by
  rw [iter_fromPair B p hp]
  exact ⟨fromPair_injective B _ p (iter_pair_valid B p hp k) hp, congrArg (fromPair B)⟩

/-- Every nonconstant periodic word is in the image of exactly one legal pair. -/
theorem periodic_digit_has_pair (B : Int) (q : Digits) (hq : q.valid B)
    (hne : (diff q).1 ≠ 0) (k : Nat) (hk : iter (step B) (k+1) q = q) :
    ∃ p, validPair B p ∧ fromPair B p = q ∧
      ∀ t, validPair B t → fromPair B t = q → t = p := by
  have hp := diff_valid B q hq hne
  have hstep : step B q = fromPair B (diff q) := by simp [step, hne]
  have he : fromPair B (iter (pairMap B) k (diff q)) = q := by
    rw [← iter_fromPair B _ hp, ← hstep]
    change iter (step B) k (iter (step B) 1 q) = q
    rw [← iter_add, Nat.add_comm]
    exact hk
  exact ⟨_, iter_pair_valid B _ hp k, he, fun t ht et =>
    fromPair_injective B t _ ht (iter_pair_valid B _ hp k) (et.trans he.symm)⟩

/-- The two nonconstant binary cycles, stated on four-slot words. -/
theorem base_two_cycles (q : Digits) (hq : q.valid 2) (hne : (diff q).1 ≠ 0) :
    (∃ k : Nat, iter (step 2) (k+1) q = q) ↔
      q = ⟨0,1,1,1⟩ ∨ q = ⟨1,0,0,1⟩ := by
  constructor
  · rintro ⟨k, hk⟩
    obtain ⟨⟨u,v⟩, hp, he, _⟩ := periodic_digit_has_pair 2 q hq hne k hk
    grind [validPair, fromPair]
  · rintro (rfl | rfl) <;> exact ⟨0, by decide⟩

#print axioms base_two_cycles
#print axioms sort4_count
#print axioms pair_of_subtraction
#print axioms step_value
#print axioms digit_return_iff
#print axioms periodic_digit_has_pair
end KaprekarBoundary
