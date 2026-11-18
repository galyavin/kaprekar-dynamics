import FiveDigit.Core

namespace FiveDigit

structure Digits where
  a : Int
  b : Int
  c : Int
  d : Int
  e : Int
  deriving DecidableEq, Repr

def Digits.valid (B : Int) (q : Digits) : Prop :=
  0 ≤ q.a ∧ q.a < B ∧ 0 ≤ q.b ∧ q.b < B ∧ 0 ≤ q.c ∧ q.c < B ∧
    0 ≤ q.d ∧ q.d < B ∧ 0 ≤ q.e ∧ q.e < B

def Digits.value (B : Int) (q : Digits) : Int :=
  (((q.a * B + q.b) * B + q.c) * B + q.d) * B + q.e

def Digits.rev (q : Digits) : Digits := ⟨q.e,q.d,q.c,q.b,q.a⟩

def Digits.entries (q : Digits) := [q.a,q.b,q.c,q.d,q.e]

/-- Merge the outer and inner pairs, then insert the middle digit. -/
def sortDigits (q : Digits) : Digits :=
  let x := max q.a q.e
  let y := min q.a q.e
  let z := max q.b q.d
  let w := min q.b q.d
  let a := max x z
  let b := max (min x z) (max y w)
  let c := min (min x z) (max y w)
  let d := min y w
  ⟨max a q.c, max b (min a q.c), max c (min b (min a q.c)),
    max d (min c (min b (min a q.c))), min d (min c (min b (min a q.c)))⟩

theorem sort_ordered (q : Digits) :
    (sortDigits q).a ≥ (sortDigits q).b ∧ (sortDigits q).b ≥ (sortDigits q).c ∧
    (sortDigits q).c ≥ (sortDigits q).d ∧ (sortDigits q).d ≥ (sortDigits q).e := by
  grind [sortDigits]

private theorem compare_count (a b : Int) (f : Int → Int) :
    f (max a b) + f (min a b) = f a + f b := by
  grind

/-- Sorting preserves the multiplicity of every digit. -/
theorem sort_multiset (q : Digits) (t : Int) :
    (sortDigits q).entries.count t = q.entries.count t := by
  let f : Int → Int := fun z : Int => if z = t then 1 else 0
  have summation : ∀ q : Digits, (q.entries.count t : Int) =
      f q.a + f q.b + f q.c + f q.d + f q.e := by
    intro q
    simp [Digits.entries, List.count_cons, f, BEq.beq, decide_eq_true_eq]
    split <;> split <;> split <;> split <;> split <;> grind
  have h := compare_count q.a q.e f
  have i := compare_count q.b q.d f
  let x := max q.a q.e
  let y := min q.a q.e
  let z := max q.b q.d
  let w := min q.b q.d
  have j := compare_count x z f
  have k := compare_count y w f
  have l := compare_count (min x z) (max y w) f
  let a := max x z
  let b := max (min x z) (max y w)
  let c := min (min x z) (max y w)
  let d := min y w
  have r := compare_count a q.c f
  have s := compare_count b (min a q.c) f
  have v := compare_count c (min b (min a q.c)) f
  have p := compare_count d (min c (min b (min a q.c))) f
  have hq := summation q
  have hs := summation (sortDigits q)
  change (_ : Int) = f (max a q.c) + f (max b (min a q.c)) +
    f (max c (min b (min a q.c))) + f (max d (min c (min b (min a q.c)))) + f (min d (min c (min b (min a q.c)))) at hs
  dsimp only [a,b,c,d,x,y,z,w] at *
  omega

def differences (q : Digits) : Pair :=
  ((sortDigits q).a - (sortDigits q).e, (sortDigits q).b - (sortDigits q).d)

def borrowed (B : Int) (p : Pair) : Digits :=
  if p.2 = 0 then ⟨p.1-1,B-1,B-1,B-1,B-p.1⟩
  else ⟨p.1,p.2-1,B-1,B-1-p.2,B-p.1⟩

def subtraction (B : Int) (q : Digits) : Digits :=
  if (differences q).1 = 0 then ⟨0,0,0,0,0⟩ else borrowed B (differences q)

theorem subtraction_value (B : Int) (q : Digits) :
    (subtraction B q).value B = (sortDigits q).value B - (sortDigits q).rev.value B := by
  have ho := sort_ordered q
  simp only [subtraction, differences, borrowed, Digits.value, Digits.rev]
  generalize sortDigits q = z at *
  rcases z with ⟨a,b,c,d,e⟩
  dsimp at *
  split
  · have : a = e := by omega
    have : b = e := by omega
    have : c = e := by omega
    have : d = e := by omega
    grind
  · split <;> grind

theorem differences_valid (B : Int) (q : Digits) (hv : q.valid B)
    (hn : (differences q).1 ≠ 0) : valid B (differences q) := by
  have ho := sort_ordered q
  have hs : (sortDigits q).valid B := by
    unfold Digits.valid sortDigits at *
    dsimp
    omega
  unfold differences at *
  unfold Digits.valid at hs
  unfold valid
  dsimp at *
  omega

theorem borrowed_valid (B : Int) (q : Pair) (hq : valid B q) : (borrowed B q).valid B := by
  grind [Digits.valid, valid, borrowed]

theorem subtraction_valid (B : Int) (hB : 2 ≤ B) (q : Digits) (hv : q.valid B) :
    (subtraction B q).valid B := by
  unfold subtraction
  split
  · simp [Digits.valid]; omega
  · exact borrowed_valid B _ (differences_valid B q hv (by assumption))

theorem borrowed_pair (B : Int) (p : Pair) (hp : valid B p) :
    differences (borrowed B p) = pairStep B p := by
  rcases p with ⟨u,v⟩
  obtain ⟨_,_,_,_⟩ := hp
  grind [differences, sortDigits, borrowed, pairStep, positiveStep]

theorem subtraction_pair (B : Int) (q : Digits) (hv : q.valid B)
    (hn : (differences q).1 ≠ 0) :
    differences (subtraction B q) = pairStep B (differences q) := by
  rw [subtraction, if_neg hn]
  exact borrowed_pair B _ (differences_valid B q hv hn)

/-- Zero outer difference characterizes words with five equal digits. -/
theorem constant_iff (q : Digits) :
    (differences q).1 = 0 ↔ q.a = q.b ∧ q.b = q.c ∧ q.c = q.d ∧ q.d = q.e := by
  grind (splits := 30) [differences, sortDigits]

/-- The diagram commutes at every iterate from a nonconstant input, including
words whose most significant digit is zero. -/
theorem iterated_subtraction (B : Int) (hB : 2 ≤ B) (n : Nat) (q : Digits)
    (hv : q.valid B) (hn : (differences q).1 ≠ 0) :
    differences (iter (subtraction B) n q) = iter (pairStep B) n (differences q) := by
  induction n generalizing q with
  | zero => rfl
  | succ n ih =>
    have stepEq := subtraction_pair B q hv hn
    have hp := (pair_valid B _ (differences_valid B q hv hn)).1
    have nextNonzero : (differences (subtraction B q)).1 ≠ 0 := by rw [stepEq]; omega
    exact (ih _ (subtraction_valid B hB q hv) nextNonzero).trans (congrArg _ stepEq)

#print axioms constant_iff
#print axioms iterated_subtraction
#print axioms sort_ordered
#print axioms sort_multiset
#print axioms subtraction_value
#print axioms subtraction_valid
#print axioms borrowed_pair
#print axioms subtraction_pair
end FiveDigit
