import Std

namespace FiveDigit

abbrev Pair := Int × Int

def valid (B : Int) (q : Pair) : Prop :=
  1 ≤ q.1 ∧ q.1 < B ∧ 0 ≤ q.2 ∧ q.2 ≤ q.1

def margins (B : Int) (q : Pair) : Pair :=
  (q.1 + q.2 - B, 2*q.1 - q.2 - B)

def core (B : Int) (q : Pair) : Prop :=
  valid B q ∧ 2*q.2 ≤ B-1 ∧ 2*(q.1-q.2) ≤ B-3 ∧
    -2 ≤ (margins B q).1 ∧ -2 ≤ (margins B q).2

def strip (B : Int) (q : Pair) : Prop :=
  valid B q ∧ B-1 ≤ 2*q.1 ∧ 2*q.2 ≤ B+1 ∧ 2*(q.1-q.2) ≤ B+2

def marginStep (q : Pair) : Pair :=
  if 1 ≤ q.1 then (q.2, q.1-3) else
  if -1 ≤ q.1 then (q.2, -2*q.1) else (q.2+1, 3)

def height (q : Pair) : Int := max q.1 q.2

def rare (m : Int) (q : Pair) : Prop :=
  q = (m,m-1) ∨ q = (m+1,0) ∨ q = (m,0)

def twoCycle (B : Int) (q : Pair) : Prop :=
  (2*q.1 = B ∧ q.2 = q.1-1) ∨ (2*(q.1-1) = B ∧ q.2 = 0)

/-- The largest digit is B-1. Merge two sorted pairs of the other four
digits. Only their maximum, minimum and third-largest value are needed. -/
def positiveStep (B u v : Int) : Pair :=
  let x := max u (B-u)
  let y := min u (B-u)
  let z := max (v-1) (B-1-v)
  let w := min (v-1) (B-1-v)
  (B-1-min y w, max x z - min (min x z) (max y w))

def pairStep (B : Int) (q : Pair) : Pair :=
  if q.2 = 0 then
    (B-1-min (q.1-1) (B-q.1), B-1-max (q.1-1) (B-q.1))
  else positiveStep B q.1 q.2

def iter {α : Type} (f : α → α) : Nat → α → α
  | 0, q => q
  | n+1, q => iter f n (f q)

end FiveDigit
