import FiveDigit.Branches

namespace FiveDigit

/-- The two integer inequalities in (36), after multiplication by two. -/
def marginBound (B initial : Int) (p : Pair) : Prop :=
  2*(margins B p).1 ≤ 2*initial-B+40 ∧ 2*(margins B p).2 ≤ 2*initial-B+40

/-- A certificate stops at the first core entry or at an exceptional cycle
state. In the high starting range only core entry is possible. -/
private inductive Entry (B initial : Int) : Nat → Pair → Prop where
  | here {n p} : core B p → (3*B ≤ 4*initial → marginBound B initial p) → Entry B initial n p
  | exceptional {n p} : twoCycle B p → 4*initial < 3*B → Entry B initial n p
  | next {n p} : ¬ core B p → Entry B initial n (pairStep B p) → Entry B initial (n+1) p

private theorem entry_time {B initial n p} (h : Entry B initial n p) :
    ∃ j ≤ n, core B (iter (pairStep B) j p) ∨ twoCycle B (iter (pairStep B) j p) := by
  induction h with
  | here hc _ => exact ⟨0, Nat.zero_le _, Or.inl hc⟩
  | exceptional ht _ => exact ⟨0, Nat.zero_le _, Or.inr ht⟩
  | next _ _ ih =>
    obtain ⟨j,hj,hit⟩ := ih
    exact ⟨j+1, by omega, hit⟩

private theorem high_capture {B initial n p} (h : Entry B initial n p)
    (high : 3*B ≤ 4*initial) : ∃ j ≤ n,
    core B (iter (pairStep B) j p) ∧
    (∀ i < j, ¬ core B (iter (pairStep B) i p)) ∧
    marginBound B initial (iter (pairStep B) j p) := by
  induction h with
  | here hc bound => exact ⟨0, Nat.zero_le _, hc, by simp, bound high⟩
  | exceptional _ low => omega
  | next notCore _ ih =>
    obtain ⟨j,hj,hc,early,bound⟩ := ih
    refine ⟨j+1, by omega, hc, ?_, bound⟩
    intro i hi
    cases i with
    | zero => exact notCore
    | succ i => exact early i (by omega)

private theorem entry_bound {B initial n p} (h : Entry B initial n p)
    (high : 3*B ≤ 4*initial) (j : Nat) (hj : j ≤ n)
    (hc : core B (iter (pairStep B) j p))
    (first : ∀ i < j, ¬ core B (iter (pairStep B) i p)) :
    marginBound B initial (iter (pairStep B) j p) := by
  induction h generalizing j with
  | here hcore hb =>
    cases j with
    | zero => exact hb high
    | succ j => exact False.elim (first 0 (by omega) hcore)
  | exceptional _ hl => omega
  | next hnc _ ih =>
    cases j with
    | zero => exact False.elim (hnc hc)
    | succ j => exact ih j (by omega) hc (fun i hi => first (i+1) (by omega))


private theorem split_entry {B initial n p}
    (hc : core B p → 3*B ≤ 4*initial → marginBound B initial p)
    (ht : twoCycle B p → 4*initial < 3*B)
    (advance : ¬ core B p → ¬ twoCycle B p → Entry B initial n p) :
    Entry B initial n p := by
  by_cases h : core B p
  · exact .here h (hc h)
  · by_cases t : twoCycle B p
    · exact .exceptional t (ht t)
    · exact advance h t

-- Split the affine branches one step at a time. Contradictory guards and
-- terminal cases are discharged over integers with no upper bound on B.
local macro "linear" : tactic =>
  `(tactic| (simp only [strip, core, valid, margins, twoCycle, marginBound,
      Prod.fst, Prod.snd] at * <;> omega))

local macro "takeStep" : tactic =>
  `(tactic|
    first
    | exfalso; linear
    | (apply split_entry
       · intros; linear
       · intro; linear
       intro hc ht
       first
       | exfalso; linear
       | (apply Entry.next hc
          apply step_cases
          · linear
          all_goals intros)))
set_option maxHeartbeats 16000000 in
set_option maxRecDepth 16000 in
private theorem certificate (B u v : Int) (hB : 30 ≤ B) (hs : strip B (u,v)) :
    Entry B u 6 (u,v) := by
  iterate 7 (all_goals takeStep)

theorem six_step_entry (B : Int) (hB : 30 ≤ B) (p : Pair) (hp : strip B p) :
    ∃ j ≤ 6, core B (iter (pairStep B) j p) ∨ twoCycle B (iter (pairStep B) j p) := by
  exact entry_time (certificate B p.1 p.2 hB hp)

theorem high_first_entry (B : Int) (hB : 30 ≤ B) (p : Pair) (hp : strip B p)
    (high : 3*B ≤ 4*p.1) : ∃ j ≤ 6,
    core B (iter (pairStep B) j p) ∧
    (∀ i < j, ¬ core B (iter (pairStep B) i p)) ∧
    marginBound B p.1 (iter (pairStep B) j p) :=
  high_capture (certificate B p.1 p.2 hB hp) high

theorem first_entry_bound (B : Int) (hB : 30 ≤ B) (p : Pair) (hp : strip B p)
    (high : 3*B ≤ 4*p.1) (j : Nat) (hj : j ≤ 6)
    (hc : core B (iter (pairStep B) j p))
    (first : ∀ i < j, ¬ core B (iter (pairStep B) i p)) :
    marginBound B p.1 (iter (pairStep B) j p) := by
  exact entry_bound (certificate B p.1 p.2 hB hp) high j hj hc first

#print axioms six_step_entry
#print axioms high_first_entry
#print axioms first_entry_bound
end FiveDigit
