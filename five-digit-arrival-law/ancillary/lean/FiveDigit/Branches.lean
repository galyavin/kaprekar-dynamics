import FiveDigit.Core

namespace FiveDigit

/-- Under the state-domain inequalities only nine of the 24 positive-coordinate
orders can occur. The names of the positive cases give the digit order on
(u, v-1, B-1-v, B-u), with ties resolved by increasing index. -/
theorem step_cases (B u v : Int) (hp : valid B (u,v)) (P : Pair → Prop)
    (p0123 : 0 < v → B ≤ 2*v → v < u → P (u-1,u+v-B+1))
    (p0132 : 0 < v → B < u+v → u = v → P (v,2*u-B))
    (p0213 : 0 < v → 2*v < B → B < u+v → P (u-1,u-v+1))
    (p0231 : 0 < v → B-1 ≤ u+v → v < u → u+v ≤ B → P (B-v,2*u-B))
    (p0312 : 0 < v → B ≤ 2*u → u+v ≤ B → B ≤ 2*v → P (v,u-v+1))
    (p2031 : 0 < v → u+v ≤ B-2 → B ≤ 2*u → P (B-v,u-v-1))
    (p2301 : 0 < v → v < u → 2*u < B → P (B-v,B-u-v-1))
    (p3021 : 0 < v → 2*u < B → B-1 ≤ u+v → 2*v < B → P (B-v,v-u+1))
    (p3201 : 0 < v → u = v → u+v ≤ B-2 → P (B-v,B-2*u))
    (z0 : v = 0 → B < 2*u → P (u-1,B-u))
    (z1 : v = 0 → 2*u ≤ B → P (B-u,u-1)) : P (pairStep B (u,v)) := by
  obtain ⟨hu,huB,hv,hvu⟩ := hp
  dsimp at *
  by_cases hz : v = 0
  · by_cases hi : B < 2*u
    · have := z0 hz hi
      grind [pairStep]
    · have := z1 hz (by omega)
      grind [pairStep]
  · have pos : 0 < v := by omega
    by_cases hh : B ≤ 2*v
    · by_cases hd : v < u
      · have := p0123 pos hh hd
        grind [pairStep, positiveStep]
      · by_cases ht : B < u+v
        · have := p0132 pos ht (by omega)
          grind [pairStep, positiveStep]
        · have := p0312 pos (by omega) (by omega) hh
          grind [pairStep, positiveStep]
    · by_cases hu2 : 2*u < B
      · by_cases hd : v < u
        · have := p2301 pos hd hu2
          grind [pairStep, positiveStep]
        · by_cases ht : B-1 ≤ u+v
          · have := p3021 pos hu2 ht (by omega)
            grind [pairStep, positiveStep]
          · have := p3201 pos (by omega) (by omega)
            grind [pairStep, positiveStep]
      · by_cases ht : B < u+v
        · have := p0213 pos (by omega) ht
          grind [pairStep, positiveStep]
        · by_cases hb : B-1 ≤ u+v
          · have := p0231 pos hb (by omega) (by omega)
            grind [pairStep, positiveStep]
          · have := p2031 pos (by omega) (by omega)
            grind [pairStep, positiveStep]

#print axioms step_cases
end FiveDigit
