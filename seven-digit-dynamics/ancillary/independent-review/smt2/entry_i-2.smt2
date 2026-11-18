(set-logic QF_LRA)
; benchmark generated from python API
(set-info :status unknown)
(declare-fun c () Real)
(declare-fun b () Real)
(declare-fun a () Real)
(assert
 (>= c 0.0))
(assert
 (<= c b))
(assert
 (<= b a))
(assert
 (< a 1.0))
(assert
 (> (/ 1.0 2.0) a))
(assert
 (>= (+ (* 2.0 a) b) 1.0))
(assert
 (let ((?x169 (- 1.0 a)))
(let ((?x170 (- ?x169 b)))
(let ((?x8 (+ a b)))
(let ((?x10 (- ?x8 1.0)))
(let ((?x14 (ite (>= ?x10 0.0) ?x10 (- ?x10))))
(let (($x20 (>= a ?x14)))
(let ((?x24 (ite $x20 ?x14 a)))
(let ((?x15 (+ b c)))
(let ((?x16 (- ?x15 1.0)))
(let ((?x19 (ite (>= ?x16 0.0) ?x16 (- ?x16))))
(let ((?x26 (ite (>= ?x24 ?x19) ?x19 ?x24)))
(let ((?x137 (- (+ ?x8 c) 1.0)))
(let ((?x21 (ite $x20 a ?x14)))
(let (($x22 (>= ?x21 ?x19)))
(let ((?x27 (ite $x22 ?x19 ?x21)))
(let ((?x29 (ite (>= ?x24 ?x27) ?x24 ?x27)))
(let ((?x23 (ite $x22 ?x21 ?x19)))
(or (and (distinct ?x23 (+ a (ite (>= (- ?x137) 0.0) (- ?x137) 0.0))) true) (and (distinct ?x29 (- a (ite (>= ?x137 0.0) ?x137 0.0))) true) (and (distinct ?x26 ?x170) true))))))))))))))))))))
(check-sat)
