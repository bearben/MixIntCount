(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 5) x1) (* (- 7) x2) (* (- 8) x3) ) (- 1932)))
(assert (<= (+ (* 10 x1) (* (- 8) x2) (* 8 x3) ) (- 3470)))
(assert (<= (+ (* 2 x1) (* (- 2) x2) (* 8 x3) ) (- 754)))
(assert (<= (+ (* 7 x1) (* 10 x2) (* (- 7) x3) ) 2404))
(assert (<= (+ (* 5 x1) (* 6 x2) (* 9 x3) ) 1982))
(assert (<= (+ (* 2 x1) (* 1 x2) (* 10 x3) ) 5095))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 5623))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 5623))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 5623))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 5623))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 5623))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 5623))