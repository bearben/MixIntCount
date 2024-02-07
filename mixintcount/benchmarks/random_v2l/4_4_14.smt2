(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* (- 9) x1) (* (- 6) x2) (* (- 7) x3) (* (- 1) x4) ) (- 1)))
(assert (<= (+ (* (- 1) x1) (* 7 x2) (* (- 9) x3) (* (- 9) x4) ) 11))
(assert (<= (+ (* 6 x1) (* (- 6) x2) (* (- 9) x3) (* (- 8) x4) ) 7))
(assert (<= (+ (* (- 9) x1) (* 7 x2) (* (- 1) x3) (* (- 7) x4) ) 14))
(assert (<= (+ (* 1 x1) (* (- 7) x2) (* (- 5) x3) (* (- 3) x4) ) 1))
(assert (<= (+ (* (- 5) x1) (* (- 10) x2) (* 7 x3) (* 8 x4) ) 10))
(assert (<= (+ (* 2 x1) (* 8 x2) (* 5 x3) (* 5 x4) ) (- 1)))
(assert (<= (+ (* (- 4) x1) (* (- 8) x2) (* (- 9) x3) (* (- 5) x4) ) (- 10)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 14))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 14))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 14))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 14))
