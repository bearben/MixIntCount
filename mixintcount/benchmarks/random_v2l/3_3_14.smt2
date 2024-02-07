(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 3) x1) (* (- 6) x2) (* 4 x3) ) 10))
(assert (<= (+ (* (- 1) x1) (* (- 6) x2) (* (- 10) x3) ) 1))
(assert (<= (+ (* 2 x1) (* (- 1) x2) (* (- 6) x3) ) 11))
(assert (<= (+ (* 2 x1) (* 0 x2) (* 4 x3) ) 4))
(assert (<= (+ (* 1 x1) (* (- 6) x2) (* (- 4) x3) ) (- 4)))
(assert (<= (+ (* (- 3) x1) (* (- 9) x2) (* (- 6) x3) ) (- 4)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 14))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 14))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 14))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 14))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 14))
