(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 5) x1) (* 0 x2) (* (- 3) x3) ) 6))
(assert (<= (+ (* 2 x1) (* (- 3) x2) (* 0 x3) ) 1))
(assert (<= (+ (* 9 x1) (* (- 7) x2) (* (- 1) x3) ) (- 12)))
(assert (<= (+ (* 7 x1) (* 7 x2) (* (- 10) x3) ) (- 14)))
(assert (<= (+ (* (- 6) x1) (* (- 8) x2) (* 1 x3) ) 34))
(assert (<= (+ (* (- 10) x1) (* 2 x2) (* (- 1) x3) ) (- 15)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 35))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 35))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 35))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 35))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 35))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 35))