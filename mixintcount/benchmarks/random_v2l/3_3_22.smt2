(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 8) x1) (* (- 4) x2) (* 3 x3) ) 4))
(assert (<= (+ (* 2 x1) (* (- 4) x2) (* 3 x3) ) 17))
(assert (<= (+ (* 0 x1) (* 1 x2) (* (- 4) x3) ) (- 6)))
(assert (<= (+ (* 5 x1) (* 1 x2) (* (- 9) x3) ) (- 12)))
(assert (<= (+ (* (- 8) x1) (* (- 6) x2) (* 8 x3) ) (- 21)))
(assert (<= (+ (* 0 x1) (* (- 5) x2) (* 1 x3) ) (- 18)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 22))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 22))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 22))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 22))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 22))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 22))