(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 6) x1) (* 5 x2) (* 6 x3) ) 8))
(assert (<= (+ (* (- 9) x1) (* 9 x2) (* 2 x3) ) 0))
(assert (<= (+ (* 2 x1) (* (- 9) x2) (* (- 2) x3) ) 10))
(assert (<= (+ (* (- 5) x1) (* (- 3) x2) (* 1 x3) ) 3))
(assert (<= (+ (* (- 3) x1) (* 7 x2) (* 7 x3) ) (- 5)))
(assert (<= (+ (* 7 x1) (* 7 x2) (* (- 3) x3) ) 11))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 15))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 15))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 15))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 15))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 15))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 15))