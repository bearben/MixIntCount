(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 5) x1) (* (- 3) x2) (* (- 7) x3) ) (- 4)))
(assert (<= (+ (* (- 6) x1) (* (- 2) x2) (* (- 8) x3) ) 35))
(assert (<= (+ (* (- 1) x1) (* 3 x2) (* (- 6) x3) ) (- 5)))
(assert (<= (+ (* 0 x1) (* (- 3) x2) (* (- 2) x3) ) 33))
(assert (<= (+ (* (- 9) x1) (* 8 x2) (* 1 x3) ) 0))
(assert (<= (+ (* 2 x1) (* 9 x2) (* 6 x3) ) 26))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 39))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 39))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 39))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 39))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 39))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 39))