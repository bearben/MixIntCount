(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 0 x1) (* (- 8) x2) (* (- 5) x3) ) 33))
(assert (<= (+ (* (- 5) x1) (* 2 x2) (* (- 2) x3) ) (- 41)))
(assert (<= (+ (* 1 x1) (* 2 x2) (* 2 x3) ) 35))
(assert (<= (+ (* (- 9) x1) (* (- 2) x2) (* (- 1) x3) ) (- 10)))
(assert (<= (+ (* 4 x1) (* 4 x2) (* (- 7) x3) ) (- 25)))
(assert (<= (+ (* (- 3) x1) (* 2 x2) (* (- 10) x3) ) (- 36)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 44))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 44))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 44))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 44))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 44))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 44))