(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 3) x1) (* 8 x2) (* (- 9) x3) ) (- 227)))
(assert (<= (+ (* (- 4) x1) (* (- 8) x2) (* (- 4) x3) ) (- 1718)))
(assert (<= (+ (* (- 10) x1) (* (- 2) x2) (* (- 4) x3) ) (- 105)))
(assert (<= (+ (* 8 x1) (* 2 x2) (* 1 x3) ) 1127))
(assert (<= (+ (* 6 x1) (* 3 x2) (* (- 2) x3) ) (- 940)))
(assert (<= (+ (* 4 x1) (* 2 x2) (* 6 x3) ) 1842))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 1995))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 1995))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 1995))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 1995))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 1995))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 1995))