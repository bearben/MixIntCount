(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* (- 9) x1) (* (- 9) x2) (* 9 x3) (* 5 x4) ) 7))
(assert (<= (+ (* (- 4) x1) (* (- 4) x2) (* (- 2) x3) (* 3 x4) ) (- 16)))
(assert (<= (+ (* 9 x1) (* (- 6) x2) (* 9 x3) (* (- 10) x4) ) (- 2)))
(assert (<= (+ (* 4 x1) (* 2 x2) (* 3 x3) (* 3 x4) ) 14))
(assert (<= (+ (* 10 x1) (* 4 x2) (* (- 8) x3) (* (- 8) x4) ) (- 10)))
(assert (<= (+ (* (- 6) x1) (* (- 7) x2) (* 9 x3) (* 10 x4) ) (- 3)))
(assert (<= (+ (* 10 x1) (* 9 x2) (* 3 x3) (* 8 x4) ) (- 12)))
(assert (<= (+ (* (- 5) x1) (* (- 8) x2) (* 7 x3) (* 3 x4) ) 17))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 17))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 17))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 17))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 17))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 17))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 17))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 17))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 17))
