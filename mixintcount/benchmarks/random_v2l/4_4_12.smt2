(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 6 x1) (* (- 10) x2) (* (- 6) x3) (* (- 10) x4) ) 11))
(assert (<= (+ (* (- 1) x1) (* (- 8) x2) (* 3 x3) (* 10 x4) ) (- 6)))
(assert (<= (+ (* (- 4) x1) (* (- 8) x2) (* 5 x3) (* 10 x4) ) (- 3)))
(assert (<= (+ (* (- 7) x1) (* 10 x2) (* 7 x3) (* 8 x4) ) 8))
(assert (<= (+ (* 5 x1) (* 8 x2) (* 10 x3) (* (- 6) x4) ) (- 5)))
(assert (<= (+ (* (- 4) x1) (* 6 x2) (* 9 x3) (* (- 10) x4) ) 5))
(assert (<= (+ (* (- 2) x1) (* 6 x2) (* 3 x3) (* 2 x4) ) 0))
(assert (<= (+ (* 5 x1) (* 9 x2) (* 6 x3) (* (- 5) x4) ) (- 9)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 12))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 12))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 12))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 12))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 12))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 12))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 12))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 12))