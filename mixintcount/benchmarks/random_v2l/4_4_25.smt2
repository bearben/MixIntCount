(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 1 x1) (* 10 x2) (* 1 x3) (* (- 6) x4) ) 25))
(assert (<= (+ (* (- 2) x1) (* 10 x2) (* 1 x3) (* (- 2) x4) ) (- 17)))
(assert (<= (+ (* 4 x1) (* 1 x2) (* (- 6) x3) (* 0 x4) ) 1))
(assert (<= (+ (* 2 x1) (* 9 x2) (* (- 4) x3) (* 10 x4) ) (- 6)))
(assert (<= (+ (* 10 x1) (* 7 x2) (* 3 x3) (* 0 x4) ) (- 6)))
(assert (<= (+ (* 3 x1) (* (- 9) x2) (* (- 1) x3) (* (- 4) x4) ) (- 18)))
(assert (<= (+ (* 10 x1) (* 5 x2) (* (- 10) x3) (* (- 2) x4) ) (- 10)))
(assert (<= (+ (* 8 x1) (* 7 x2) (* 5 x3) (* (- 2) x4) ) 1))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 25))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 25))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 25))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 25))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 25))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 25))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 25))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 25))
