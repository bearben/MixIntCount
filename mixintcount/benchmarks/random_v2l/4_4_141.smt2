(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 3 x1) (* 0 x2) (* (- 4) x3) (* 7 x4) ) 74))
(assert (<= (+ (* 8 x1) (* 8 x2) (* 2 x3) (* 5 x4) ) 81))
(assert (<= (+ (* 9 x1) (* (- 10) x2) (* (- 7) x3) (* (- 2) x4) ) (- 18)))
(assert (<= (+ (* (- 5) x1) (* 6 x2) (* 10 x3) (* 9 x4) ) 127))
(assert (<= (+ (* 1 x1) (* 3 x2) (* 2 x3) (* 0 x4) ) 117))
(assert (<= (+ (* (- 2) x1) (* 0 x2) (* (- 3) x3) (* 9 x4) ) 90))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 4 x3) (* 5 x4) ) (- 5)))
(assert (<= (+ (* (- 4) x1) (* 2 x2) (* (- 9) x3) (* 8 x4) ) (- 36)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 141))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 141))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 141))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 141))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 141))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 141))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 141))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 141))
