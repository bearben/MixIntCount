(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 1 x1) (* 5 x2) (* 6 x3) (* 4 x4) ) (- 289)))
(assert (<= (+ (* 5 x1) (* (- 10) x2) (* (- 10) x3) (* 9 x4) ) 269))
(assert (<= (+ (* (- 8) x1) (* (- 3) x2) (* (- 2) x3) (* 0 x4) ) 92))
(assert (<= (+ (* 8 x1) (* (- 4) x2) (* (- 8) x3) (* 6 x4) ) (- 270)))
(assert (<= (+ (* (- 10) x1) (* 8 x2) (* (- 10) x3) (* 1 x4) ) 289))
(assert (<= (+ (* (- 7) x1) (* (- 1) x2) (* (- 2) x3) (* 10 x4) ) (- 203)))
(assert (<= (+ (* (- 9) x1) (* 0 x2) (* (- 9) x3) (* 5 x4) ) 327))
(assert (<= (+ (* (- 9) x1) (* 0 x2) (* (- 8) x3) (* (- 8) x4) ) 63))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 354))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 354))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 354))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 354))