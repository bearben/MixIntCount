(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 4 x1) (* 2 x2) (* (- 2) x3) ) (- 978)))
(assert (<= (+ (* 0 x1) (* (- 9) x2) (* (- 4) x3) ) (- 89)))
(assert (<= (+ (* 6 x1) (* 8 x2) (* (- 3) x3) ) (- 548)))
(assert (<= (+ (* 0 x1) (* 9 x2) (* 1 x3) ) 895))
(assert (<= (+ (* (- 9) x1) (* 7 x2) (* 6 x3) ) (- 747)))
(assert (<= (+ (* (- 4) x1) (* (- 4) x2) (* (- 8) x3) ) (- 520)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 1122))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 1122))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 1122))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 1122))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 1122))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 1122))
