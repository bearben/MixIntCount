(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 1 x1) (* 4 x2) (* 3 x3) ) (- 124)))
(assert (<= (+ (* 1 x1) (* 2 x2) (* 2 x3) ) 118))
(assert (<= (+ (* (- 5) x1) (* (- 9) x2) (* (- 3) x3) ) (- 9)))
(assert (<= (+ (* (- 10) x1) (* 1 x2) (* 6 x3) ) 143))
(assert (<= (+ (* (- 4) x1) (* 0 x2) (* (- 7) x3) ) (- 109)))
(assert (<= (+ (* (- 3) x1) (* (- 10) x2) (* 7 x3) ) (- 47)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 177))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 177))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 177))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 177))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 177))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 177))
