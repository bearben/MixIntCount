(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 10 x1) (* (- 7) x2) (* (- 8) x3) ) (- 179)))
(assert (<= (+ (* 9 x1) (* (- 10) x2) (* 7 x3) ) 142))
(assert (<= (+ (* (- 1) x1) (* (- 6) x2) (* (- 1) x3) ) 27))
(assert (<= (+ (* (- 6) x1) (* 1 x2) (* 8 x3) ) 21))
(assert (<= (+ (* (- 1) x1) (* 8 x2) (* 4 x3) ) 143))
(assert (<= (+ (* 3 x1) (* 6 x2) (* (- 10) x3) ) 175))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 199))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 199))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 199))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 199))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 199))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 199))