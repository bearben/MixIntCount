(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 6 x1) (* (- 6) x2) (* 0 x3) ) 191))
(assert (<= (+ (* (- 8) x1) (* (- 6) x2) (* 2 x3) ) 153))
(assert (<= (+ (* (- 3) x1) (* (- 4) x2) (* 6 x3) ) (- 330)))
(assert (<= (+ (* (- 8) x1) (* (- 2) x2) (* 4 x3) ) 196))
(assert (<= (+ (* 7 x1) (* (- 10) x2) (* (- 6) x3) ) 224))
(assert (<= (+ (* (- 7) x1) (* (- 7) x2) (* (- 3) x3) ) (- 142)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 354))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 354))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 354))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 354))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 354))
