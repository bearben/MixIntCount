(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 4 x1) (* 8 x2) (* 2 x3) (* 6 x4) ) 160))
(assert (<= (+ (* (- 7) x1) (* 7 x2) (* 8 x3) (* 5 x4) ) (- 227)))
(assert (<= (+ (* 3 x1) (* (- 5) x2) (* 10 x3) (* 4 x4) ) 97))
(assert (<= (+ (* 0 x1) (* (- 9) x2) (* 8 x3) (* 10 x4) ) 114))
(assert (<= (+ (* 4 x1) (* (- 10) x2) (* (- 6) x3) (* 4 x4) ) 36))
(assert (<= (+ (* 3 x1) (* (- 4) x2) (* 8 x3) (* (- 10) x4) ) 50))
(assert (<= (+ (* 7 x1) (* 2 x2) (* (- 7) x3) (* (- 4) x4) ) 4))
(assert (<= (+ (* 9 x1) (* 0 x2) (* (- 8) x3) (* 4 x4) ) (- 119)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 251))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 251))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 251))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 251))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 251))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 251))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 251))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 251))