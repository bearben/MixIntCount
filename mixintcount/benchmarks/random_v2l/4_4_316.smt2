(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 8 x1) (* 9 x2) (* (- 8) x3) (* 3 x4) ) (- 285)))
(assert (<= (+ (* 3 x1) (* (- 10) x2) (* 8 x3) (* (- 9) x4) ) (- 3)))
(assert (<= (+ (* (- 9) x1) (* 0 x2) (* (- 4) x3) (* (- 2) x4) ) 177))
(assert (<= (+ (* 4 x1) (* (- 10) x2) (* 4 x3) (* (- 3) x4) ) (- 126)))
(assert (<= (+ (* (- 3) x1) (* (- 5) x2) (* (- 3) x3) (* (- 1) x4) ) 287))
(assert (<= (+ (* (- 10) x1) (* 2 x2) (* 9 x3) (* (- 8) x4) ) 199))
(assert (<= (+ (* 5 x1) (* 10 x2) (* (- 8) x3) (* 7 x4) ) 316))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* (- 1) x3) (* (- 2) x4) ) 30))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 316))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 316))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 316))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 316))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 316))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 316))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 316))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 316))