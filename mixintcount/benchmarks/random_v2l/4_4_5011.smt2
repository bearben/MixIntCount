(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* (- 2) x1) (* 4 x2) (* 1 x3) (* 3 x4) ) (- 1998)))
(assert (<= (+ (* (- 10) x1) (* (- 10) x2) (* 5 x3) (* (- 10) x4) ) (- 1664)))
(assert (<= (+ (* 1 x1) (* (- 7) x2) (* 9 x3) (* 1 x4) ) (- 2251)))
(assert (<= (+ (* (- 3) x1) (* (- 1) x2) (* (- 2) x3) (* 5 x4) ) (- 2810)))
(assert (<= (+ (* 9 x1) (* (- 6) x2) (* (- 8) x3) (* (- 1) x4) ) 2410))
(assert (<= (+ (* 0 x1) (* 6 x2) (* 0 x3) (* 1 x4) ) 3818))
(assert (<= (+ (* (- 8) x1) (* 8 x2) (* 10 x3) (* 3 x4) ) 2558))
(assert (<= (+ (* (- 9) x1) (* 1 x2) (* 0 x3) (* 5 x4) ) 3666))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 5011))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 5011))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 5011))
