(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* (- 9) x1) (* (- 9) x2) (* 0 x3) (* 5 x4) ) 1595))
(assert (<= (+ (* 1 x1) (* 3 x2) (* (- 3) x3) (* 2 x4) ) 260))
(assert (<= (+ (* 1 x1) (* 10 x2) (* 6 x3) (* (- 3) x4) ) 3480))
(assert (<= (+ (* (- 9) x1) (* (- 9) x2) (* 6 x3) (* 5 x4) ) 3959))
(assert (<= (+ (* 6 x1) (* 5 x2) (* (- 1) x3) (* 0 x4) ) 1994))
(assert (<= (+ (* (- 2) x1) (* (- 2) x2) (* (- 2) x3) (* 10 x4) ) (- 1754)))
(assert (<= (+ (* (- 6) x1) (* 9 x2) (* 0 x3) (* 2 x4) ) 2600))
(assert (<= (+ (* (- 3) x1) (* (- 9) x2) (* (- 6) x3) (* 4 x4) ) 3903))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 3981))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 3981))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 3981))
