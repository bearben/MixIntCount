(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* (- 1) x1) (* 7 x2) (* 6 x3) ) 222))
(assert (<= (+ (* 9 x1) (* (- 9) x2) (* (- 6) x3) ) 366))
(assert (<= (+ (* 8 x1) (* 8 x2) (* (- 10) x3) ) (- 12)))
(assert (<= (+ (* (- 1) x1) (* 8 x2) (* 3 x3) ) 474))
(assert (<= (+ (* (- 2) x1) (* 4 x2) (* (- 3) x3) ) 341))
(assert (<= (+ (* (- 3) x1) (* (- 2) x2) (* (- 3) x3) ) (- 484)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 501))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 501))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 501))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 501))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 501))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 501))