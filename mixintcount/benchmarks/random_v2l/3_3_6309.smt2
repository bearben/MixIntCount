(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 10 x1) (* (- 8) x2) (* 1 x3) ) (- 4481)))
(assert (<= (+ (* 10 x1) (* 7 x2) (* (- 8) x3) ) 3830))
(assert (<= (+ (* 6 x1) (* 4 x2) (* 0 x3) ) 439))
(assert (<= (+ (* 5 x1) (* (- 2) x2) (* (- 7) x3) ) 1097))
(assert (<= (+ (* 7 x1) (* 4 x2) (* 10 x3) ) (- 2373)))
(assert (<= (+ (* 4 x1) (* 2 x2) (* 8 x3) ) 4441))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 6309))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 6309))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 6309))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 6309))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 6309))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 6309))