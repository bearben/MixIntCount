(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(assert (<= (+ (* 7 x1) (* (- 3) x2) (* 10 x3) ) 262))
(assert (<= (+ (* 6 x1) (* 4 x2) (* 0 x3) ) (- 43)))
(assert (<= (+ (* 4 x1) (* 0 x2) (* (- 2) x3) ) (- 35)))
(assert (<= (+ (* 7 x1) (* 2 x2) (* 0 x3) ) 149))
(assert (<= (+ (* 1 x1) (* 0 x2) (* (- 4) x3) ) (- 51)))
(assert (<= (+ (* (- 5) x1) (* 7 x2) (* (- 8) x3) ) (- 89)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) ) 281))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) ) 281))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) ) 281))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) ) 281))
