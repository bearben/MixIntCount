(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(assert (<= (+ (* 1 x1) (* 3 x2) (* 3 x3) (* (- 7) x4) ) (- 208)))
(assert (<= (+ (* 5 x1) (* 8 x2) (* (- 8) x3) (* 9 x4) ) 255))
(assert (<= (+ (* 2 x1) (* (- 6) x2) (* 10 x3) (* 5 x4) ) (- 156)))
(assert (<= (+ (* 6 x1) (* 1 x2) (* 3 x3) (* (- 10) x4) ) (- 63)))
(assert (<= (+ (* (- 5) x1) (* (- 7) x2) (* 8 x3) (* (- 3) x4) ) (- 137)))
(assert (<= (+ (* (- 3) x1) (* (- 3) x2) (* 8 x3) (* 2 x4) ) 56))
(assert (<= (+ (* (- 9) x1) (* (- 9) x2) (* 6 x3) (* 3 x4) ) (- 72)))
(assert (<= (+ (* 3 x1) (* (- 3) x2) (* 8 x3) (* 5 x4) ) 20))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 281))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) ) 281))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) ) 281))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) ) 281))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) ) 281))