(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(declare-fun x5 () Int)
(assert (<= (+ (* 6 x1) (* 1 x2) (* (- 1) x3) (* (- 10) x4) (* 3 x5) ) 20))
(assert (<= (+ (* 6 x1) (* (- 7) x2) (* 4 x3) (* 2 x4) (* 9 x5) ) 111))
(assert (<= (+ (* 2 x1) (* (- 5) x2) (* 1 x3) (* 5 x4) (* (- 9) x5) ) 91))
(assert (<= (+ (* (- 10) x1) (* (- 4) x2) (* 7 x3) (* (- 10) x4) (* 10 x5) ) 41))
(assert (<= (+ (* 6 x1) (* 4 x2) (* 7 x3) (* (- 10) x4) (* (- 8) x5) ) (- 7)))
(assert (<= (+ (* 8 x1) (* 8 x2) (* 10 x3) (* (- 6) x4) (* 6 x5) ) 82))
(assert (<= (+ (* (- 6) x1) (* 1 x2) (* 6 x3) (* 9 x4) (* (- 9) x5) ) (- 32)))
(assert (<= (+ (* 6 x1) (* 1 x2) (* 7 x3) (* (- 4) x4) (* (- 5) x5) ) 14))
(assert (<= (+ (* 6 x1) (* (- 6) x2) (* (- 9) x3) (* 1 x4) (* (- 6) x5) ) 109))
(assert (<= (+ (* 3 x1) (* 9 x2) (* 2 x3) (* (- 1) x4) (* 7 x5) ) 10))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) (* 0 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 1 x5) ) 125))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* (- 1) x5) ) 125))