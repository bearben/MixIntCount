(set-logic QF_LIA)
(set-info :ran/ran2.smt2)
(set-info :smt-lib-version 2.0)
(declare-fun x1 () Int)
(declare-fun x2 () Int)
(declare-fun x3 () Int)
(declare-fun x4 () Int)
(declare-fun x5 () Int)
(assert (<= (+ (* 4 x1) (* 5 x2) (* 7 x3) (* (- 9) x4) (* 3 x5) ) (- 39)))
(assert (<= (+ (* (- 10) x1) (* 8 x2) (* (- 5) x3) (* 6 x4) (* (- 6) x5) ) (- 39)))
(assert (<= (+ (* (- 9) x1) (* (- 10) x2) (* 10 x3) (* (- 5) x4) (* 10 x5) ) 6))
(assert (<= (+ (* (- 5) x1) (* (- 10) x2) (* 7 x3) (* (- 5) x4) (* 8 x5) ) 25))
(assert (<= (+ (* (- 3) x1) (* 1 x2) (* 10 x3) (* 5 x4) (* 1 x5) ) 34))
(assert (<= (+ (* (- 4) x1) (* (- 8) x2) (* (- 10) x3) (* (- 10) x4) (* (- 8) x5) ) (- 57)))
(assert (<= (+ (* (- 9) x1) (* (- 8) x2) (* 1 x3) (* (- 6) x4) (* 8 x5) ) 18))
(assert (<= (+ (* (- 8) x1) (* 9 x2) (* 1 x3) (* 10 x4) (* (- 8) x5) ) 21))
(assert (<= (+ (* (- 3) x1) (* (- 3) x2) (* (- 3) x3) (* (- 9) x4) (* 0 x5) ) 58))
(assert (<= (+ (* (- 5) x1) (* 5 x2) (* 6 x3) (* (- 6) x4) (* (- 3) x5) ) (- 7)))
(assert (<= (+ (* 1 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* (- 1) x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 1 x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* (- 1) x2) (* 0 x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 1 x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* (- 1) x3) (* 0 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 1 x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* (- 1) x4) (* 0 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* 1 x5) ) 63))
(assert (<= (+ (* 0 x1) (* 0 x2) (* 0 x3) (* 0 x4) (* (- 1) x5) ) 63))
