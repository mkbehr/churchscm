;;; Discrete distributions

;; Simulates a fair coin flip: returns 1 with probability 1/2 and 0
;; with probability 1/2.
(define flip-spec
  (make-pspec
   (lambda () (random 2))
   (lambda ()
     (lambda (x)
       (mass->logmass
        (cond ((= x 1) 1/2)
              ((= x 0) 1/2)
              (else 0)))))))
(define flip (pspec->operator flip-spec))

;; Bernoulli distribution: returns 1 with probability p and 0 with
;; probability 1-p.
(define bernoulli-spec
  (make-pspec
   (lambda (p)
     (if (<= (random 1.0) p)
         1
         0))
   (lambda (p)
     (lambda (x)
       (mass->logmass
        (cond ((= x 1) p)
              ((= x 0) (- 1 p))
              (else 0)))))))
(define bernoulli
  (pspec->operator bernoulli-spec))

;;; "Continuous" distributions

;; These distributions are in fact implemented as discrete
;; distributions, based on the machine's precision.


(define (logdensity->logmass logdensity x)
  (let ((epsilon
         ;; This calculation adapted from gjs's scmutils: find the
         ;; smallest number that makes a difference when added to x
         (if (= x (+ 1. x))
             (let loop ((e 1.0))
               (if (= x (+ e x))
                   (loop (* e 2))
                   e))
             (let loop ((e 1.0))
               (if (= x (+ (/ e 2) x))
                   e
                   (loop (/ e 2)))))))
    (+ logdensity (log epsilon))))
(define density->logdensity mass->logmass)

;; Uniform distribution: returns a value between a inclusive and b
;; exclusive. Assumes a < b.
(define cont-uniform-spec
  (make-pspec
   (lambda (a b)
     (+ (* (random 1.0)
           (- b a))
        a))
   (lambda (a b)
     (lambda (x)
       (if (and (>= x a) (< x b))
           (logdensity->logmass
            (- (log (- b a)))
            x)
           -inf)))))
(define cont-uniform
  (pspec->operator cont-uniform-spec))

;; Normal distribution with mean of mean and standard deviation of
;; stdev.
(define normal-spec
  (make-pspec
   (lambda (mean stdev)
     ;; draw from standard normal distribution using marsaglia polar
     ;; method
     ;; note: if it is desired to draw as few random numbers as
     ;; possible, this can be rewritten to generate two standard
     ;; normal deviates instead of one, the other being
     ;; (* v (sqrt (/ (* -2 (log s)) s)))
     (let lp ((u (- 1 (random 2.0)))
              (v (- 1 (random 2.0))))
       (let ((s (+ (square u) (square v))))
         (if (>= s 1.)
             (lp (- 1 (random 2.0)) (- 1 (random 2.0)))
             (let ((z (* u (sqrt (/ (* -2 (log s)) s)))))
               (+ mean (* z stdev)))))))
   (lambda (mean stdev)
     (lambda (x)
       (logdensity->logmass
        (density->logdensity
         (/
          (exp (-
                (/ (square (- x mean))
                   (* 2 (square stdev)))))
          (* stdev (sqrt (* 2 *pi*)))))
        x)))))
(define normal
  (pspec->operator normal-spec))